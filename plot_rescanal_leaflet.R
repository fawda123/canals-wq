library(sf)
library(leaflet)
library(dplyr)

# import residential canal polygons with distance zonal stats
canals <- st_read(
  'data/ResCanal_Polygons_10212026_WMDs_with_DistZonalStats.shp',
  quiet = TRUE
)

# fix the one invalid geometry and project to WGS84 for leaflet
canals <- canals |>
  st_make_valid() |>
  st_transform(crs = 4326)

# IWR water quality data and stations inside canals, both from extract_iwr_wq.R
iwr <- readRDS('data/iwrwq.rds')
canalsta <- readRDS('data/canalsta.rds') |>
  st_transform(crs = 4326)

# summarize water quality data by station and parameter -------------------

canaliwr <- iwr |>
  filter(sta %in% canalsta$sta) |>
  mutate(
    date = as.Date(paste(year, month, day, sep = '-')),
    result = suppressWarnings(as.numeric(result))
  )

# one row per station and parameter, ordered so the latest result is last
parsta <- canaliwr |>
  arrange(sta, masterCode, date) |>
  group_by(sta, masterCode) |>
  summarise(
    n = n(),
    datemin = min(date, na.rm = TRUE),
    datemax = max(date, na.rm = TRUE),
    min = min(result, na.rm = TRUE),
    med = median(result, na.rm = TRUE),
    max = max(result, na.rm = TRUE),
    last = last(result),
    .groups = 'drop'
  )

# html table of parameters for each station popup
fmt <- function(x) formatC(x, digits = 3, format = 'g')
partab <- parsta |>
  mutate(
    row = paste0(
      '<tr><td>', masterCode, '</td><td>', n, '</td><td>',
      format(datemin, '%Y-%m-%d'), ' to ', format(datemax, '%Y-%m-%d'), '</td><td>',
      fmt(min), '</td><td>', fmt(med), '</td><td>', fmt(max), '</td><td>', fmt(last), '</td></tr>'
    )
  ) |>
  group_by(sta) |>
  summarise(tab = paste(row, collapse = ''))

# overall date range and sample count for each station
stasum <- canaliwr |>
  group_by(sta) |>
  summarise(
    nobs = n(),
    datemin = min(date, na.rm = TRUE),
    datemax = max(date, na.rm = TRUE)
  )

# station points with popup text
stapts <- canalsta |>
  left_join(distinct(iwr, sta, staname), by = 'sta') |>
  left_join(stasum, by = 'sta') |>
  left_join(partab, by = 'sta') |>
  mutate(
    pop = paste0(
      '<b>Station:</b> ', sta, '<br>',
      '<b>Name:</b> ', staname, '<br>',
      '<b>WBID:</b> ', wbid, '<br>',
      '<b>Canal ID:</b> ', CHAWQ_ID, '<br>',
      '<b>Dates:</b> ', format(datemin, '%Y-%m-%d'), ' to ', format(datemax, '%Y-%m-%d'), '<br>',
      '<b>Samples:</b> ', nobs, '<br><br>',
      '<table style="font-size:11px; border-collapse:collapse;" border="1" cellpadding="2">',
      '<tr><th>Parameter</th><th>n</th><th>Dates</th><th>Min</th><th>Median</th><th>Max</th><th>Last</th></tr>',
      tab, '</table>'
    )
  )

# map ---------------------------------------------------------------------

# color polygons by mean distance
pal <- colorNumeric(
  palette = 'viridis',
  domain = canals$DIST_MEAN,
  na.color = 'grey'
)

# canal popup text
pop <- paste0(
  '<b>Canal ID:</b> ', canals$CHAWQ_ID, '<br>',
  '<b>Waterbody:</b> ', canals$WATERBODY_, '<br>',
  '<b>County:</b> ', canals$County, '<br>',
  '<b>WMD:</b> ', canals$WMDNAME, '<br>',
  '<b>Receiving water:</b> ', canals$TYP_RECWAT, '<br>',
  '<b>Form:</b> ', canals$TYP_FORM, '<br>',
  '<b>Connection:</b> ', canals$TYP_CONNEC, '<br>',
  '<b>Area (acres):</b> ', round(canals$area_acres, 2), '<br>',
  '<b>Mean distance:</b> ', round(canals$DIST_MEAN, 1), '<br>',
  '<b>Median distance:</b> ', round(canals$DIST_MED, 1), '<br>',
  '<b>Max distance:</b> ', round(canals$DIST_MAX, 1)
)

m <- leaflet() |>
  addProviderTiles(providers$Esri.WorldGrayCanvas) |>
  addPolygons(
    data = canals,
    fillColor = ~pal(DIST_MEAN),
    fillOpacity = 0.8,
    color = 'black',
    weight = 0.5,
    popup = pop,
    label = ~WATERBODY_,
    highlightOptions = highlightOptions(weight = 2, color = 'red', bringToFront = FALSE),
    group = 'Canals'
  ) |>
  addCircleMarkers(
    data = stapts,
    radius = 5,
    color = 'black',
    weight = 1,
    fillColor = 'orange',
    fillOpacity = 0.9,
    popup = ~pop,
    popupOptions = popupOptions(maxWidth = 600),
    label = ~sta,
    group = 'IWR stations'
  ) |>
  addLegend(
    data = canals,
    pal = pal,
    values = ~DIST_MEAN,
    title = 'Mean distance',
    position = 'bottomright'
  ) |>
  addLayersControl(
    overlayGroups = c('Canals', 'IWR stations'),
    options = layersControlOptions(collapsed = FALSE)
  )

m
