library(sf)
library(leaflet)
library(dplyr)

# leaflet map of canal polygons and IWR stations --------------------------
#
# canals: canal polygons (sf) with the shapefile attributes
# stations: IWR stations (sf) with sta, wbid and CHAWQ_ID
# wq: IWR results with sta, staname, masterCode, units, year, month, day and
#   result, used for the station popups
# colvar: canal column used for the polygon fill
# pal: leaflet palette function for colvar. NULL uses viridis on a numeric
#   column.
# legtitle: legend title
# popextra: optional html added to the end of each canal popup, in the same
#   order as canals
canalmap <- function(canals, stations, wq, colvar = 'DIST_MEAN', pal = NULL,
                     legtitle = 'Mean distance', popextra = NULL) {

  # project to WGS84 for leaflet
  canals <- st_transform(canals, crs = 4326)
  stations <- st_transform(stations, crs = 4326)

  wq <- wq |>
    filter(sta %in% stations$sta) |>
    mutate(date = as.Date(paste(year, month, day, sep = '-')))

  # one row per station and parameter, ordered so the latest result is last
  parsta <- wq |>
    arrange(sta, masterCode, date) |>
    group_by(sta, masterCode, units) |>
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
        '<tr><td>', masterCode, '</td><td>', units, '</td><td>', n, '</td><td>',
        format(datemin, '%Y-%m-%d'), ' to ', format(datemax, '%Y-%m-%d'), '</td><td>',
        fmt(min), '</td><td>', fmt(med), '</td><td>', fmt(max), '</td><td>', fmt(last), '</td></tr>'
      )
    ) |>
    group_by(sta) |>
    summarise(tab = paste(row, collapse = ''))

  # overall date range and sample count for each station
  stasum <- wq |>
    group_by(sta) |>
    summarise(
      nobs = n(),
      datemin = min(date, na.rm = TRUE),
      datemax = max(date, na.rm = TRUE)
    )

  # station points with popup text
  stapts <- stations |>
    left_join(distinct(wq, sta, staname), by = 'sta') |>
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
        '<tr><th>Parameter</th><th>Units</th><th>n</th><th>Dates</th><th>Min</th><th>Median</th><th>Max</th><th>Last</th></tr>',
        tab, '</table>'
      )
    )

  # polygon fill
  if (is.null(pal))
    pal <- colorNumeric(palette = 'viridis', domain = canals[[colvar]], na.color = 'grey')
  canals$fillval <- canals[[colvar]]

  # canal popup text
  pop <- paste0(
    '<b>Canal ID:</b> ', canals$CHAWQ_ID, '<br>',
    '<b>Waterbody:</b> ', canals$WATERBODY_, '<br>',
    '<b>County:</b> ', canals$County, '<br>',
    '<b>WMD:</b> ', canals$WMDNAME, '<br>',
    '<b>Receiving water:</b> ', canals$TYP_RECWAT, '<br>',
    '<b>Form:</b> ', canals$TYP_FORM, '<br>',
    '<b>Connection:</b> ', canals$TYP_CONNEC, '<br>',
    '<b>Scale:</b> ', canals$TYP_SCALE, '<br>',
    '<b>Area (acres):</b> ', round(canals$area_acres, 2), '<br>',
    '<b>Mean distance:</b> ', round(canals$DIST_MEAN, 1), '<br>',
    '<b>Median distance:</b> ', round(canals$DIST_MED, 1), '<br>',
    '<b>Max distance:</b> ', round(canals$DIST_MAX, 1)
  )
  if (!is.null(popextra))
    pop <- paste0(pop, popextra)

  # base maps for the layers control. The first one added is shown by default.
  basemaps <- c('Esri gray', 'Esri imagery')

  leaflet() |>
    addProviderTiles(providers$Esri.WorldGrayCanvas, group = 'Esri gray') |>
    addProviderTiles(providers$Esri.WorldImagery, group = 'Esri imagery') |>
    addPolygons(
      data = canals,
      fillColor = ~pal(fillval),
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
      values = ~fillval,
      title = legtitle,
      position = 'bottomright'
    ) |>
    addLayersControl(
      baseGroups = basemaps,
      overlayGroups = c('Canals', 'IWR stations'),
      options = layersControlOptions(collapsed = FALSE)
    )

}

# full map of all canals and stations -------------------------------------

# skipped when the script is sourced while rendering index.qmd, which only
# needs canalmap()
if (!isTRUE(getOption('knitr.in.progress'))) {

  # residential canal polygons, with the one invalid geometry fixed
  canals <- st_read(
    'data/ResCanal_Polygons_10212026_WMDs_with_DistZonalStats.shp',
    quiet = TRUE
  ) |>
    st_make_valid()

  # IWR water quality data and stations inside canals, both from extract_iwr_wq.R
  iwr <- readRDS('data/iwrwq.rds')
  canalsta <- readRDS('data/canalsta.rds')

  m <- canalmap(canals, canalsta, iwr)
  m

}
