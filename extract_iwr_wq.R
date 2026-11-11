# extract statewide water quality data from the FDEP IWR run database
# adapted from C:/proj/tidalcreek-stats/Creek_select_tbeptools.Rmd
# IWR runs are available here: https://publicfiles.dep.state.fl.us/dear/iwr/
# requires the 64-bit Microsoft Access ODBC driver

library(RODBC)
library(dplyr)
library(sf)

# zipped IWR database, update when the copy finishes
# should match the station run below (Run 68)
zippth <- 'C:/Users/mbeck/Desktop/iwr2026_Run68_2026-08-21.zip'

# extracted location, skips unzipping if the .accdb files are already there
pth <- tools::file_path_sans_ext(zippth)
if(length(list.files(pth, pattern = 'rawDataDB.*\\.accdb$', recursive = TRUE)) == 0)
  unzip(zippth, exdir = pth)

# main run database with the parameter lookup tables, unzip it if missing
if(length(list.files(pth, pattern = '^iwr.*\\.accdb$', recursive = TRUE)) == 0){
  mainfl <- grep('^iwr.*\\.accdb$', basename(unzip(zippth, list = TRUE)$Name), value = TRUE)
  unzip(zippth, files = mainfl, exdir = pth)
}
maindb <- list.files(pth, pattern = '^iwr.*\\.accdb$', full.names = TRUE, recursive = TRUE)
stopifnot(length(maindb) == 1)

# .accdb files to query, recursive in case the zip has a nested folder
dbs <- list.files(pth, pattern = 'rawDataDB.*\\.accdb$', full.names = TRUE, recursive = TRUE)
stopifnot(length(dbs) > 0)

# parameters to extract, all years
mcode <- c("CHLAC","COLOR","COND","DO","DOSAT","NO3O2","ORGN","SALIN","TKN","TN","TP","TSS","TURB") %>%
  paste0("masterCode = '", ., "'") %>%
  paste0(collapse = ' or ')

# query builder, table name in each .accdb matches the file name
qry_fun <- function(tab){

  out <- paste0("select sta, year, month, day, time, depth, dUnits, masterCode, param, result, rCode from ", tab,
                " where ", mcode)
  return(out)

}

# loop through the accdb files
iwr <- vector('list', length(dbs))
for(i in seq_along(dbs)){

  db <- dbs[i]
  cat(basename(db), '\n')

  # setup the connection
  con <- odbcConnectAccess2007(db)

  # build the query
  qry <- basename(db) %>%
    gsub('DB|\\.accdb', '', .) %>%
    qry_fun()

  # retrieve the data
  iwr[[i]] <- sqlQuery(con, qry, stringsAsFactors = FALSE)

  odbcClose(con)

}
iwr <- bind_rows(iwr)

# units for each master parameter code, from the main run database
# param is the STORET code for each result, IWR maps source units to these codes
con <- odbcConnectAccess2007(maindb)
units <- sqlFetch(con, 'Master parameter codes', stringsAsFactors = FALSE) %>%
  select(masterCode = `Master parameter Code`, units = Units) %>%
  mutate(units = enc2utf8(units)) %>%
  filter(masterCode %in% unique(iwr$masterCode))
odbcClose(con)
stopifnot(!any(duplicated(units$masterCode)))

iwr <- left_join(iwr, units, by = 'masterCode')

# depth in meters, dUnits is the unit for depth
# f is assumed to be feet, other codes such as 64 are set to NA
iwr <- iwr %>%
  mutate(
    depth_m = case_when(
      dUnits == 'm' ~ depth,
      dUnits %in% c('ft', 'f') ~ depth * 0.3048,
      dUnits == 'in' ~ depth * 0.0254,
      TRUE ~ NA_real_
    )
  )

# depth units and how many results have a depth
print(count(iwr, dUnits, hasdepth = !is.na(depth)))

# check that values are comparable across STORET codes within each parameter
unitchk <- iwr %>%
  group_by(masterCode, units, param) %>%
  summarise(
    n = n(),
    p10 = quantile(result, 0.1, na.rm = TRUE),
    med = median(result, na.rm = TRUE),
    p90 = quantile(result, 0.9, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(masterCode, desc(n))
print(as.data.frame(unitchk))

# station locations and wbids, extracted from Stations_Run68.zip
# the shapefile repeats some stations as exact duplicates, so keep distinct rows
stas <- st_read('data/Stations_Run68/IWR_stations_Run68.shp', quiet = TRUE) %>%
  st_drop_geometry() %>%
  select(sta = STATION_ID, wbid = WATERBODY_, staname = STATION_NA, lat = LATITUDE, lon = LONGITUDE) %>%
  distinct()
stopifnot(!any(duplicated(stas$sta)))

# add station info to the water quality data
iwr <- left_join(iwr, stas, by = 'sta')

# report stations in the water quality data without a location
nomatch <- unique(iwr$sta[is.na(iwr$lat)])
cat(length(nomatch), 'of', length(unique(iwr$sta)), 'stations have no location\n')

head(iwr)

# save
saveRDS(iwr, 'data/iwrwq.rds')

# intersect stations with residential canal polygons ----------------------

# canal polygons, fix the one invalid geometry
canals <- st_read('data/ResCanal_Polygons_10212026_WMDs_with_DistZonalStats.shp', quiet = TRUE) %>%
  st_make_valid()

# unique stations with water quality data and a location, projected to match canals
iwrsta <- iwr %>%
  filter(!is.na(lat)) %>%
  distinct(sta, wbid, lat, lon) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
  st_transform(crs = st_crs(canals))

# stations inside a canal polygon, with the canal id
canalsta <- st_join(iwrsta, canals['CHAWQ_ID'], join = st_intersects, left = FALSE)

cat(length(unique(canalsta$sta)), 'of', nrow(iwrsta), 'stations fall inside a canal polygon\n')
cat(length(unique(canalsta$CHAWQ_ID)), 'of', nrow(canals), 'canal polygons have at least one station\n')

# save canal stations for mapping
saveRDS(canalsta, 'data/canalsta.rds')

# breakdown by parameter --------------------------------------------------

# station to canal lookup
canalid <- canalsta %>%
  st_drop_geometry() %>%
  distinct(sta, CHAWQ_ID)

# statewide counts and canal station counts for each parameter, all years
parsum <- iwr %>%
  mutate(incanal = sta %in% canalid$sta) %>%
  group_by(masterCode) %>%
  summarise(
    nobs = n(),
    nsta = n_distinct(sta),
    yrmin = min(year, na.rm = TRUE),
    yrmax = max(year, na.rm = TRUE),
    nobs_canal = sum(incanal),
    nsta_canal = n_distinct(sta[incanal]),
    yrmin_canal = if(any(incanal)) min(year[incanal], na.rm = TRUE) else NA,
    yrmax_canal = if(any(incanal)) max(year[incanal], na.rm = TRUE) else NA
  )

# number of canal polygons with data for each parameter
parcanal <- iwr %>%
  distinct(masterCode, sta) %>%
  inner_join(canalid, by = 'sta') %>%
  group_by(masterCode) %>%
  summarise(ncanal = n_distinct(CHAWQ_ID))

parsum <- parsum %>%
  left_join(parcanal, by = 'masterCode') %>%
  mutate(ncanal = coalesce(ncanal, 0L)) %>%
  arrange(desc(ncanal))

print(parsum, n = Inf, width = Inf)
