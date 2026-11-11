# Multivariate grouping of canal water quality
#
# Canals that drain to lakes, reservoirs, rivers or creeks are dropped to
# focus on estuarine and marine canals. Each canal polygon with IWR stations is summarized as the median of all
# results for each parameter. The medians are transformed and scaled, then
# grouped with Ward's hierarchical clustering. The number of groups is chosen
# from average silhouette width, and a PCA shows how the groups separate.
#
# Requires data/iwrwq.rds and data/canalsta.rds from extract_iwr_wq.R.

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cluster)
library(ggrepel)

# settings ----------------------------------------------------------------

# parameters used for clustering. Only canals with all of them are kept.
# With data since 2015 and the receiving water filter, these five give 191
# canals. Adding TURB, COLOR or TSS drops that number considerably. DOSAT
# and COND are left out because they are nearly the same as DO and SALIN
# (r > 0.94).
wqvars <- c('SALIN', 'DO', 'TN', 'TP', 'CHLAC')

# canal attributes used for clustering, from the canal shapefile. DIST_MEAN
# is the mean distance within the canal, likely in meters (the source
# projection is in meters). DIST_MED and DIST_MAX are left out because they
# are nearly the same as DIST_MEAN (r > 0.97).
canalvars <- c('area_acres', 'DIST_MEAN')

vars <- c(wqvars, canalvars)

# parameters that are log10 transformed because they are right skewed
logvars <- c('TN', 'TP', 'CHLAC', 'TURB', 'NO3O2', 'TKN', 'COLOR', 'TSS', 'area_acres', 'DIST_MEAN')

# first year of data to include. NA uses the full record, which goes back to
# the 1930s for some canals. 2015 keeps roughly the last decade, which gives
# a current picture at the cost of fewer canals (191 vs 342 with all years).
yrmin <- 2015

# receiving water types (TYP_RECWAT) to drop
dropwat <- c('Lake/Reservoir', 'River/Creek')

# number of groups. NA picks the k with the highest average silhouette width.
k <- NA
kmax <- 8

# data --------------------------------------------------------------------

iwr <- readRDS('data/iwrwq.rds')
canalsta <- readRDS('data/canalsta.rds')

canals <- st_read(
  'data/ResCanal_Polygons_10212026_WMDs_with_DistZonalStats.shp',
  quiet = TRUE
) |>
  st_drop_geometry()

cat('Canals dropped by receiving water type\n')
print(table(canals$TYP_RECWAT[canals$TYP_RECWAT %in% dropwat]))

canals <- canals |>
  filter(!TYP_RECWAT %in% dropwat) |>
  select(CHAWQ_ID, WATERBODY_, County, WMDNAME, TYP_RECWAT, TYP_FORM, TYP_CONNEC, TYP_SCALE, area_acres, DIST_MEAN)

# results at stations inside canals, tagged with the canal ID
canalwq <- iwr |>
  filter(masterCode %in% wqvars) |>
  inner_join(st_drop_geometry(canalsta)[, c('sta', 'CHAWQ_ID')], by = 'sta') |>
  filter(CHAWQ_ID %in% canals$CHAWQ_ID)
if (!is.na(yrmin))
  canalwq <- filter(canalwq, year >= yrmin)

# median of each parameter in each canal, one row per canal
canalmed <- canalwq |>
  group_by(CHAWQ_ID, masterCode) |>
  summarise(med = median(result, na.rm = TRUE), .groups = 'drop') |>
  pivot_wider(names_from = masterCode, values_from = med) |>
  select(CHAWQ_ID, all_of(wqvars)) |>
  left_join(canals[, c('CHAWQ_ID', canalvars)], by = 'CHAWQ_ID')

cat('Canals with any of the parameters:', nrow(canalmed), '\n')
canalmed <- canalmed[complete.cases(canalmed), ]
cat('Canals with all of the parameters:', nrow(canalmed), '\n')

# transform and scale -----------------------------------------------------

trans <- canalmed
for (v in intersect(vars, logvars)) {
  # half the smallest positive value is added so zeros can be logged
  offset <- min(trans[[v]][trans[[v]] > 0]) / 2
  trans[[v]] <- log10(trans[[v]] + offset)
}
scl <- scale(as.matrix(trans[, vars]))
rownames(scl) <- canalmed$CHAWQ_ID

# clustering --------------------------------------------------------------

dis <- dist(scl)
hc <- hclust(dis, method = 'ward.D2')

# average silhouette width for each k
silk <- tibble(k = 2:kmax) |>
  rowwise() |>
  mutate(silwidth = mean(silhouette(cutree(hc, k), dis)[, 'sil_width'])) |>
  ungroup()
print(silk)

if (is.na(k))
  k <- silk$k[which.max(silk$silwidth)]
cat('Number of groups:', k, '\n')

grp <- cutree(hc, k)

# groups are numbered by increasing median salinity so the numbers mean
# something and stay stable if k changes
ord <- tapply(canalmed$SALIN, grp, median) |> rank(ties.method = 'first')
grp <- factor(ord[as.character(grp)], levels = 1:k)

canalclust <- canalmed |>
  mutate(grp = grp) |>
  left_join(select(canals, -all_of(canalvars)), by = 'CHAWQ_ID')

# group summaries ---------------------------------------------------------

cat('\nCanals in each group\n')
print(table(canalclust$grp))

# median of the canal medians in each group, in original units
cat('\nGroup profiles (median of canal medians)\n')
grpprof <- canalclust |>
  group_by(grp) |>
  summarise(n = n(), across(all_of(vars), median))
print(grpprof, width = Inf)

# canal attributes by group
for (v in c('WMDNAME', 'County', 'TYP_RECWAT', 'TYP_FORM', 'TYP_CONNEC', 'TYP_SCALE')) {
  cat('\n', v, 'by group\n')
  print(table(canalclust[[v]], canalclust$grp, useNA = 'ifany'))
}

# pca ---------------------------------------------------------------------

pca <- prcomp(scl)
pcvar <- round(100 * summary(pca)$importance['Proportion of Variance', ], 1)
cat('\nPCA variance explained (%)\n')
print(pcvar)
cat('\nPCA loadings\n')
print(round(pca$rotation, 2))

# plots -------------------------------------------------------------------

# fixed categorical order. Groups also get a shape, because past three
# colors a scatter can't rely on color alone for color-blind readers.
cols <- c('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948')[1:k]
shps <- c(16, 17, 15, 18, 3, 4, 8, 1)[1:k]

# silhouette width by k
p1 <- ggplot(silk, aes(k, silwidth)) +
  geom_line(linewidth = 0.6, colour = 'grey40') +
  geom_point(size = 2.5, colour = 'grey20') +
  geom_point(data = filter(silk, k == !!k), size = 4, colour = cols[1]) +
  scale_x_continuous(breaks = silk$k) +
  labs(x = 'Number of groups', y = 'Average silhouette width') +
  theme_minimal()

# pca biplot with groups
scrs <- as_tibble(pca$x[, 1:2]) |> mutate(grp = grp)
lds <- as_tibble(pca$rotation[, 1:2], rownames = 'var')
lscl <- 0.8 * max(abs(scrs[, c('PC1', 'PC2')])) / max(abs(lds[, c('PC1', 'PC2')]))
p2 <- ggplot(scrs, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = 'grey85') +
  geom_vline(xintercept = 0, colour = 'grey85') +
  geom_point(aes(colour = grp, shape = grp), size = 2, alpha = 0.7) +
  geom_segment(
    data = lds, aes(x = 0, y = 0, xend = PC1 * lscl, yend = PC2 * lscl),
    arrow = arrow(length = unit(0.2, 'cm')), colour = 'grey20'
  ) +
  geom_label_repel(
    data = lds, aes(x = PC1 * lscl * 1.1, y = PC2 * lscl * 1.1, label = var),
    size = 3.2, label.size = 0, fill = 'white', alpha = 0.8, min.segment.length = Inf
  ) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shps) +
  labs(
    x = paste0('PC1 (', pcvar[1], '%)'), y = paste0('PC2 (', pcvar[2], '%)'),
    colour = 'Group', shape = 'Group'
  ) +
  coord_equal() +
  theme_minimal()

# canal medians by group, one panel per parameter. Logged parameters are
# plotted as log10 values so skewed panels stay readable.
p3dat <- canalclust |>
  select(CHAWQ_ID, grp, all_of(vars)) |>
  pivot_longer(all_of(vars), names_to = 'var', values_to = 'val') |>
  filter(!(var %in% logvars & val <= 0)) |>
  mutate(
    val = ifelse(var %in% logvars, log10(val), val),
    var = ifelse(var %in% logvars, paste0('log10 ', var), var),
    var = factor(var, levels = ifelse(vars %in% logvars, paste0('log10 ', vars), vars))
  )
p3 <- ggplot(p3dat, aes(grp, val)) +
  geom_boxplot(outlier.shape = NA, colour = 'grey30', width = 0.6) +
  geom_jitter(aes(colour = grp), width = 0.15, size = 1, alpha = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = cols) +
  facet_wrap(~var, scales = 'free_y') +
  labs(x = 'Group', y = 'Canal median (water quality) or value (area, distance)') +
  theme_minimal()

print(p1)
print(p2)
print(p3)

# save --------------------------------------------------------------------

wqclust <- list(
  settings = list(wqvars = wqvars, canalvars = canalvars, logvars = intersect(vars, logvars), yrmin = yrmin, dropwat = dropwat, k = k),
  canalclust = canalclust,
  grpprof = grpprof,
  silk = silk,
  hc = hc,
  pca = pca
)
saveRDS(wqclust, 'data/wqclust.rds')
