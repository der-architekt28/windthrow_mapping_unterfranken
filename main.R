# load packages

library(lidR)
library(ITSMe)
library(terra)
library(fgr)
library(sf)
library(dplyr)
library(purrr)
library(data.table)
library(raster)
library(ggplot2)
library(scales)

####### 1. Preliminary Tasks #######

# load point cloud

pc <- readLAS("C:/Users/l-reg/Documents/Studium_Wuerzburg/LIDAR/data/500_5544_5545_merged.laz")

#plot(pc, color = "Intensity", bg = "black")

# assign UTM reference system

st_crs(pc) <- 25832

# ground classification

ws <- seq(3, 12, 3)
th <- seq(0.1, 1.0, length.out = length(ws))

pc_gclass <- classify_ground(pc, algorithm = pmf(ws = ws, th = th))



##### 1: DTM from pc #####

# create dtm from point cloud

dtm_tin <- rasterize_terrain(las = pc, res = 0.5, algorithm = tin())
plot(dtm_tin)

# DSM creation
dsm <- rasterize_canopy(las=pc, res=0.5, algorithm = p2r())
plot(dsm)

# CHM
chm <- dsm - dtm_tin
plot(chm)

#  smoothing

kernel <- matrix(1,5,5)
schm <- terra::focal(x = dsm, w = kernel, fun = median, na.rm = TRUE)
plot(schm)

# detect trees

ttops <- locate_trees(las = chm, algorithm = lmf(ws = 2.5))
ttops
plot(schm)
plot(ttops, cex=0.05, add = TRUE)


##### 2: external DEM #####

# load external DEM (Geodata Bayern)

dem <- rast("C:/Users/l-reg/Documents/Studium_Wuerzburg/LIDAR/data/kahlstein_dem.tif")

# normalize point cloud

pc_norm <- pc_gclass - dem
pc_above_ground <- filter_poi(pc_norm, Z > 0.4)


# create canopy height model

dsm2 <- rasterize_canopy(las=pc_above_ground, res=0.3, algorithm = p2r())

# dem_res <- resample(dem, dsm2, method="near")

# CHM
# chm2 <- dsm2 - dem_res
# plot(chm2)

#  apply smoothing

kernel <- matrix(1,5,5)
schm2 <- terra::focal(x = dsm2, w = kernel, fun = median, na.rm = TRUE)


# detect tree tops

ttops2 <- locate_trees(las = schm2, algorithm = lmf(ws = 6))


# Segment trees using the dalponte algorithm

seg2 <- segment_trees(las = pc_above_ground, algorithm = dalponte2016(chm = schm2, treetops = ttops2))

length(unique(seg2$treeID) |> na.omit())

# # Visualize using intensity values as colors
# plot(seg2, color = "treeID", bg = "white")




####### 2. Parameter Retrieval #######

# (A) average spacing between trees

# create buffer

b <- st_buffer(ttops2, 15)
plot(b$geometry, col = NA, border = "grey")
plot(ttops2, cex=0.05, add = TRUE)

b_rad <- 10
nn_val <- b_rad + 1


# preparing treetops dataframe

ttops_xy <- st_zm(ttops2, drop = TRUE, what = "ZM")
stopifnot(!anyDuplicated(ttops_xy$treeID))

# extract all neighbouring trees

next_nbs <- st_is_within_distance(
  ttops_xy,
  ttops_xy,
  dist = b_rad,
  sparse = TRUE
)

# remove reference tree from the list

next_nbs <- lapply(
  seq_along(next_nbs),
  function(i) setdiff(next_nbs[[i]], i)
)

# get coordinates

nb_xy <- st_coordinates(ttops_xy)[, 1:2, drop = FALSE]

# calculate mean distance to next neighbours

mean_nb_dist <- map_dbl(
  seq_along(next_nbs),
  function(i) {

    nb_id <- next_nbs[[i]]

    if (length(nb_id) == 0) {
      return(nn_val)
    }

    dx <- nb_xy[nb_id, 1] - nb_xy[i, 1]
    dy <- nb_xy[nb_id, 2] - nb_xy[i, 2]

    distances <- sqrt(dx^2 + dy^2)

    mean(distances)
  }
)

# add distances to the treetops dataframe

ttops2 <- ttops2 %>%
  mutate(
    n_neighbours = lengths(next_nbs),
    mean_neighbour_distance = mean_nb_dist,
    no_neighbours = n_neighbours == 0
  )


# (B) aerodynamic profile

# for each segment: length, width, number of points and point spacing

pts <- seg2@data


setDT(pts)

tree_stats <- pts[
  !is.na(treeID),
  .(
    n_pts = .N,

    X_min = min(X, na.rm = TRUE),
    X_max = max(X, na.rm = TRUE),
    Y_min = min(Y, na.rm = TRUE),
    Y_max = max(Y, na.rm = TRUE),
    Z_min = min(Z, na.rm = TRUE),
    Z_max = max(Z, na.rm = TRUE),

    # ground margin in Z direction is added again

    X_range = max(X, na.rm = TRUE) - min(X, na.rm = TRUE),
    Y_range = max(Y, na.rm = TRUE) - min(Y, na.rm = TRUE),
    Z_range = max(Z, na.rm = TRUE) - min(Z, na.rm = TRUE) + 0.4
  ),
  by = treeID
]

# match stats to dataframe

idx <- match(ttops2$treeID, tree_stats$treeID)

ttops2 <- cbind(
  ttops2,
  tree_stats[idx, !"treeID"]
)


# LAI estimation via point spatial density (Indirabai et al. 2020)

pt_density <- 25.07      # point cloud related parameter

d <- sqrt(1/pt_density)

# add spatial density

ttops2 <- ttops2 %>%
  mutate(
    pt_spacing = d
  )

# main LAI function

calc_LAI <- function(x, y, z, n, d) {
  ((n * d * (x + y)) / (x * y)) * (1 / z)
}

ttops2 <- ttops2 %>%
  mutate(
    LAI = pmap_dbl(
      list(
        X_range,
        Y_range,
        Z_range,
        n_pts,
        pt_spacing
      ),
      calc_LAI
    )
  )


# estimate zero plane displacement and aerodynamic roughness via Raupach model

# zero plane displacement

# calc_displ <- function(LAI, h) {
#
#   lambda <- 0.5 * LAI
#   displ <- h * (1 - ((1-exp(-sqrt(2*7.5*lambda))) / sqrt(2*7.5*lambda)))
#
#   return(displ)
# }
#
# ttops2 <- ttops2 %>%
#   mutate(
#     zero_plane_displ = pmap_dbl(
#       list(
#         LAI,
#         Z_range
#       ),
#       calc_displ
#     )
#   )
#
# # aerodynamic roughness z_0
#
# calc_z0 <- function(h, displ) {
#
#   z0 <- (h - displ) * exp(1)^(-0.4 * 0.31 + 0.193)
#
# }
#
# ttops2 <- ttops2 %>%
#   mutate(
#     aero_z0 = pmap_dbl(
#       list(
#         Z_range,
#         zero_plane_displ
#       ),
#       calc_z0
#     )
#   )


# zero plane displacement (Raupach 1994 / Floors et al. 2021)

calc_displ <- function(LAI, h) {
  lambda <- 0.5 * LAI
  a      <- sqrt(2 * 7.5 * lambda)
  displ  <- h * (1 - (1 - exp(-a)) / a)
  return(displ)
}


# aerodynamic roughness z_0 (Raupach 1994 / Floors et al. 2021)

calc_z0_raupach <- function(h, d, LAI,
                            kappa = 0.40,
                            CS    = 0.003,
                            CR    = 0.30,
                            cmax  = 0.30,
                            Psi_h = 0.193) {

  lambda <- 0.5 * LAI
  term   <- pmin(CS + CR * lambda, cmax)
  z0     <- (h - d) * exp(-kappa * term - Psi_h)
  return(z0)
}

# add both parameters to dataframe

ttops2 <- ttops2 %>%
  mutate(

    zero_plane_displ = calc_displ(LAI, Z_range),

    aero_z0 = calc_z0_raupach(
      h   = Z_range,
      d   = zero_plane_displ,
      LAI = LAI
    )
  )



# (C) diameter at breast height (DBH)

# # add XY coordinates to dataframe as separate columns
#
# ttop_coos <- st_coordinates(ttops_xy)
#
# ttops_xy <- ttops_xy %>%
#   mutate(
#     X_top  = ttop_coos[, 1],
#     Y_top  = ttop_coos[, 2]
#   )
#
# # define radius of virtual cylinder
#
# cyl_rad <- 1
#
# # matching point cloud with treetop dataframe
#
# pts <- seg2@data
# top_id <- match(pts$treeID, ttops_xy$treeID)
#
# # calculate distances
#
# cyl_dist <- (pts$X - ttops_xy$X_top[top_id])^2 +
#   (pts$Y - ttops_xy$Y_top[top_id])^2
#
# # creating virtual cylinders and filtering point cloud
#
# cyl_keep <- !is.na(top_id) & cyl_dist <= cyl_rad^2
# seg2_cylinder <- filter_poi(seg2, cyl_keep)
#
# # filter point cloud to approximate tree breast height
#
# seg2_bh <- filter_poi(seg2_cylinder, Z > 1 & Z < 5)
#
# seg2_ground <- filter_poi(seg2, Z > 1.5 & Z < 3)




# first step: construction of a virtual cylinder along the stem
# -> removal of most of the tree crown

# get the lowermost section of the trees

ttrunk <- filter_poi(
  seg2,
  Z > 0.4 & Z < 1.3
)

# calculate tree trunk centers

stem_pts <- as.data.table(ttrunk@data)

ttrunk_center <- stem_pts[
  !is.na(treeID),
  .(
    X_base = mean(X, na.rm = TRUE),
    Y_base = mean(Y, na.rm = TRUE),
    n_ttrunk = .N
  ),
  by = treeID
]

# get coordinates of the tree tops

top_coords <- st_coordinates(ttops2)

tops_df <- data.frame(
  treeID = ttops2$treeID,
  X_top  = top_coords[, 1],
  Y_top  = top_coords[, 2],
  Z_top  = top_coords[, 3]
)

# construct axis from trunk center to tree top

axis_df <- merge(
  as.data.frame(ttrunk_center),
  tops_df,
  by = "treeID"
)

stem_z <- stem_pts[
  !is.na(treeID),
  .(Z_base = mean(Z, na.rm = TRUE)),
  by = treeID
]

axis_df <- merge(
  axis_df,
  as.data.frame(stem_z),
  by = "treeID"
)

# define assumed maximum stem radius

proj_stem_rad <- 0.5

# match segmented trees with their respective stem axis

pts <- as.data.table(seg2@data)
id_axis <- match(pts$treeID, axis_df$treeID)

# vector trunk center - tree top

vx <- axis_df$X_top[id_axis] - axis_df$X_base[id_axis]
vy <- axis_df$Y_top[id_axis] - axis_df$Y_base[id_axis]
vz <- axis_df$Z_top[id_axis] - axis_df$Z_base[id_axis]

# trunk center - individual point

wx <- pts$X - axis_df$X_base[id_axis]
wy <- pts$Y - axis_df$Y_base[id_axis]
wz <- pts$Z - axis_df$Z_base[id_axis]

# axis length and projection to main axis

v2 <- vx^2 + vy^2 + vz^2
dot <- wx * vx + wy * vy + wz * vz
t <- dot / v2

# distance to axis

d2 <- wx^2 + wy^2 + wz^2 - dot^2 / v2
d2 <- pmax(d2, 0)

# keep every point within virtual cylinder

stem_cyl_pts <- !is.na(id_axis) &
  is.finite(t) &
  v2 > 0 &
  t >= 0 & t <= 1 &
  d2 <= proj_stem_rad^2

# filtering point cloud

pc_stem_cyl <- filter_poi(seg2, stem_cyl_pts)


# second task: clipping the stem segments in vertical sections
# goal: fit circles at different heights to obtain new tree-specific radius values

# fit horizontal circles

# filter cylinder pc from 0.5 to 2.0

pc_stem_05_20 <- filter_poi(pc_stem_cyl, Z > 0.5 & Z <= 2.0)

id_05_20 <- unique(pc_stem_05_20$treeID)
id_05_20 <- id_05_20[!is.na(id_05_20)]

circle_05_20 <- purrr::map_dfr(id_05_20, function(id) {

  tree_pc_05_20 <- filter_poi(pc_stem_05_20, treeID == id)

  fit <- tryCatch(
    fit_circle(tree_pc_05_20, num_iterations = 100, inlier_threshold = 0.5),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      treeID = id,
      radius_05_20 = NA_real_
    ))
  }

  tibble(
    treeID = id,
    radius_05_20 = fit$radius
  )
})

ttops2 <- ttops2 %>%
  left_join(circle_05_20, by = "treeID")

ttops2

# filter cylinder pc from 2.0 to 3.5

pc_stem_20_35 <- filter_poi(pc_stem_cyl, Z > 2.0 & Z <= 3.5)

id_20_35 <- unique(pc_stem_20_35$treeID)
id_20_35 <- id_20_35[!is.na(id_20_35)]

circle_20_35 <- purrr::map_dfr(id_20_35, function(id) {

  tree_pc_20_35 <- filter_poi(pc_stem_20_35, treeID == id)

  fit <- tryCatch(
    fit_circle(tree_pc_20_35, num_iterations = 100, inlier_threshold = 0.5),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      treeID = id,
      radius_20_35 = NA_real_
    ))
  }

  tibble(
    treeID = id,
    radius_20_35 = fit$radius
  )
})

ttops2 <- ttops2 %>%
  left_join(circle_20_35, by = "treeID")

ttops2

# filter cylinder pc from 3.5 to 5.0

pc_stem_35_50 <- filter_poi(pc_stem_cyl, Z > 3.5 & Z <= 5.0)

id_35_50 <- unique(pc_stem_35_50$treeID)
id_35_50 <- id_35_50[!is.na(id_35_50)]

circle_35_50 <- purrr::map_dfr(id_35_50, function(id) {

  tree_pc_35_50 <- filter_poi(pc_stem_35_50, treeID == id)

  fit <- tryCatch(
    fit_circle(tree_pc_35_50, num_iterations = 100, inlier_threshold = 0.5),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      treeID = id,
      radius_35_50 = NA_real_
    ))
  }

  tibble(
    treeID = id,
    radius_35_50 = fit$radius
  )
})

ttops2 <- ttops2 %>%
  left_join(circle_35_50, by = "treeID")

ttops2

# filter cylinder pc from 5.0 to 6.5

pc_stem_50_65 <- filter_poi(pc_stem_cyl, Z > 5.0 & Z <= 6.5)

id_50_65 <- unique(pc_stem_50_65$treeID)
id_50_65 <- id_50_65[!is.na(id_50_65)]

circle_50_65 <- purrr::map_dfr(id_50_65, function(id) {

  tree_pc_50_65 <- filter_poi(pc_stem_50_65, treeID == id)

  fit <- tryCatch(
    fit_circle(tree_pc_50_65, num_iterations = 100, inlier_threshold = 0.5),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      treeID = id,
      radius_50_65 = NA_real_
    ))
  }

  tibble(
    treeID = id,
    radius_50_65 = fit$radius
  )
})

ttops2 <- ttops2 %>%
  left_join(circle_50_65, by = "treeID")

ttops2


# calculate mean radius from the respective stem sections

calc_new_rad <- function(r1, r2, r3, r4) {
  mean(c(r1, r2, r3, r4), na.rm = TRUE)
}

ttops2 <- ttops2 %>%
  mutate(
    new_radius = pmap_dbl(
      list(
        radius_05_20,
        radius_20_35,
        radius_35_50,
        radius_50_65
      ),
      calc_new_rad
    )
  )

ttops2

# remove every new radius greater than 3

ttops2 <- ttops2 %>%
  mutate(
    new_radius = if_else(new_radius > 3, NA_real_, new_radius)
  )


# cylindric filtering: second iteration with tree-specific stem radius

new_stem_rad_df <- ttops2

if (inherits(new_stem_rad_df, "sf")) {
  new_stem_rad_df <- sf::st_drop_geometry(new_stem_rad_df)
}

# only take relevant columns
# for entries with unlikely radius, the previous value of 0.5 is used

new_stem_rad_df <- new_stem_rad_df |>

  dplyr::select(treeID, new_radius) |>
  dplyr::mutate(radius = dplyr::coalesce(new_radius, 0.5)) |>
  dplyr::select(treeID, radius)

# add radius to axis table

axis_df <- axis_df |>
  dplyr::select(-dplyr::any_of("radius")) |>
  dplyr::left_join(new_stem_rad_df, by = "treeID")

# second iteration of the filtering method

pts <- as.data.table(seg2@data)
id_axis <- match(pts$treeID, axis_df$treeID)

# vector trunk center - tree top

vx <- axis_df$X_top[id_axis] - axis_df$X_base[id_axis]
vy <- axis_df$Y_top[id_axis] - axis_df$Y_base[id_axis]
vz <- axis_df$Z_top[id_axis] - axis_df$Z_base[id_axis]

# trunk center - individual point

wx <- pts$X - axis_df$X_base[id_axis]
wy <- pts$Y - axis_df$Y_base[id_axis]
wz <- pts$Z - axis_df$Z_base[id_axis]

# axis length and projection to main axis

v2 <- vx^2 + vy^2 + vz^2
dot <- wx * vx + wy * vy + wz * vz
t <- dot / v2

# distance to axis

d2 <- wx^2 + wy^2 + wz^2 - dot^2 / v2
d2 <- pmax(d2, 0)

# keep every point within virtual cylinder

est_stem_rad <- axis_df$radius[id_axis]

new_stem_cyl_pts <- !is.na(id_axis) &
  is.finite(t) &
  is.finite(d2) &
  is.finite(est_stem_rad) &
  v2 > 0 &
  t >= 0 & t <= 1 &
  d2 <= est_stem_rad^2

# filtering point cloud

new_pc_stem_cyl <- filter_poi(seg2, new_stem_cyl_pts)

# define tree breast height
# setting of a larger section due to low point density

breast_height <- 1.3
bh_low <- breast_height - 0.5
bh_high <- breast_height + 0.5

# filter cylinder point cloud for breast height

pc_stem_bh <- filter_poi(new_pc_stem_cyl, Z > bh_low & Z < bh_high)

id_bh <- unique(pc_stem_bh$treeID)
id_bh <- id_bh[!is.na(id_bh)]

# fit a circle on the section

circle_bh <- purrr::map_dfr(id_bh, function(id) {

  tree_pc_bh <- filter_poi(pc_stem_bh, treeID == id)

  fit <- tryCatch(
    fit_circle(tree_pc_bh, num_iterations = 120, inlier_threshold = 0.5),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      treeID = id,
      dbh = NA_real_
    ))
  }

  # multiply the radius with 2 to obtain the diameter

  tibble(
    treeID = id,
    dbh = fit$radius * 2
  )
})

# add to dataframe

ttops2 <- ttops2 %>%
  left_join(circle_bh, by = "treeID")

ttops2



# circular fitting not succesfull for every tree segment
# --> additional use of alternative model

# DBH estimation by empirical model (Bhebhe et al. 2025)

# delete every dbh greater than 5

ttops2 <- ttops2 %>%
  mutate(
    dbh = if_else(dbh > 5, NA_real_, dbh)
  )


# define dbh function

calc_dbh <- function(H, a1, a2) {
  a1 * (H - 1.3)^a2
}

# fit only on valid data

fit_dbh <- ttops2 %>%
  filter(!is.na(dbh), Z_range > 1.3)

fit_dbh_fun <- nls(
  dbh ~ calc_dbh(Z_range, a1, a2),
  data = fit_dbh,
  start = list(a1 = 1, a2 = 0.3)
)

# replace NA-values

ttops2 <- ttops2 %>%
  mutate(
    dbh = if_else(
      is.na(dbh) & Z_range > 1.3,
      predict(fit_dbh_fun, newdata = cur_data()),
      dbh
    )
  )




# (D) biomass (stem weight) estimation (Bhebhe et_al. 2025)

calc_biomass <- function(DBH, H) {

  wood_dens <- (0.6 + 0.9) / 2
  biomass <- exp(-2.187 + 0.916 * log(wood_dens * ((DBH*100)^2) * H))

  return(biomass)
}

ttops2 <- ttops2 %>%
  mutate(
    AGB = pmap_dbl(
      list(
        dbh,
        Z_range
      ),
      calc_biomass
    )
  )



####### 3. Final CWS Calculation (Gardiner, 2004) #######

calc_CWS_break <- function(D, dbh, d, h, z0 ) {

  k = 0.4
  MOR = 105000
  p = 1.2
  G = 2.0
  f_knot = 0.9
  f_CW = 1

  CWS_break = (1/(k*D)) * ((pi*MOR*dbh^3)/( 32*p*G * (d-1.3)))^(1/2) * (f_knot/f_CW)^(1/2) * log((h-d)/ z0)

  return(CWS_break)
}

ttops2 <- ttops2 %>%
  mutate(
    CWS_break = pmap_dbl(
      list(
        mean_neighbour_distance,
        dbh,
        zero_plane_displ,
        Z_range,
        aero_z0
      ),
      calc_CWS_break
    )
  )

ttops2



calc_CWS_overturn <- function(D, SW, d, h, z0 ) {

  k = 0.4
  C_reg = 190
  p = 1.2
  G = 2.0
  f_CW = 1

  CWS_overturn = (1/(k*D)) * ((C_reg * SW)/( p*G *d))^(1/2) * (1/f_CW)^(1/2) * log((h-d)/ z0)

  return(CWS_overturn)

}

ttops2 <- ttops2 %>%
  mutate(
    CWS_overturn = pmap_dbl(
      list(
        mean_neighbour_distance,
        AGB,
        zero_plane_displ,
        Z_range,
        aero_z0
      ),
      calc_CWS_overturn
    )
  )



####### 4. Data Plot #######

# retrieve crown metrics for each tree

tcrowns <- lidR::crown_metrics(
  seg2,
  func = NULL,
  geom = "convex",
  attribute = "treeID"
)

# get tree attributes from dataframe

tree_attributes <- ttops2 |>
  sf::st_drop_geometry() |>
  dplyr::select(
    treeID,
    dplyr::any_of(c("mean_neighbour_distance",
                    "Z_range",
                    "LAI",
                    "zero_plane_displ",
                    "aero_z0",
                    "dbh",
                    "AGB",
                    "CWS_break",
                    "CWS_overturn"))
  ) |>
  dplyr::distinct(treeID, .keep_all = TRUE)

# attach tree attributes to crown shapes

tcrowns <- tcrowns |>
  dplyr::left_join(tree_attributes, by = "treeID")


# remove unplausible stem weight values

tcrowns <- tcrowns %>%
  mutate(
    AGB = ifelse(AGB > 30000, NA, AGB)
  )


# plots of the critical wind speeds

plot_cws_fun <- function(data, attribute,
                         title = NULL,
                         outline = "grey30",
                         strength = NULL) {

  stopifnot(attribute %in% names(data))

  # prepare values, check for validity

  x <- data[[attribute]]
  x_valid <- x[is.finite(x)]

  if (length(x_valid) == 0) {
    stop("Error: No valid values.")
  }

  middle <- median(x_valid)

  if (is.null(strength)) {
    strength <- diff(range(x_valid)) / 20
  }

  strength <- max(strength, .Machine$double.eps)

  # define median-centered transformation

  med_trans <- scales::trans_new(
    name = "median_centered",
    transform = function(x) {
      sign(x - middle) * log1p(abs(x - middle) / strength)
    },
    inverse = function(x) {
      middle + sign(x) * strength * expm1(abs(x))
    }
  )

  # calculate transformed values and median of the data

  tx <- med_trans$transform(x_valid)
  tx_range <- range(tx)

  mid_pos <- scales::rescale(0, from = tx_range)

  # set positions of the color transitions

  values <- c(
    0,
    mid_pos * 0.35,
    mid_pos * 0.70,
    mid_pos,
    mid_pos + (1 - mid_pos) * 0.35,
    mid_pos + (1 - mid_pos) * 0.70,
    1
  )

  # final plot

  ggplot(data) +
    geom_sf(
      aes(fill = .data[[attribute]]),
      colour = outline,
      linewidth = 0.15
    ) +
    scale_x_continuous(
      breaks = seq(9.007, 9.013, by = 0.002)
    ) +
    scale_y_continuous(
      breaks = seq(50.056, 50.062, by = 0.002)
    ) +
    scale_fill_gradientn(
      colours = c(
        "red",
        "orangered",
        "darkorange",
        "orange",
        "gold",
        "khaki",
        "lightyellow"
      ),
      values = values,
      na.value = "grey85",
      name = NULL,
      trans = med_trans,
      breaks = scales::breaks_pretty(n = 8),
      labels = scales::label_number()
    ) +
    guides(
      fill = guide_colorbar(
        barheight = grid::unit(8, "cm"),
        nbin = 256
      )
    ) +
    labs(title = title) +
    coord_sf() +
    theme_minimal()
}

# plots for CWS_break and CWS_overturn

cws_break_plot <- plot_cws_fun(
  tcrowns,
  attribute = "CWS_break",
  title = "CWS for breaking [m/s-1]",
  strength = 3
)

cws_overturn_plot <- plot_cws_fun(
  tcrowns,
  attribute = "CWS_overturn",
  title = "CWS for overturning [m/s-1]",
  strength = 3
)


# plots of the parameters used for calculation

# plot function

plot_params_fun <- function(data, attribute,
                            title = NULL,
                            outline = "grey30",
                            strength = NULL) {

  stopifnot(attribute %in% names(data))

  # prepare values, check for validity

  x <- data[[attribute]]
  x_valid <- x[is.finite(x)]

  if (length(x_valid) == 0) {
    stop("Error: No valid values.")
  }

  middle <- median(x_valid)

  if (is.null(strength)) {
    strength <- diff(range(x_valid)) / 20
  }

  strength <- max(strength, .Machine$double.eps)

  # define median-centered transformation

  med_trans <- scales::trans_new(
    name = "median_centered",
    transform = function(x) {
      sign(x - middle) * log1p(abs(x - middle) / strength)
    },
    inverse = function(x) {
      middle + sign(x) * strength * expm1(abs(x))
    }
  )

  # calculate transformed values and median of the data

  tx <- med_trans$transform(x_valid)
  tx_range <- range(tx)

  mid_pos <- scales::rescale(0, from = tx_range)

  # set positions of the color transitions

  values <- c(
    0,
    mid_pos * 0.35,
    mid_pos * 0.70,
    mid_pos,
    mid_pos + (1 - mid_pos) * 0.35,
    mid_pos + (1 - mid_pos) * 0.70,
    1
  )

  # final plot

  ggplot(data) +
    geom_sf(
      aes(fill = .data[[attribute]]),
      colour = outline,
      linewidth = 0.15
    ) +
    scale_x_continuous(
      breaks = seq(9.007, 9.013, by = 0.002)
    ) +
    scale_y_continuous(
      breaks = seq(50.056, 50.062, by = 0.002)
    ) +
    scale_fill_viridis_c(
      values = values,
      na.value = "grey85",
      name = NULL,
      trans = med_trans,
      breaks = scales::breaks_pretty(n = 8),
      labels = scales::label_number()
    ) +
    guides(
      fill = guide_colorbar(
        barheight = grid::unit(8, "cm"),
        nbin = 256
      )
    ) +
    labs(title = title) +
    coord_sf() +
    theme_minimal()
}


# parameter plots

dbh_plot <- plot_params_fun(
  tcrowns,
  attribute = "dbh",
  title = "Diameter at Breast Height (DBH) [m]",
  strength = 3
)

sw_plot <- plot_params_fun(
  tcrowns,
  attribute = "AGB",
  title = "Stem Weight (SW) [kg]",
  strength = 3
)

lai_plot <- plot_params_fun(
  tcrowns,
  attribute = "LAI",
  title = "Leaf Area Index (LAI)",
  strength = 3
)

spacing_plot <- plot_params_fun(
  tcrowns,
  attribute = "mean_neighbour_distance",
  title = "Mean Distance to Neighbours [m]",
  strength = 3
)

height_plot <- plot_params_fun(
  tcrowns,
  attribute = "Z_range",
  title = "Tree Height [m]",
  strength = 3
)

disp_plot <- plot_params_fun(
  tcrowns,
  attribute = "zero_plane_displ",
  title = "Zero Plane Displacement (d) [m]",
  strength = 3
)

z0_plot <- plot_params_fun(
  tcrowns,
  attribute = "aero_z0",
  title = "Aerodynamic Rougness (z0) [m]",
  strength = 3
)
