

# Created 22.07.2020


# Population weights of the exposure 


###########################################################################################

library(raster)
library(sf)
library(FNN)
library(dplyr)
library(viridis)
library(ggplot2)
library(spdep)
library(lwgeom)


# Read gridded World population files. They can be downloaded from https://sedac.ciesin.columbia.edu/data/collection/gpw-v4
r_1 = raster("data/gpw_v4_population_count_rev11_2015_30_sec_2.asc")
r_2 = raster("data/gpw_v4_population_count_rev11_2015_30_sec_3.asc")

# read the land of England file
england.shape = readRDS("data/EnglandLand")
england.shape = spTransform(england.shape, crs(r_2))

t_0 = Sys.time()
s_1 = crop(r_1, england.shape)
s_2 = crop(r_2, england.shape)
t_1 = Sys.time()
t_1 - t_0 # 40 secs

plot(s_1)
plot(s_2)

# and combine rasters
s_3 <- mosaic(s_1, s_2, fun = mean)
plot(s_3)
plot(england.shape, add = TRUE)

# keep only England
s_4 <- mask(s_3, england.shape)
plot(s_4)

# covert it to polygon
s_5 <- rasterToPolygons(s_4)
plot(s_5[1:10000,])
# this conversion is slightly problematic in the boundaries.

# read LSOA file (only the boundaries)
LSOA <- readRDS("data/LSOA")

s_5 <- spTransform(s_5, crs(LSOA))
s_6 <- st_as_sf(s_5)
s_6$ID <- 1:nrow(s_6)

# make sure the obejcts are valid
LSOA <- LSOA %>% st_set_precision(1000000) %>% lwgeom::st_make_valid()
s_6 <- s_6 %>% st_set_precision(1000000) %>% lwgeom::st_make_valid()

# compute the intersection
t_0 = Sys.time()
s_7 <- st_intersection(s_6, LSOA) 
t_1 = Sys.time()
t_1 - t_0 # 2.3mins
s_7


sum(is.na(s_7$layer))
plot(s_7$geometry[1:1000])

sum(st_area(LSOA))
sum(st_area(s_7))


length(unique(s_7$lsoa11cd)) - nrow(LSOA)
# there are 8 LSOAs missing from the intersection

length(unique(s_7$ID)) - nrow(s_6)
# but all the grid cells are there

# lets check which LSOAs are missing
england.shape <- spTransform(england.shape, crs(LSOA))

plot(LSOA$geometry[!(LSOA$lsoa11cd %in% s_7$lsoa11cd)], border = "transparent", 
     col = "red")
plot(england.shape, add = TRUE)

par(ask = TRUE)
for(i in 1:8)plot(LSOA$geometry[!(LSOA$lsoa11cd %in% s_7$lsoa11cd)][i], border = "transparent", 
                  col = "red")

# add them
missingLSOAs <- LSOA$geometry[!(LSOA$lsoa11cd %in% s_7$lsoa11cd)]
missingLSOAs <- st_as_sf(missingLSOAs)
missingLSOAs$lsoa11cd <- LSOA$lsoa11cd[!(LSOA$lsoa11cd %in% s_7$lsoa11cd)]
# and I will need to assign them ID grid cell (the closest) and the corresponding population

# first union

missingLSOAs$layer = NA
missingLSOAs$ID = NA
missingLSOAs <- as(missingLSOAs, "Spatial")
missingLSOAs@data <- missingLSOAs@data[,c("layer", "ID", "lsoa11cd")]
missingLSOAs <- st_as_sf(missingLSOAs)

# 
t_0 = Sys.time()
s_8 <- rbind(s_7, missingLSOAs)
t_1 = Sys.time()
t_1 - t_0 # 18s


length(unique(s_8$lsoa11cd)) - nrow(LSOA)
# perfect match

length(unique(s_8$ID)) - nrow(s_6)
# the extra 1 is the NA

# Now I need to impute the NA in the ID and layer. I will do by using the NN
s_9 <- as(s_8, "Spatial")

na.points <- as.matrix(coordinates(s_9)[is.na(s_9$ID),])
notna.points <- as.matrix(coordinates(s_9)[!is.na(s_9$ID),])

NN <- get.knnx(notna.points, na.points, k = 1) $nn.index

s_9$layer[is.na(s_9$ID)] <- s_9$layer[NN]
s_9$ID[is.na(s_9$ID)] <- s_9$ID[NN]


s_10 <- st_as_sf(s_9)

length(unique(s_10$lsoa11cd)) - nrow(LSOA)
# perfect match
length(unique(s_10$ID)) - nrow(s_6)
# perfect match

# saveRDS(s_10, file = "data/PopEngland1kmLSOA")

# load the file with the population of LSOAs and put together with the World population
LSOApop <- readRDS("data/LSOApop")
LSOApop <- as.data.frame(LSOApop)
LSOApop$geometry <- NULL
colnames(LSOApop)[2] <- "LSOApop"
s_11 <- left_join(s_10, LSOApop, by = c("lsoa11cd" = "lsoa11cd")) 
colnames(s_11)[1:2] <- c("grid1kmpop", "ID1km")

# check the sums if they are close
sum(LSOApop$LSOApop) - # 55,977,178
  sum(s_5$layer) # 56,540,128

# ok the discrepancy is not big, plus have in mind that the world pop is for 2020, which an increase makes sense


# And now implement the algorithm for informing the LSOA weights from the world population
s_11$area <- as.numeric(st_area(s_11))
s_11$w <- s_11$grid1kmpop/s_11$area

w_agg <- aggregate(s_11$w, by = list(ID1km = s_11$ID1km), sum)
head(w_agg)
max(w_agg$ID1km) == nrow(s_6) # correct
colnames(w_agg)[2] <- "w_agg"
summary(w_agg)
sum(is.na(w_agg$w_agg))

s_12 <- left_join(s_11, w_agg, by = c("ID1km" = "ID1km")) 
sum(is.na(s_12$w_agg))
s_12$w_tilde <- s_12$w/s_12$w_agg
sum(is.na(s_12$w_tilde))
# the NAs here are cells with no population
s_12$w_tilde[is.na(s_12$w_tilde)] = 0

sum(s_12$w_tilde) - nrow(s_6) # it is not identical since some grid cells are 0s. They equal to the zero grid cells
sum(s_6$layer == 0)


# get the population by ij based on the grid
s_12$popij <- s_12$grid1kmpop*s_12$w_tilde
sum(s_12$popij) == sum(s_5$layer) # perfect match!

# Now calculate population weights based on pij. The normalization should be based on LSOAs now.

P_agg <- aggregate(s_12$popij, by = list(lsoa11cd = s_12$lsoa11cd), sum)
length(P_agg$lsoa11cd) == nrow(LSOApop) # correct
colnames(P_agg)[2] <- "P_agg"
summary(P_agg)
sum(is.na(P_agg$P_agg))
sum(P_agg$P_agg) - sum(s_5$layer) # as the gridded population since we havent used LSOA population yet.

# merge back
s_13 <- left_join(s_12, P_agg, by = c("lsoa11cd" = "lsoa11cd"))
s_13$v <- s_13$popij/s_13$P_agg
sum(s_13$v) == nrow(LSOApop) # perfect!

s_13$Popijtilde <- s_13$v*s_13$LSOApop
sum(s_13$Popijtilde) - sum(LSOApop$LSOApop)

# and store
# saveRDS(s_13, file = "data/PopEngland1kmLSOA")


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################











