############################################
### Script for Microlophus Bergmann paper ##
############################################

############################################
### Setup ##################################
############################################
## set working environment to folder containing script 
## (if script it being run in R studio)
if(Sys.getenv("RSTUDIO") == "1")
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# define var ending function
var_endings <- function(vars, fcat){
  
  pasted <- c(rep("0", length(vars)))
  for(j in 1:length(vars)){
    pasted[j] <- paste(vars[j], fcat, sep="")
    
  }
  
  return(pasted)
}
#

#load packages#
library(Rcpp)
library(rgeos)
library(proj4)
library(rgdal)
library(dplyr)
library(raster)
library(sp)
library(WorldClimTiles)
library("sf")

########################################
### Load in climate data ###############
########################################

##make SA map with high res tiles##
#first, get country boundaries for countries that have microlophus
peru <- getData("GADM", country = "PER", level = 0)
chi <- getData("GADM", country = "CHILE", level = 0)
ecu <- getData("GADM", country = "ECUADOR", level = 0)
col <- getData("GADM", country = "COLOMBIA", level = 0)
arg <- getData("GADM", country = "ARGENTINA", level = 0)

#get tiles associated with country boundaries
tile_nms <- unique(c(tile_name(peru),tile_name(chi),tile_name(ecu),tile_name(col),tile_name(arg)))

# get world clim data
tiles <- tile_get(tile_nms, "bio")

#get altitudinal data
tile_alt <- tile_get(tile_nms, "alt")

# merge tiles
alt_merge <- tile_merge(tile_alt)
merged <- tile_merge(tiles)

#stack combined raster of altitude and worldclim
comb_rast <- stack(merged, alt_merge)

#######################################
### Get species' range data ###########
#######################################

# NOTE: obtain these shape files from the IUCN Red List website and add them to the 'shape_files' folder
# https://www.iucnredlist.org/

####create vector of shape locations####
file_shx <- c("./shape_files/Microlophus_albemarlensis.shx","./shape_files/Microlophus_atacamensis.shx",
              "./shape_files/Microlophus_bivittatus.shx","./shape_files/Microlophus_delanonis.shx",
              "./shape_files/Microlophus_duncanensis.shx",
              "./shape_files/Microlophus_grayii.shx","./shape_files/Microlophus_habelii.shx",
              "./shape_files/Microlophus_indefatigabilis.shx","./shape_files/Microlophus_occipitalis.shx",
              "./shape_files/Microlophus_pacificus.shx","./shape_files/Microlophus_peruvianus.shx",
              "./shape_files/Microlophus_stolzmanni.shx","./shape_files/Microlophus_theresiae.shx",
              "./shape_files/Microlophus_theresioides.shx","./shape_files/Microlophus_thoracicus.shx",
              "./shape_files/Microlophus_tigris.shx","./shape_files/Microlophus_yanezi.shx",
              "./shape_files/Microlophus_tarapacensis.shx","./shape_files/Microlophus_quadrivittatus.shx",
              "./shape_files/Microlophus_koepckeorum.shx","./shape_files/Microlophus_indefatigabilis.shx")

# initialize data holders for succeeding loop
shp_list <- list()
rast_data <- list()
shape_order_spec <- c(rep("0", length(file_shx)))
ext_list <- data.frame(matrix(nrow=21,ncol=4))
colnames(ext_list) <- c("long_min","long_max","lat_min","lat_max")

# also want range centroids for lat/lon correlation
centroids <- data.frame(matrix(nrow=21,ncol=2))
colnames(centroids) <- c("long_cent","lat_cent")

# generate data list items for each shape
for (i in 1:length(file_shx)){
  
  #extract species name from shape file name, pass to vector
  shape_order_spec[i] <- substr(file_shx[i], 15, (nchar(file_shx[i])-4))
  
  # read shape, check projection, make a spatial object
  shp <- st_read(file_shx[i])
  shp <- st_transform(shp, crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  shp <- as_Spatial(shp)
  
  # add shape to list, then extract data from high res raster layers
  shp_list[[i]] <- shp
  rast_data[[i]] <- extract(comb_rast,shp)
  
  # get spatial extent of shape (min/max lat and lon)
  ext_list[i,] <- as.vector(extent(shp))
  centroids[i,] <- c(extent(gCentroid(shp))[1],extent(gCentroid(shp))[3])
}

##########################################################
### summarize species data, turn into a dataframe ########
##########################################################

# names of all the bioclim variables
coln <- c("annual_mean_temp", "mean_diurnal_range","isothermality","temperature_seasonality", 
          "max_temp_warmest_month", "min_temp_coldest_month", 
          "ann_temperature_range", "mean_temp_wettest_quarter", 
          "mean_temp_driest_quarter", "mean_temp_warmest_quarter", 
          "mean_temp_coldest_quarter", "ann_precipitaiton", 
          "precip_wettest_month", "precip_driest_month",
          "precip_seasonality", "precip_wettest_quarter", 
          "precip_driest_quarter", "precip_warmest_quarter", "precip_coldest_quarter", "altitude")

# generate empty dataframe to hold summarized data
rast_sum <-data.frame(matrix(nrow=length(file_shx), ncol=121))

# run through each shape's list item, make it a df, extract summary stats
for(i in 1:length(rast_data)){
  
  dat <- rast_data[[i]]
  # check if shape file contained >1 polygon (max in data set was 2)
  if(length(dat)==2)
  {
    dat <- data.frame(rbind(dat[[1]],dat[[2]]))
    
  }else{
    dat <- data.frame(dat)
  }
  
  # get summary stats
  gr_cel <- length(dat[,1]) # number of grid cells
  means <- as.numeric(lapply(data.frame(dat),mean,na.rm=TRUE)) # means (for all 19 bioclim variables)
  stdevs <- as.numeric(lapply(data.frame(dat),sd,na.rm=TRUE)) # standard deviations
  sterr <- as.numeric(stdevs)/sqrt(gr_cel) # standard errors
  medians <- as.numeric(lapply(data.frame(dat),median,na.rm=TRUE)) # medians
  mins <- as.numeric(lapply(data.frame(dat),min,na.rm=TRUE)) # minimum values
  maxs <- as.numeric(lapply(data.frame(dat),max,na.rm=TRUE)) # max values
  
  # put all summary stats into one large df
  rast_sum[i,] <- c(gr_cel,means,stdevs,sterr,medians,mins,maxs)
}

# generate column names using the var_endings function
cnames <- c("grid_cells", var_endings(coln, "_mean"), var_endings(coln, "_stdev"), var_endings(coln,"_sterr"), var_endings(coln, "_median"), var_endings(coln, "_min"), var_endings(coln,"_max"))

#apply column names
colnames(rast_sum) <- cnames

# add species names back in
rast_sum$species <- shape_order_spec

# combine all data
rast_sum <- cbind(rast_sum,ext_list,centroids)


######################################################################
### get climate data for species with points only ####################
######################################################################

## get data for barringtonensis and heterolepis (no shapes available)
sps <- c("Microlophus barringtonensis","Microlophus heterolepis")

# load points file
micro <- read.csv("./barringtonensis_and_heterolepis.csv", header=TRUE)

# grab just species of interest from points file
stitch <- rbind(micro[micro$species=="Microlophus barringtonensis",], micro[micro$species=="Microlophus heterolepis",])

# remove all duplicate points
op <- unique(stitch[c("species","longitude", "latitude", "loc")])

# remove non lat lon columns for data acquisition
opt2 <- op[,-c(1,4)]

# get bioclim and altitudinal data for each point
r_points <- data.frame(extract(comb_rast,opt2))

# add other variables back in
r_points$species <- op$species
r_points$lon <- op$longitude
r_points$lat <- op$latitude

# generate data holders
p_sum <-data.frame(matrix(nrow=length(sps), ncol=121))

p_ext <- data.frame(matrix(nrow=2, ncol=4))
colnames(p_ext) <- c("long_min","long_max","lat_min","lat_max")

p_cent <- data.frame(matrix(nrow=length(sps),ncol=2))
colnames(p_cent) <- c("long_cent","lat_cent")

# loop through both species
for(i in 1:length(sps)){
  
  # get points for species from df
  subs <- r_points[r_points$species==sps[i],]
  
  # extract max/min extent of range
  p_ext[i,] <- c(min(subs$lon),max(subs$lon),min(subs$lat),max(subs$lat))
  
  # average all points to get 'center' of distribution
  p_cent[i,] <- c(mean(subs$lon),mean(subs$lat))
  
  # remove lat, lon, species
  subs$species <- NULL
  subs$lon <- NULL
  subs$lat <- NULL
  
  # generate summary stats
  gr_cel <- length(subs[,1]) # number of points
  means <- as.numeric(lapply(data.frame(subs),mean,na.rm=TRUE))# means for those points (across all 19 variables)
  stdevs <- as.numeric(lapply(data.frame(subs),sd,na.rm=TRUE))# standard deviations
  sterr <- as.numeric(stdevs)/sqrt(gr_cel) # standard errors
  medians <- as.numeric(lapply(data.frame(subs),median,na.rm=TRUE))# medians
  mins <- as.numeric(lapply(data.frame(subs),min,na.rm=TRUE))# minimum values
  maxs <- as.numeric(lapply(data.frame(subs),max,na.rm=TRUE))# maximum values
  
  p_sum[i,] <- c(gr_cel,means,stdevs,sterr,medians,mins,maxs)# put all data into one large df
}

# adjust column names
colnames(p_sum) <- cnames

#add species names in
p_sum$species <- c("Microlophus_barringtonensis", "Microlophus_heterolepis")

#bind column with geographical data
p_sum <- cbind(p_sum,p_ext,p_cent)

## put it all together (shape and point data)
all <- rbind(rast_sum,p_sum)

# write to files
write.csv(all, "microlophus_climdat_all_hres.csv", row.names = FALSE)
saveRDS(all, "microlophus_climdat_all_hres.rds")

##############################################################################################
### get climate data for specific collection locations for species with large ranges #########
##############################################################################################

### get data set points for species that have large ranges ###
lr_points <- read.csv("large_range_species_points.csv")

# save species column, remove for lat/long data acquisition
spec_lr <- lr_points$species

# convert data to numeric
spec_ds <- lr_points[,-1]
spec_ds[,1] <- as.numeric(spec_ds[,1])
spec_ds[,2] <- as.numeric(spec_ds[,2])

# rename large range points columns
colnames(spec_ds) <- c("longitude","latitude")

# extract data from 0.5 degree worldclim data
lr <- data.frame(extract(comb_rast,spec_ds))
colnames(lr) <- coln

# add species names
lr$species <- spec_lr

write.csv(lr, "microlophus_climdat_large_range.csv", row.names = FALSE)
saveRDS(lr, "microlophus_climdat_large_range.rds")
##########################