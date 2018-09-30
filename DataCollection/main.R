list.of.packages <- c("acs", "tigris", "leaflet", "mapview", "stringr", "plotKML", "sp", "rgeos", "httr", "purrr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(acs)
library(tigris)
library(leaflet)
library(mapview)
library(stringr)
library(plotKML)
library(sp)
library(rgeos)
library(httr)
library(purrr)

source("chicagoMedPropVal.R")

stateName = "IL" #2 letter state code
countyName = "Cook" #County name
wealthMeasurement = "MedianValue"
wealthCutoff = 225900*2 #Set to NA to remove wealth cutoff
myspan = 5 #number of years to average data over. Use 1, 3, or 5
myendyearNew = 2010 #Last year over which to average data
tractCodes= "*"#c(20000:89900, 110000:179900, 190000:259900, 830000:849900) #Set to "*" to include all tracts
initialLatitude = 41.5;
initialLongitude = -87.93063;
horizontalLength = 10; #miles
verticalLength = 10;
amenityList = c("food")#c("food", "parks", "coffee", "music", "nightlife")

#shapefile contains the map of the specified area
shapefile = getBlockGroup(stateName, countyName, myendyearNew, myspan)
if(class(tractCodes) == "integer")
{
  shapefile = shapefile[as.numeric(substring(shapefile$GEOID,6,11)) %in% tractCodes,]
}

centers = findCenters(shapefile)
df_temp1 = geo_join(shapefile, centers, "GEOID", "GEOID")

#wealthPolyFrame contains the US Census wealth data
wealthPolyFrame = getWealthData(stateName, countyName, myendyearNew, myspan, wealthMeasurement)
if(class(tractCodes) == "integer")
{
  wealthPolyFrame = wealthPolyFrame[as.numeric(substring(wealthPolyFrame$GEOID,6,11)) %in% tractCodes,]
}

df_temp2 = geo_join(df_temp1, wealthPolyFrame, "GEOID", "GEOID")
wealth_df = df_temp2[!is.na(df_temp2[[wealthMeasurement]]),]
if(!is.na(wealthCutoff))
{
  #wealth_df = wealth_df[wealth_df[[wealthMeasurement]]>wealthCutoff,]
  wealth_df[wealth_df[[wealthMeasurement]]<wealthCutoff,"latitude"] = NA
  wealth_df[wealth_df[[wealthMeasurement]]<wealthCutoff,"longitude"] = NA
}


returnList = yelpDataGET(amenityList, initialLatitude, initialLongitude, verticalLength, horizontalLength)
amenitiesPositions = returnList[[1]]
finalLatitude = returnList[[2]]
finalLongitude = returnList[[3]]
#shapefile_adj = shapefile[shapefile$Population>300,]
#amenitiesPolyFrame = toAmenDensityPolyFrame(amenitiesPositions, shapefile_adj)
#amenity_df = geo_join(shapefile_adj, amenitiesPolyFrame, "GEOID", "GEOID")
#amenity_df$amenities = amenity_df$amenities/amenity_df$Population
#amenity_df = amenity_df[!is.na(amenity_df$amenities),]

cat("Generating Maps. ")
wealthPositionMap = makePositionMap(wealth_df, initialLatitude, finalLatitude, initialLongitude, finalLongitude, "#0000FF", TRUE)
#amenityPositionMap = makePositionMap(amenitiesPositions, initialLatitude, finalLatitude, initialLongitude, finalLongitude, "#FF0000", TRUE)
wealthDensityMap = makeDensityMap(wealth_df, wealthMeasurement, initialLatitude, finalLatitude, initialLongitude, finalLongitude, "#000000", "#0000FF", TRUE)
#amenityDensityMap = makeDensityMap(amenity_df, "amenities", initialLatitude, finalLatitude, initialLongitude, finalLongitude, "#000000", "#FF0000", TRUE)
cat("done\n")

#mypalA = colorNumeric(
#  palette = colorRamp(c("#000000", "#FF0000"), interpolate="spline"),
#  domain = amenity_df$amenities
#)
#amenityDensityMap %>% addLegend(pal = mypalA, 
#                     values = amenity_df$amenities, 
#                     position = "bottomright", 
#                     title = "Amenities Per Person") %>% addScaleBar(position="bottomleft")

