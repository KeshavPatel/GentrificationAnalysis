findCenters = function(shapefile)
{
  tic = proc.time()
  cat("Finding Centers of Each Block Group. ")
  
  #initialize data frame that will contain the centroids
  centers = data.frame("GEOID"=shapefile$GEOID, "latitude"=0, "longitude"=0)
  for(i in 1:length(shapefile))
  {
    center_i = gCentroid(shapefile[i,]) #finds centroid
    pos_i = coordinates(center_i) #get latitude and longitude
    centers[shapefile[i,]$GEOID==shapefile$GEOID,2] = pos_i[2]
    centers[shapefile[i,]$GEOID==shapefile$GEOID,3] = pos_i[1]
    
  }
  toc = proc.time() - tic
  cat(paste0("done -- it took ", round(toc[3]), " seconds.\n"))
  return(centers)
}

getBlockGroup = function(state, county, endYear, yearSpan)
{
  tic = proc.time()
  cat("Getting US Census Geography and Population Data. ")
  #api.key.install(key="169f8fc7980039c827d1568b5217da952c27a11e") #old R package
  census_api_key(key="169f8fc7980039c827d1568b5217da952c27a11e")
  
  #format city and county inputs into objects the shapefile package recognizes
  cityCodes=geo.lookup(state=state, county=county)
  stateCodeStr = str_pad(cityCodes[2,1],2,pad="0")
  countyCodeStr = str_pad(cityCodes[2,3],3,pad="0")
  
  shapefile = block_groups(state=stateCodeStr, county=countyCodeStr)
  #mytableNew = acs.lookup(endyear=endYear, table.number="B01003") #old R package
  #myvarsNew = mytableNew[1] #old R package
  #mygeo <- geo.make(state=as.numeric(stateCodeStr), county=as.numeric(countyCodeStr), tract="*", block.group="*") #old R package
  #mydataNew <- acs.fetch(endyear=endYear, span=yearSpan, geography=mygeo, variable=myvarsNew) #old R package
  
  #get census data
  tarr <- get_acs(geography = "block group", variables = "B01003_001", 
                  state=state, county=county, geometry = TRUE)
  
  #acsgeoidNew <- paste0(as.character(mydataNew@geography$state),
  #                      str_pad(as.character(mydataNew@geography$county), 3, pad='0'),
  #                      str_pad(as.character(mydataNew@geography$tract), 6, pad='0'), as.character(mydataNew@geography$blockgroup)) #old R package
  
  #format data and combine with shapefile
  mydatadfNew <- data.frame(tarr$GEOID, tarr$estimate)
  colnames(mydatadfNew)=c("GEOID", "Population")
  shapefile = geo_join(shapefile, mydatadfNew, "GEOID", "GEOID")
  toc = proc.time() - tic
  cat(paste0("done -- it took ", round(toc[3]), " seconds.\n"))
  return(shapefile)
}

getWealthData = function(state, county, endYear, yearSpan, measurement)
{
  tic = proc.time()
  cat("Getting US Census Wealth Data. ")
  #cityCodes=geo.lookup(state=state, county=county) #old R package
  #stateCodeStr = str_pad(cityCodes[2,1],2,pad="0") #old R package
  #countyCodeStr = str_pad(cityCodes[2,3],3,pad="0") #old R package
  #mytableNew = acs.lookup(endyear=endYear, table.number="B25077") #old R package
  #myvarsNew = mytableNew[1] #old R package
  #mygeo <- geo.make(state=as.numeric(stateCodeStr), county=as.numeric(countyCodeStr), tract="*", block.group="*") #old R package
  #mydataNew <- acs.fetch(endyear=endYear, span=yearSpan, geography=mygeo, variable=myvarsNew) #old R package
  
  #get census data
  tarr <- get_acs(geography = "block group", variables = "B25077_001", 
                  state=state, county=county, geometry = TRUE)
  
  #acsgeoidNew <- paste0(as.character(mydataNew@geography$state),
  #                      str_pad(as.character(mydataNew@geography$county), 3, pad='0'),
  #                      str_pad(as.character(mydataNew@geography$tract), 6, pad='0'), as.character(mydataNew@geography$blockgroup)) #old R package
  
  #format data
  mydatadfNew <- data.frame(tarr$GEOID, tarr$estimate)
  colnames(mydatadfNew)=c("GEOID", measurement)
  toc = proc.time() - tic
  cat(paste0("done -- it took ", round(toc[3]), " seconds.\n"))
  return(mydatadfNew)
}

toAmenDensityPolyFrame = function(amen_df, blockGroupShapes)
{
  tic = proc.time()
  cat(paste0("\rPlacing Amenities on the map"))
  
  #initialize dataframe that holds the number of amenities in each block group
  amenDensity = data.frame("GEOID"=blockGroupShapes$GEOID, "amenities"=0)
  for(amen in 1:nrow(amen_df))
  {
    #get longitude and latitude and store it as a SpatialPoint for use in the next line  
    sp = SpatialPoints(coords=data.frame("longitude"=amen_df[amen,]$longitude, "latitude"=amen_df[amen,]$latitude), proj4string=CRS(getCRS(blockGroupShapes[1,])))
    shape = over(sp, blockGroupShapes) #determines what block group the SpatialPoint is in
    if(!is.na(shape[1,1]))
    {
      amenDensity[amenDensity$GEOID==shape$GEOID,2] = amenDensity[amenDensity$GEOID==shape$GEOID,2] + 1
    }
    cat(paste0("\rPlacing Amenities on the map: ", floor(amen/nrow(amen_df)*100), "% done. Expect this to take ", round((proc.time()[3]-tic[3])/amen*(nrow(amen_df)-amen)), " more seconds"))
  }
  toc = proc.time() - tic
  cat(paste0("\rPlacing Amenities on the map. done -- it took ", round(toc[3]), " seconds.                    \n"))
  return(amenDensity)
}

makePositionMap = function(positions, initLat, finLat, initLng, finLng, circleColor, addLegend, df)
{
  mypopupA <- paste0("Lat: ", positions$latitude, "<br>", "Lng: ", positions$longitude)
  mymapA = leaflet(options = leafletOptions(zoomControl = FALSE)) %>% 
    #addPolygons(data = df) %>%
    addCircles(positions$longitude, positions$latitude, weight = 5, radius=10, color=circleColor, stroke = TRUE, fillOpacity = 1, popup = mypopupA) #%>% 
  setView(map=mymapA,lng=-87.76, lat=41.935, zoom = 12.4)#fitBounds(initLng, initLat, finLng, finLat)
  if(addLegend)
  {
    mymapA %>% 
      addScaleBar(position="bottomright")
  }
  return(mymapA)
}

makeDensityMap = function(df, measurement, initLat, finLat, initLng, finLng, lowerBoundColor, upperBoundColor, addLegend)
{
  mypopupA <- paste0("GEOID: ", df$GEOID, "<br>", paste0(measurement, ": ", df[[measurement]]))
  mypalA = colorNumeric(
    palette = colorRamp(c(lowerBoundColor, upperBoundColor), interpolate="spline"),
    domain = df[[measurement]]
  )
  mymapA = leaflet(options = leafletOptions(zoomControl = FALSE)) %>%
    addProviderTiles("CartoDB.Positron") %>%
    addPolygons(data = df, 
                fillColor = ~mypalA(df[[measurement]]), 
                color = "#000000", # you need to use hex colors
                fillOpacity = 1, 
                weight = 1, 
                smoothFactor = 0.2,
                popup = mypopupA) %>%
      setView(-87.76, 41.935, zoom = 12.4)#fitBounds(initLng, initLat, finLng, finLat)
  if(addLegend)
  {
    mymapA %>% addLegend(pal = mypalA, 
                         values = df[[measurement]], 
                         position = "bottomright", 
                         title = measurement)  
    mymapA %>% addScaleBar(position="bottomleft")
  }
  return(mymapA)
}

yelpDataGET = function(amenList, initialLat, initialLong, vertLength, horizLength)
{
  tic = proc.time()
  cat(paste0("\rSearching for amenities"))
  ID = "Yhbe-GZ70MtSzRG64-QsXw"
  key = "ZR8GNlsaz_QU3V0L0pwSQjfoFW3EG5P2bWoyzJ6I9CoVsZOR9ntjaZscy_DngjkglFiHOnxNxpQsGiO0ccYaUPfnV1uOUDDmx9ONaGR-6ENMGf3kviPrOihZ8ykgW3Yx"
  
  #create emppty data frame with latitudes and longitudes of the amenities
  amen_df = data.frame("latitude"=c(), "longitude"=c())
  finalLat = initialLat
  finalLong = initialLong
  
  #for each category given, record the number of entries found (not including repeats)
  hitsVec = c()
  iterNum = 0
  iterTot = (horizLength+1)*(vertLength+1)*(length(amenList))
  for(type in amenList)
  {
    numHits = 0
    for(latIter in 0:vertLength)
    {
      for(longIter in 0:horizLength)
      {
        iterNum = iterNum+1
        #center of the search radius
        latitude  = initialLat  - (latIter / 4000) * (180 / pi);
        longitude = initialLong + (longIter / 4000) * (180 / pi) / cos(latitude * pi/180);
        url = paste0("https://api.yelp.com/v3/businesses/search?term=", type, "&limit=50&radius=1600&latitude=", latitude, "&longitude=", longitude) 
        off = 0
        numEntries = 0
        repeat
        {
          urlOff = paste0(url, "&offset=", off)
          foodJSON=GET(urlOff, add_headers("Authorization"= paste("Bearer", key))) #get Yelp data as a JSON
          contentList = content(foodJSON)$business #get the data frame out of the JSON
          
          if(off==0)
          {
            contentTotal = content(foodJSON)$total #the total number of entries available at this lat, long, and radius
            numEntries = min(contentTotal, 1000)
            if(numEntries == 0 | length(contentList) == 0) break
            numHits = numHits + numEntries
          }
          if(length(contentList) == 0) break
          for(bus in 1:length(contentList))
          {
            tempList = list("latitude"=contentList[[bus]]$coordinates$latitude, "longitude"=contentList[[bus]]$coordinates$longitude)
            amen_df = rbind(amen_df, tempList) #add the latitude and longitude to the data frame
          }
          off = off+50
          if(off > numEntries) break
        }
        cat(paste0("\rSearching for amenities: ", floor(100*iterNum/iterTot), "% done; Expect this to take ", round((proc.time()[3]-tic[3])/iterNum*(iterTot-iterNum)), ' more seconds'))
      }
    }
    hitsVec = c(hitsVec, numHits)
  }
  amen_df = unique(amen_df) #removes duplicates
  out = list(amen_df, latitude, longitude)
  toc = proc.time() - tic
  cat(paste0("\rSearching for amenities. done -- it took ", round(toc[3]), " seconds.                    \n"))
  for(amen in 1:length(amenList))
  {
    print(paste0("found ", hitsVec[amen], " locations for amenity: ", amenList[amen]))
  }
  return(out)
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x[x != 0], probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 3 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- (qnt[1] - H)
    y[x > (qnt[2] + H)] <- (qnt[2] + H)
    return(y)
}
