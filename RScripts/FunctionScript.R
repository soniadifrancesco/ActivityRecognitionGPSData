# Function script

getDataset = function(directory){
  print("Data acquisition...")
  # Set directory
  setwd(directory)
  # Variables names
  variablesName <- c("Latitude", "Longitude", "Altitude", "Bearing", "Accuracy", "Speed", "TimeStamp")
  # Files Reading
  participants = list.files()
  df=NULL
  for(p in participants){
    pDirectory = paste0(directory,"\\",p)
    setwd(pDirectory)
    filesList = list.files()
    sessions=NULL; i = 1
    for(f in filesList){
      sessionData <- read.csv(file = f , header = FALSE)
      sessionData <- sessionData[-(1:4),]
      if(nrow(sessionData)>0){
        sessionData[12] <- NULL
        sessionData[11] <- NULL
        sessionData[10] <- NULL
        sessionData[2] <- NULL
        sessionData[1] <- NULL
        names(sessionData) <- variablesName
        sessionData$TimeStamp <- as.POSIXct(strptime(sessionData$TimeStamp, "%Y-%m-%dT%H:%M:%S")) 
        sessionData = sessionData[order(sessionData$TimeStamp),]
        sessionData = unique(sessionData)
        sessionData$sessionid = i
        sessions <- rbind(sessions,sessionData)
        i = i+1
      }
    }
    sessions$patient = p
    df = rbind(df,sessions)
  }
  
  # Set type of variables
  df$patient = as.factor(df$patient)
  df$sessionid = as.factor(df$sessionid)
  df$Latitude = as.numeric(as.character(df$Latitude))
  df$Longitude = as.numeric(as.character(df$Longitude))
  df$Altitude = as.numeric(as.character(df$Altitude))
  df$Bearing = as.numeric(as.character(df$Bearing))
  df$Speed = as.numeric(as.character(df$Speed))
  df$Accuracy = as.numeric(as.character(df$Accuracy))
  df$X = c(1:nrow(df))
  print("Data acquisition completed.")
  return (df)
}

getDatasetGPX = function(directory){
  # Set directory
  setwd(directory)
  # Files Reading
  partecipants <- list.files()
  df=NULL
  for(p in partecipants){
    pDirectory = paste0(directory,"\\",p)
    setwd(pDirectory)
    filesList = list.files()
    sessions=NULL; i = 1
    print(p)
    for(f in filesList){
      print(f)
      gpx = readGPX(f)
      tracks = gpx$tracks
      Longitude = tracks[[1]][[1]]$lon
      Latitude = tracks[[1]][[1]]$lat
      TimeStamp = tracks[[1]][[1]]$time
      sessionData = data.frame(Latitude, Longitude, TimeStamp)
      sessionData$TimeStamp = as.POSIXct(strptime(sessionData$TimeStamp, "%Y-%m-%dT%H:%M:%SZ"))
      sessionData$sessionid = i
      sessions = rbind(sessions,sessionData)
      i=i+1
    }
    sessions$patient = p
    df = rbind(df,sessions)
  }
  df$X = c(1:nrow(df))
  return (df)
}

dateBasedCleaning = function(df, sfd){
  datePatients = unique(sfd[,c("patient","date")])
  df$date = as.Date(df$TimeStamp)
  output = NULL
  for(i in 1:nrow(datePatients)){
    tmp = NULL
    tmp = df[df$patient==datePatients$patient[i] & df$date==datePatients$date[i],]
    output = rbind(output,tmp)
  }
  output$date = NULL
  return(output)
}

dateBasedCleaningSF = function(df, sfd){
  df$date = as.Date(df$TimeStamp)
  datePatients = unique(df[,c("patient","date")])
  output = NULL
  for(i in 1:nrow(datePatients)){
    tmp = NULL
    tmp = sfd[sfd$patient==datePatients$patient[i] & sfd$date==datePatients$date[i],]
    output = rbind(output,tmp)
  }
  return(output)
}

timeBasedMethod <- function(df, timeThreshold, distanceThreshold){
  print("Geolocations visited detection in progress...")
  df$geoVisited = FALSE
  
  # temporary data frame with information:
  # 1) time difference between consecutive GPS points (i.e., dwell time)
  # 2) distance difference between consecutive GPS points
  tmp = getMediumSpeed(df)
  
  # A geolocation visited was identified when: 
  # 1) the time difference between two consecutive GPS data points exceeded a time threshold 
  # 2) their distance was less than a distance threshold  
  tmp.geovisited=tmp[tmp$dwellTime >= timeThreshold & tmp$deltaDistance<=distanceThreshold,]
  
  for(i in 1:dim(tmp.geovisited)[1]){
    df[df$X==(tmp.geovisited$xPrev[i]),]$geoVisited=TRUE
  }
  geo.visited = df[df$geoVisited==TRUE,]
  geo.visited = geo.visited[,c("patient","sessionid","Latitude","Longitude")]
  return(geo.visited)
}

getMediumSpeed <- function(df){
  dfTemp = data.frame(df, 
                      xPrev = c(NA, df$X[-nrow(df)]),
                      patientPrev = c(NA, as.character(df$patient[-dim(df)[1]])), 
                      sessionIdPrev = c(NA, df$sessionid[-dim(df)[1]]),    
                      timeStampPrev = c(NA, as.character(df$TimeStamp[-dim(df)[1]])),
                      latPrev = c(NA, df$Latitude[-dim(df)[1]]),
                      longPrev = c(NA, df$Longitude[-dim(df)[1]]))
  # NA raw removal
  dfTemp = dfTemp[dfTemp$patient==dfTemp$patientPrev & 
                    dfTemp$sessionid==dfTemp$sessionIdPrev & 
                    !is.na(dfTemp$patientPrev) &
                    !is.na(dfTemp$latPrev) &
                    !is.na(dfTemp$longPrev), ]
  
  dfTemp$timeStampPrev = as.POSIXct(strptime(dfTemp$timeStampPrev, "%Y-%m-%d %H:%M:%S"))
  dfTemp$dwellTime = as.numeric(difftime(dfTemp$TimeStamp,dfTemp$timeStampPrev,units="secs"))
  dfTemp$deltaDistance = getDistance(dfTemp$Latitude, dfTemp$Longitude, dfTemp$latPrev, dfTemp$longPrev)
  dfTemp$mediumSpeed = (dfTemp$deltaDistance/dfTemp$dwellTime)
  
  return(dfTemp)
}

getDistance <- function(lat1, long1, lat2, long2){
  # Generic function that given two GPS data points, returns their distance in meters
  ER = 6371000 # Earth medium radius in m
  phi1 = lat1 * pi / 180
  phi2 = lat2 * pi / 180
  lam1 = long1 * pi / 180
  lam2 = long2 * pi / 180
  x = (lam2 - lam1) * cos((phi1 + phi2)/2)
  y = (phi2 - phi1)
  distance = ER* sqrt((x*x) + (y*y))
  return(distance)
}

densityBasedMethod <- function(df,timeThreshold,radius ){
  # Patient and sessions
  ps=df[,c("patient","sessionid")]
  ps=ps[!duplicated(ps),]
  geo.visited = NULL
  # For each patient, for each session
  for(row in 1:dim(ps)[1]){
    df.tmp = df[df$patient == ps$patient[row] & df$sessionid == ps$sessionid[row],]
    # For each point 
    j = 1
    while(j <=nrow(df.tmp)){
      # take 10 minutes points
      tenMinPoints = getDeltaTimeMinutesPoints(df.tmp, j, timeThreshold)
      # center calculation
      tenMinPoints = getCenter(tenMinPoints)
      # distance from the center 
      tenMinPoints$distance = getDistance(tenMinPoints$latC, tenMinPoints$longC, tenMinPoints$Latitude, tenMinPoints$Longitude)
      # points in a circle?
      inCircle = arePointsInCircle(tenMinPoints,radius)
      if(inCircle==TRUE){
        geoVisited = data.frame(patient = ps$patient[row],
                                     sessionid = ps$sessionid[row],
                                     Latitude = tenMinPoints$latC[1],
                                     Longitude = tenMinPoints$longC[1])
        geo.visited = rbind(geo.visited,
                            geoVisited)
        
        xEnd = tenMinPoints$X[nrow(tenMinPoints)]
        j = nrow(df.tmp[df.tmp$X<=xEnd,])
      }
      j = j+1
    }
  }
  return(geo.visited)
}

getDeltaTimeMinutesPoints <- function(df, start, timeThreshold){
  if(start+1 > nrow(df)){
    return(df[start,])
  }
  tStart = df$TimeStamp[start]
  tEnd = tStart + timeThreshold
  if(tEnd <= df$TimeStamp[nrow(df)]){
    tmp = df[df$TimeStamp>=tStart & df$TimeStamp<=tEnd,]
  }else{
    tmp = df[start:nrow(df),]
  }
  return(tmp)
}

getCenter <- function(df){
  df$latC = mean(df$Latitude)
  df$longC = mean(df$Longitude)
  return(df)
}

arePointsInCircle <- function(df, radius){
  if(sum(df$distance <= radius) == dim(df)[1]){
    return (TRUE)
  }
  return (FALSE)
}

spaceClustering <- function(geo.visited, distanceThreshold){
  # Create a new variable
  print("Space clustering in progress...")
  geo.visited$placeID = "" 
  # Give an ID (x) to each row
  geo.visited$X = c(1:nrow(geo.visited))
  geo.visited$X = as.factor(geo.visited$X)
  geo.visited$CenterLat = 0
  geo.visited$CenterLong = 0
  patients = unique(geo.visited[,c("patient")])
  for(p in patients){
    clusterID = 1
    tmp = geo.visited[geo.visited$patient==p,]
    center = c(tmp$Latitude[1],tmp$Longitude[1])
    while(nrow(tmp)>0){
      tmp$distance = getDistance(center[1],center[2], tmp$Latitude,tmp$Longitude)
      cluster = tmp[tmp$distance<=distanceThreshold,]
      if(nrow(cluster)==1){
        geo.visited[geo.visited$X==tmp$X[1],]$placeID = paste0("",clusterID)
        geo.visited[geo.visited$X==tmp$X[1],]$CenterLat = center[1]
        geo.visited[geo.visited$X==tmp$X[1],]$CenterLong = center[2]
        tmp = tmp[-1,]
        center = c(tmp$Latitude[1],tmp$Longitude[1])
        clusterID = clusterID+1
      }else{
        cluster$latC = mean(cluster$Latitude)
        cluster$longC = mean(cluster$Longitude)
        if(cluster$latC[1]==center[1] & cluster$longC[1]==center[2]){
          for(z in cluster$X){
            geo.visited[geo.visited$X==z,]$placeID = paste0("",clusterID)
            geo.visited[geo.visited$X==z,]$CenterLat = center[1]
            geo.visited[geo.visited$X==z,]$CenterLong = center[2]
            tmp = tmp[!(tmp$X==z),]
          }
          clusterID = clusterID+1
          if(nrow(tmp)>0){
            center = c(tmp$Latitude[1],tmp$Longitude[1])
          }
        }else{
          center = c(cluster$latC[1],cluster$longC[1])
        }
      }
    }
  }
  clusters= unique(geo.visited[,c("patient","placeID","CenterLat","CenterLong")])
  return(clusters)
}

assignPlaceID <- function(df, places, distanceThreshold){
  print("Classification in progress...")
  df$X = c(1:nrow(df))
  df$placeID=""
  patients <- unique(places[,c("patient")])
  for(p in patients){
    places.tmp = places[places$patient==p,]
    df.tmp = df[df$patient==p,]
    df.tmp$distance = distanceThreshold
    for(i in 1:nrow(places.tmp)){
      location = places.tmp[i,]
      df.tmp$newDistance = getDistance(location$CenterLat[1], location$CenterLong[1], df.tmp$Latitude, df.tmp$Longitude)
      df.tmp$diffDistance = df.tmp$distance - df.tmp$newDistance
      temp = df.tmp[df.tmp$diffDistance>0,]
      if(nrow(temp)>0){
        df.tmp[df.tmp$diffDistance>0,]$placeID = location$placeID[1]
        df.tmp[df.tmp$diffDistance>0,]$distance = df.tmp[df.tmp$diffDistance>0,]$newDistance
      }
    }
    df[df$patient==p,]$placeID = df.tmp$placeID
  }
  df$X = NULL
  return(df)
}

getPlaceList <- function(df,places){
  print("Places visited detection...")
  placeList = NULL
  df$date = as.Date(df$TimeStamp)
  ps = unique(df[,c("patient","sessionid","date")])
  for(i in 1:nrow(ps)){
    p = ps$patient[i]
    s = ps$sessionid[i]
    d = ps$date[i]
    tmp = df[df$patient==p & df$sessionid==s & df$date==d,]
    tmp$x = c(1:nrow(tmp))
    tmp = tmp[tmp$placeID!="",]
    jStart = 1
    nrow(tmp)
    j = 1
    while(j <= (nrow(tmp)-1)){
      cond1 = ((tmp$x[j]+1)!=tmp$x[j+1])
      cond2 = ((j+1)==nrow(tmp)) 
      if(cond1 | cond2){
        intervalTime = interval(tmp$TimeStamp[jStart],tmp$TimeStamp[j])
        place = data.frame(patient = p,
                           sessionid = s,
                           placeID=tmp$placeID[j],
                           latitude = places[places$placeID==tmp$placeID[j] & places$patient==p,]$CenterLat,
                           longitude = places[places$placeID==tmp$placeID[j] & places$patient==p,]$CenterLong,
                           intervalTime = intervalTime)
        placeList = rbind(placeList,place)
        jStart = j+1
      }
      j = j+1
    }
  }
  return(placeList)
}

place_activity_recognition = function(places.visited, key, filesDirectory){
  #places.visited = mergePlacesVisited(places.visited, 10*60)
  places.visited$placeType=""; places.visited$activityType=""; places.visited$activityCategory=""; 
  places.visited$duration = as.duration(places.visited$intervalTime); places.visited$sessionid=NULL
  # Find the home
  print("Home recognition in progress...")
  activities = getHome(places.visited)
  # Find the work place
  print("Work recognition in progress...")
  activities = getWorkPlace(activities, 4*60*60)
  # Find other activities 
  print("Other activity recognition in progress...")
  activities = getOtherActivities(activities,key, filesDirectory)
  # Clear activities with duration less than 5 min
  activities$duration = as.numeric(activities$duration)
  activities = activities[activities$duration>=(5*60),]
  # Dataset cleaning
  activities$duration=NULL
  activities$placeID=NULL
  return(activities)
}

mergePlacesVisited = function(places.visited, threshold){
  ps = unique(places.visited[,c("patient","sessionid")])
  for(i in 1:nrow(ps)){
    tmp = places.visited[places.visited$patient==ps$patient[i] & places.visited$sessionid==ps$sessionid[i],]
    if(nrow(tmp)>1){
      jStart = -1
      for(j in 1:(nrow(tmp)-1)){
        loc1 = as.numeric(tmp$placeID[j])
        loc2 = as.numeric(tmp$placeID[j+1])
        end1 = int_end(tmp$intervalTime[j])
        start2 = int_start(tmp$intervalTime[j+1])
        diff = as.numeric(difftime(start2, end1, units = "secs"))
        if((loc1==loc2) & (diff<=threshold)){
          if(jStart==-1){
            jStart = j
          }
          tmp[jStart:(j+1),]$intervalTime = interval(int_start(tmp$intervalTime[jStart]),int_end(tmp$intervalTime[j+1]))
        }else{
          jStart = -1
        }
      }
      places.visited[places.visited$patient==ps$patient[i] & places.visited$sessionid==ps$sessionid[i],] = tmp
    }
  }
  places.visited$duration = as.numeric(as.duration(places.visited$intervalTime))
  places.visited = unique(places.visited)
  return(places.visited)
}

getHome <- function(places.visited){
  nightPlaces = getNightPlaces(places.visited)
  activities = homeFinder(nightPlaces, places.visited)
  return(activities)
}

getNightPlaces <- function(places.visited){
  patients = unique(places.visited[,c("patient")])
  nightPlaces = NULL
  for(p in patients){
    tmp = places.visited[places.visited$patient == p,]
    for(j in 1:nrow(tmp)){
      start = int_start(tmp$intervalTime[j])
      hour(start)=21;minute(start)=0;second(start)=0
      end = start + (12*60*60)
      interval = interval(start = start, end = end)
      if(int_overlaps(tmp$intervalTime[j],interval)){
        place = data.frame(patient = p, placeID = tmp$placeID[j])
        nightPlaces=rbind(nightPlaces,place)
      }
    }
  }
  return(nightPlaces)
}

homeFinder <- function(nightPlaces, places.visited){
  nightPlaces$patient = as.character(nightPlaces$patient)
  places.visited$patient = as.character(places.visited$patient)
  nightPlaces$timeSum = 0
  pl = unique(nightPlaces[,c("patient","placeID")])

  maxTimes = data.frame(patient = pl$patient, max=0)
  
  for(i in 1:nrow(pl)){
    p = pl$patient[i]
    l = pl$placeID[i]
    nightPlaces[nightPlaces$patient==p & nightPlaces$placeID==l,]$timeSum = sum(places.visited[places.visited$patient==p & places.visited$placeID==l,]$duration)
    maxTimes[maxTimes$patient==p,]$max = max(nightPlaces[nightPlaces$patient==p,]$timeSum)
  }
  
  activities = places.visited
  
  for(i in 1:nrow(maxTimes)){
    p = maxTimes$patient[i]
    t = maxTimes$max[i]
    l = nightPlaces[nightPlaces$patient==p & nightPlaces$timeSum==t,]$placeID[1]
    activities[activities$patient==p & activities$placeID==l,]$activityType = "home"
    activities[activities$patient==p & activities$placeID==l,]$activityCategory = "home"
    activities[activities$patient==p & activities$placeID==l,]$placeType = "residential"
  }
  return(activities)
}

getWorkPlace = function(activities, time){
  activities$date = as.Date(int_start(activities$intervalTime))
  timeThreshold = time
  pl = unique(activities[,c("patient", "placeID", "date")])
  wp = unique(activities[,c("patient", "placeID")])
  wp$totTime = 0
  wp$days = 0
  for(i in 1:nrow(pl)){
    p = pl$patient[i]
    l = pl$placeID[i]
    d = pl$date[i]
    tmp = activities[activities$patient==p 
                     & activities$placeID==l 
                     & activities$date==d
                     & activities$activityType!="home",
                     ]
    if(nrow(tmp)>0){
      totDuration = sum(tmp$duration)
      if(totDuration>timeThreshold){
        wp[wp$patient==p & wp$placeID==l,]$days = wp[wp$patient==p & wp$placeID==l,]$days[1] + 1
        wp[wp$patient==p & wp$placeID==l,]$totTime = wp[wp$patient==p & wp$placeID==l,]$totTime[1] + totDuration
      }
    }
  }
  wp = wp[wp$days>=3,]
  if(nrow(wp)>0){
    patients = unique(wp[,c("patient")])
    for(p in 1:length(patients)){
      item = wp[wp$patient==patients[p],]
      max = item[item$totTime == max(item$totTime),]
      activities[activities$patient==patients[p] & activities$placeID==max$placeID[1],]$placeType = "work"
      activities[activities$patient==patients[p] & activities$placeID==max$placeID[1],]$activityType = "working"
      activities[activities$patient==patients[p] & activities$placeID==max$placeID[1],]$activityCategory = "employment"
    }
  }
  return(activities)
}

getOtherActivities=function(activities, key, filesDirectory){
  
  activities$patient=as.character(activities$patient)
  activities$placeID=as.character(activities$placeID)
  tmp = activities[activities$placeType!="work" & activities$placeType!="residential", ]
  places = unique(tmp[,c("patient","placeID","latitude","longitude")])
  places$patient = as.character(places$patient)
  places$placeID = as.character(places$placeID)

  for(j in 1:nrow(places)){
    place = places[j,]
    # get POIs from Google Place API and Geonames
    placePOI = NULL
    placePOIs = getGooglePlacesPOI(place$latitude, place$longitude, key)
    placePOIs = rbind(placePOIs,
                      getGeonamePOI(place$latitude, place$longitude))
    # get residential POIs from OpenStreetMap

    if(is.null(placePOIs)){
      min = 100
    }else{
      min = min(placePOIs$distance)*2
    }
    placePOIs = rbind(placePOIs,
                      getResidentialPOIs(place$latitude, place$longitude, min))
    if(!is.null(placePOIs)){
      # If there are POIs, it order by distance and take only POIs within 50 meters
      placePOIs = placePOIs[order(placePOIs$distance),]
      #placePOIs = placePOIs[placePOIs$distance<=50,]
      if(nrow(placePOIs)>0){
        # Check the validity of POI and map the activity type and category
        POI = getValidPOI(placePOIs, filesDirectory)
        activities[activities$patient==place$patient[1] & activities$placeID==place$placeID[1],]$placeType = as.character(POI$placeType)
        activities[activities$patient==place$patient[1] & activities$placeID==place$placeID[1],]$activityType = as.character(POI$activityType)
        activities[activities$patient==place$patient[1] & activities$placeID==place$placeID[1],]$activityCategory = as.character(POI$activityCategory)
      }
    }
  }
  return(activities)
}

getGooglePlacesPOI = function(latitude, longitude, key){
  # Variables
  type = getGooglePlaceType()
  i = 1
  POIs = NULL
  # Request submission
  response = googlePlacesAPI_nearbyResearch(key[i], latitude, longitude, type)
  # if status OK
  if(response$status=="OK"){
    if(length(response)>1){
      POIs = getGooglePOIs(response)
    }
  }else if(response$status=="OVER_QUERY_LIMIT"){
    # Change of key and request submission again
    i = i + 1
    response = googlePlacesAPI_nearbyResearch(key[i], latitude, longitude, type)
    
    if(response$status=="OK"){
      POIs = getGooglePOIs(response)
    }else{
      print(response$status)
    }
  }else{
    print(response$status)
  }
  if(!is.null(POIs)){
    POIs$distance = getDistance(latitude, longitude, POIs$placeLatitude, POIs$placeLongitude)
    POIs$flag = "google"
  }
  return(POIs)
}

googlePlacesAPI_nearbyResearch = function(key, latitude, longitude, type){
  url=paste("https://maps.googleapis.com/maps/api/place/nearbysearch/xml?location=",latitude,",",longitude,"&radius=100&types=",type,"&key=",key,sep="")
  # submit the request
  data=getURL(url)
  # convert XML to list
  response=xmlToList(data)
  return(response)
}

getGooglePOIs = function(response){
  POIs = NULL
  for(j in 2:length(response)){
    if(length(response[[j]])>4){
      POIs = rbind(POIs,
                   newPOI(response[[j]]$name, 
                          response[[j]]$type, 
                          (response[[j]]$geometry$location$lat), 
                          (response[[j]]$geometry$location$lng)))
    }
  }
  return(POIs)
}

newPOI = function(placeName, placeType, placeLatitude, placeLongitude){
  if(!is.null(placeName)){
    POI = data.frame(name = as.character(placeName),
                     placeType = as.character(placeType),
                     placeLatitude= as.numeric(placeLatitude),
                     placeLongitude= as.numeric(placeLongitude))
  }else{
    POI = data.frame(name = "no name",
                     placeType = as.character(placeType),
                     placeLatitude= as.numeric(placeLatitude),
                     placeLongitude= as.numeric(placeLongitude))
  }
  return(POI)
}

getGooglePlaceType = function(){
  return("accounting|amusement_park|aquarium|art_gallery|bakery|bank|bar|beauty_salon|bicycle_store|book_store|bowling_alley|cafe|campground|car_dealer|car_rental|car_repair|car_wash|casino|church|city_hall|clothing_store|convenience_store|courthouse|dentist|department_store|doctor|electrician|electronics_store|embassy|finance|florist|funeral_home|furniture_store|grocery_or_supermarket|gym|hair_care|hardware_store|hindu_temple|home_goods_store|hospital|insurance_agency|jewelry_store|laundry|lawyer|library|liquor_store|local_government_office|locksmith|lodging|meal_delivery|meal_takeaway|mosque|movie_rental|movie_theater|museum|night_club|painter|park|pet_store|pharmacy|physiotherapist|place_of_worship|plumber|police|post_office|real_estate_agency|restaurant|roofing_contractor|rv_park|school|shoe_store|shopping_mall|spa|stadium|store|synagogue|travel_agency|university|veterinary_care|zoo")
}

getGeonamePOI = function(latitude, longitude){
  POIs = NULL
  url = paste("api.geonames.org/findNearbyPOIsOSM?lat=",latitude,"&lng=",longitude,"&radius=",0.1,"&username=universityofmanchest") 
  url = URLencode(url)
  # submit the request
  output=getURL(url)
  # convert XML to list
  response=xmlToList(output)
  if(!is.null(response)){
    for(j in 1: length(response)){
      if(!is.null(response[j]$poi)){
        POIs = rbind(POIs,
                     newPOI(response[j]$poi$name,
                            response[j]$poi$typeName,
                            (response[j]$poi$lat),
                            (response[j]$poi$lng)))
      }
    }
    POIs$distance = getDistance(latitude, longitude, POIs$placeLatitude, POIs$placeLongitude)
    POIs$flag = "geonames"
  }
  return(POIs)
}

getResidentialPOIs <- function(latitude, longitude, squareSize){
  # Output variable
  POIs = NULL
  
  src = osmsource_api()
  centerBox = center_bbox(longitude, latitude, squareSize, squareSize)
  map = get_osm(centerBox, source = src)
  ontology = data.frame(key = c("building","building"), value = c("residential", "house"))
  ontology$key = as.character(ontology$key)
  ontology$value = as.character(ontology$value)
  
  for(j in 1:nrow(map$ways$tags)){
    tag = map$ways$tags[j,]
    item = ontology[ontology$key==as.character(tag$k) & ontology$value==as.character(tag$v),]
    if(nrow(item)==1){
      POIs = rbind(POIs,
                   newPOI("building",
                          item$value,
                          latitude,
                          longitude))
    }
  }
  if(!is.null(POIs)){
    POIs$distance = 0
    POIs$flag = "osm"
  }
  return(POIs)
}

getValidPOI = function(placePOIs, filesDirectory){
  placePOIs$placeType = as.character(placePOIs$placeType)
  # Ontology acquisition
  setwd(filesDirectory)
  osmOntology = read.xlsx("Ontology.xlsx", sheetIndex = 1, header = TRUE)
  osmOntology$value = as.character(osmOntology$value)
  googleOntology = read.xlsx("Ontology.xlsx", sheetIndex = 2, header = TRUE)
  googleOntology$place_type = as.character(googleOntology$place_type)
  valid = FALSE
  j = 1
  
  while((!valid) & (j<=nrow(placePOIs))){
    POI = placePOIs[j,]
    if(POI$flag=="google"){
      item = googleOntology[googleOntology$place_type==POI$placeType,]
      if(nrow(item)==1){
        POI$activityType=as.character(item$activity)
        POI$activityCategory= as.character(item$category)
        valid = TRUE
      }
    }else{
      item = osmOntology[osmOntology$value==POI$placeType,]
      if(nrow(item)>0){
        POI$activityType=as.character(item$activity[1])
        POI$activityCategory= as.character(item$category[1])
        valid = TRUE
      }
    }
    j = j+1
  }
  if(!valid){
    POI$placeType=""
    POI$activityType=""
    POI$activityCategory=""
  }
  
  return(POI)
}

getActivityList = function(activities){
  pd = unique(activities[,c("patient","date")])
  activityList = NULL
  for(i in 1:nrow(pd)){
    p = pd$patient[i]
    d = pd$date[i]
    tmp = activities[activities$patient==p & activities$date==d,]
    if(nrow(tmp)>1){
      for(j in 1:(nrow(tmp)-1)){
        at1 = tmp$activityType[j];   
        at2 = tmp$activityType[j+1];
        if(at1!=at2){
          activityList = rbind(activityList,
                               tmp[j,])
        }
      }
    }else{
      activityList = rbind(activityList,
                           tmp)
    }
    
  }
  activityList = activityList[,c("patient","date","placeType","activityType","activityCategory")]
  return(activityList)
}

getDailyActivities = function(activities){
  output = unique(activities[,c("patient","date","activityType", "activityCategory")])
  return(output)
}

getDailyCategoriesActivities = function(activities){
  output = unique(activities[,c("patient","date", "activityCategory")])
  return(output)
}

getPerformance = function(sfd, activityList){
  meanRecall = 0; meanPrecision = 0
  activityList = unique(activityList[,c("patient","date","activityCategory")])
  patients = unique(sfd$patient)
  performance = data.frame(patient=patients, recall=0, precision=0, n_activities=0)
  for(p in patients){
    sfd.tmp = sfd[sfd$patient==p,]
    activityList.tmp = activityList[activityList$patient==p,]
    confusionMatrix = getConfusionMatrix(sfd.tmp, activityList.tmp)
    recall = (confusionMatrix[1,1]/(confusionMatrix[1,1]+confusionMatrix[2,1]))
    precision = (confusionMatrix[1,1]/(confusionMatrix[1,1]+confusionMatrix[1,2]))
    n_activities = confusionMatrix[1,1] + confusionMatrix[2,1]
    performance[performance$patient==p, ]$recall = recall
    performance[performance$patient==p, ]$precision = precision
    performance[performance$patient==p, ]$n_activities = n_activities
    meanRecall = meanRecall + (recall*(n_activities)) 
    meanPrecision = meanPrecision + (precision*(n_activities))
  }
  meanRecall = meanRecall/sum(performance$n_activities)
  meanPrecision = meanPrecision / sum(performance$n_activities)
  sdRecall = wheightedSD(meanRecall,length(patients),performance$n_activities, performance$recall)
  sdPrecision = wheightedSD(meanPrecision,length(patients),performance$n_activities, performance$precision)
  mean = c("mean (sd)", 
           paste0(meanRecall," (",sdRecall,")"),
           paste0(meanPrecision," (",sdPrecision,")"),
           paste0(mean(performance$n_activities), " (", sd(performance$n_activities), " )"))
  performance$patient = as.character(performance$patient)
  performance = rbind(performance, mean)
  return(performance)
}


getConfusionMatrix = function(sfd,activityList){
  matrix = matrix(0,nrow = 2, ncol = 2)
  TP = 0; FN = 0; FP = 0
  dates = unique(sfd$date)
  for(j in 1:length(dates)){
    sfd.tmp = sfd[sfd$date==dates[j],]
    activityList.tmp = activityList[activityList$date==dates[j],]
    categories = unique(sfd.tmp[,c("activityCategory")])
    k = 1
    while(k<=length(categories)){
      sfd.c = sfd.tmp[sfd.tmp$activityCategory==categories[k],]
      act.c = activityList.tmp[activityList.tmp$activityCategory==categories[k],]
      TP = TP + nrow(act.c)
      if(nrow(sfd.c)>=nrow(act.c)){
        FN = FN + (nrow(sfd.c)-nrow(act.c))
      }
      k = k+1
    }
  }
  dates = unique(activityList$date)
  for(j in 1:length(dates)){
    sfd.tmp = sfd[sfd$date==dates[j],]
    activityList.tmp = activityList[activityList$date==dates[j],]
    categories = unique(activityList.tmp[,c("activityCategory")])
    k = 1
    while(k<=length(categories)){
      sfd.c = sfd.tmp[sfd.tmp$activityCategory==categories[k],]
      act.c = activityList.tmp[activityList.tmp$activityCategory==categories[k],]
      if(nrow(act.c)>nrow(sfd.c)){
        FP = FP + (nrow(act.c)-nrow(sfd.c))
      }
      k = k+1
    }
  }
  matrix[1,1]=TP
  matrix[1,2]=FP
  matrix[2,1]=FN
  return(matrix)
}

wheightedSD = function(mean, M, w, x){
  sd = 0
  for(i in 1:length(x)) {
    diff = (x[i]-mean)^2
    sd = sd + (w[i]*diff)
  }
  den = ((M-1)/M)*sum(w)
  sd = sd/den
  sd = sqrt(sd)
  return(sd)
}

getDataSummary = function(df, sfd){
  df$date = as.Date(df$TimeStamp)
  df$patient=as.character(df$patient)
  patients = unique(df$patient)
  summary = data.frame(participant = patients)
  summary$participant = as.character(summary$participant)
  summary$days = 0
  summary$hours = 0
  summary$activities = 0
  sfd = getListActivitiesSFD(sfd)
  for(i in 1:length(patients)){
    p = summary$participant[i]
    tmp = df[df$patient==p,]
    summary$days[i] = length(unique(tmp$date))
    summary$hours[i] = getTotalHoursParticipant(tmp)
    summary$activities[i] = nrow(sfd[sfd$patient==p,])
  }
  sum = c("TOTAL",
          sum(summary$days),
          sum(summary$hours),
          sum(summary$activities))
  mean = c("mean (sd)",
           paste0(mean(summary$days), " (",sd(summary$days),")"),
           paste0(mean(summary$hours), " (",sd(summary$hours),")"),
           paste0(mean(summary$activities), " (",sd(summary$activities),")"))
  
  summary$days = as.character(summary$days)
  summary$hours = as.character(summary$hours)
  summary$activities = as.character(summary$activities)
  summary = rbind(summary,
                  sum,
                  mean)
  return(summary)
}

getTotalHoursParticipant = function(df){
  dates = unique(df[,c("date")])
  hours = 0
  for(j in 1:length(dates)){
    d = dates[j]
    tmp = df[df$date==d,]
    tmp = tmp[order(tmp$TimeStamp),]
    h = as.numeric(difftime(tmp$TimeStamp[nrow(tmp)],tmp$TimeStamp[1],units="hours"))
    hours = hours + h
  }
  return(hours)
}

getListActivitiesSFD = function(sfd){
  dates = unique(sfd[,c("patient","date")])
  output = NULL
  for(k in 1:nrow(dates)){
    d = dates$date[k]
    p = dates$patient[k]
    tmp = sfd[sfd$patient==p & sfd$date==d,]
    tmp = tmp[,c("patient","date","activityCategory")]
    tmp = rbind(tmp,
                data.frame(patient=tmp$patient[1],
                           date = tmp$date[1],
                           activityCategory = "home"))
    output = rbind( output,
                    unique(tmp))
  }
  return(output)
}

plotMap = function(coordinates){
  # costant
  variation = 0.0009
  rangeLat = range(coordinates$latitude)
  rangeLong = range(coordinates$longitude)
  location = c(rangeLong[1]-variation,rangeLat[1]-variation,rangeLong[2]+variation,rangeLat[2]+variation)
  map = get_map(location = location, source = "osm")
  map1 = ggmap(map)
  map1 = map1 + geom_point(data = coordinates, aes(x=longitude, y=latitude), size = 5, col = "red") 
  map1
}


getNPV = function(activityList){
  pd = unique(activityList[,c("patient","date")])
  activityList = activityList[!(activityList$activityCategory=="home"),]
  NPV = pd
  NPV$NPV = 0
  activityList = 
  for(j in 1:nrow(pd)){
    p = pd$patient[j]
    d = pd$date[j]
    tmp = activityList[activityList$patient==p & activityList$date==d,]
    NPV[NPV$patient==p & NPV$date==d,]$NPV = nrow(tmp)
  }
  return(NPV)
}

getOutOfHomeHours = function(activities){
  pd = unique(activities[,c("patient","date")])
  outOfHomeHours = pd
  outOfHomeHours$hours = 0
  for(j in 1:nrow(pd)){
    p = pd$patient[j]
    d = pd$date[j]
    tmp = activities[activities$patient==p & activities$date==d,]
    if(tmp$activityCategory[1]=="home"){
      start = int_end(tmp$intervalTime[1])
    }else{
      start = int_start(tmp$intervalTime[1])
    }
  }
}

activitiesAcquisition = function(datasetType, pipeline){
  activitiesLabel = paste0(datasetType," activities ",pipeline, " method.xlsx")
  activities = read.xlsx(activitiesLabel, sheetIndex = 1)
  dates <- t(data.frame(strsplit(as.character(activities$intervalTime),"--")))
  start= as.POSIXct((strptime(dates[,1], "%Y-%m-%d %H:%M:%S GMT")))
  end= as.POSIXct((strptime(dates[,2], "%Y-%m-%d %H:%M:%S GMT")))
  activities$intervalTime = interval(start, end)
  activities$NA.=NULL
  activities$activityCategory=as.character(activities$activityCategory)
  activities$method = pipeline
  activities = activities[,c("patient","latitude","longitude","activityType","activityCategory","method","date")]
  return(activities)
}


plotMapComparison <- function(df){

  m = unique(df$method)
  # costant
  variation = 0.0025
  
  rangeLat = range(df$latitude)
  rangeLong = range(df$longitude)
  
  location = c(rangeLong[1]-variation,rangeLat[1]-variation,rangeLong[2]+variation,rangeLat[2]+variation)
  print(location)
  map = get_map(location = location, source = "osm")
  map = ggmap(map)
  
  map1 =  map + geom_point(data = df[df$method==m[1],], aes(x=longitude, y=latitude), size = 6, colour = "red" ) + ggtitle("Time-based pipeline")
  map2 =  map + geom_point(data = df[df$method==m[2],], aes(x=longitude, y=latitude), size = 6, colour = "red" ) + ggtitle("Density-based pipeline")
  map3 =  map + geom_point(data = df[df$method==m[3],], aes(x=longitude, y=latitude), size = 6, colour = "red" ) + ggtitle("Combination of pipelines")
  
  fileName = paste0("Comparison_",".png")
  png(filename = fileName, width = 1600, height = 1600, units = "px" )
  multiplot(map1,map2,map3, cols = 3)
  dev.off()
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


