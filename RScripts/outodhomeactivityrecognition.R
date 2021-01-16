# Out of home activity recognition
rm(list=ls())

# **** This R Script implements methods described in:
# S. Difrancesco et al., "Out-of-Home Activity Recognition from GPS Data in Schizophrenic Patients," 
# 2016 IEEE 29th International Symposium on Computer-Based Medical Systems (CBMS), Dublin, 2016, 
# pp. 324-328, doi: 10.1109/CBMS.2016.54.

# Before run the code, specify which pipeline you want to run:
# (a) = time-based 
# (b) = density-based
# (c) = timePlusDensity-based
pipeline = "timePlusDensity-based"

# Directories
# Directory of R scripts - Function script
functionDir = "C:\\Users\\sonia\\Documents\\ActivityRecognitionGPSData\\R Scripts\\FunctionScript.R"
# Directory of GPS data folder - Inside of this folder I created a folder for each participant of the study 
directory <- "C:\\Users\\sonia\\Documents\\ActivityRecognitionGPSData\\GPS_data"
# Directory of other files: ontology and 'social functioning' files 
filesDirectory <- "C:\\Users\\sonia\\Documents\\ActivityRecognitionGPSData\\Files"
# sourcing file with functions
source(functionDir)
# libraries
library(osmar)
library(ggmap)
library(ggplot2)
library(xlsx)
library(chron)
library(lubridate)
library(plotKML)
library(RCurl)
library(XML)

# key Google Place API
key= c("AIzaSyDWs7eStEfQRGG8tuNDheo2SJR8ooPjr14",
       "AIzaSyCdyXaICjKXqefkUUzebnw7A6wDvcQac7",
       "AIzaSyCl_RMVOmZfOVdfj8Umn9RRytSHdSIIV3k",
       "AIzaSyAoZSLAoWUxPlOBGE3EeHUlxJ9arrPJt90")

# GPS data acquisition
df = getDataset(directory)

if(pipeline=="time-based"){
  # Pipeline 1: time-based method
  # ______________________________
  # Variables setting
  timeThreshold=60*10
  distanceThreshold = 50
  radius = 50
  # Step 1: Geolocations visited detection
  # ______________________________________
  # 'timeBasedMethod' is a function that identify the geolocations visited 
  geo.visited = timeBasedMethod(df,timeThreshold,distanceThreshold) 
}else if(pipeline=="density-based"){
  # Pipeline 2: density-based method
  # ______________________________
  # Variables setting
  timeThreshold=60*10
  radius = 50
  # Step 1: Geolocations visited detection
  # ______________________________________
  # 'densityBasedMethod' is a function that identify the geolocations visited 
  geo.visited = densityBasedMethod(df,timeThreshold,radius) 
}else if(pipeline=="timePlusDensity-based"){
  timeThreshold=60*10
  distanceThreshold = 50
  radius = 50
  geo.visited = timeBasedMethod(df,timeThreshold,distanceThreshold) 
  geo.visited = rbind(geo.visited,
                      densityBasedMethod(df,timeThreshold,radius) )
}
# Step 2: Places visited identification
# ______________________________________
# 'spaceClustering' is a function that cluster geolocations visited in places
# 'assignPlaceID' is the function that classify GPS data points with place ID 
places = spaceClustering(geo.visited, radius)    
df.places = assignPlaceID(df,places,radius)
places.visited = getPlaceList(df.places, places)
# Step 3: Type of Place and Activity recognition
# ______________________________________________
activities = place_activity_recognition(places.visited, key, filesDirectory)


# Step 4: extraction of NPV and OOH-hours
activityList = getActivityList(activities)
getNPV = getNPV(activityList)
