
library(tidyverse)

# get the file path
RPROJ <- list(PROJHOME = normalizePath(getwd()))
attach(RPROJ)
path.p<-file.path(PROJHOME, 'RespirometryFiles')

# bring in the oxygen files
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
#basename above removes the subdirectory name from the file
#add file names that include the subdirectory name (note, these are the same for this example, but I often have lots of subfolders for different Runs)
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) #list all csv file names in the folder and subfolders

#Load Sample Info
#Sample.Info <- read.csv(file=paste0(path.p,"/../Panama MetaData/Nubbin_Sample_Info_T0_Panama_QC.csv"), header=T) #read sample.info data
Sample.Info <- read.csv(file=paste0(path.p,"/../metadata.csv"), header=T) #read sample.info data


# load surface area data
SA <- read.csv(file=paste0(path.p,"/../Burn Volumes.csv"), header=T) #read sample.info data
# add 610 ml to the NAs in volume (the blanks)
#Calculat the volume of water
SA$Volume<-610-SA$Initial.Volume..ml.
#Sample.Info$Volume[which(is.na(Sample.Info$Volume))]<-610
# add 0's for the "not blanks"
SA$Blank[is.na(SA$Blank)]<-0

# joint the sample info and surface area and volume measurements
Sample.Info<-left_join(Sample.Info, SA)

# make start and stop times real times
Sample.Info$Start.time <- as.POSIXct(Sample.Info$Start.time,format="%H:%M:%S", tz = "") #convert time from character to time
Sample.Info$Stop.Time <- as.POSIXct(Sample.Info$Stop.Time,format="%H:%M:%S", tz = "") #convert time from character to time

#Add names for photosynthesis or respiration for for loop
#PR<-c('Photo','Resp')


# read in the photo data
Photo.R<-read.csv('Output/Photo.R.csv')

# read in the files that have messy respiration rates

#Bad<-read.csv('Quality Control.csv')
Bad<-read.csv('QC_noblanks.csv')

# pull out the bad ones
bad.rows<-which(Photo.R$ID %in% Bad$ID)
Photo.R<-Photo.R[-bad.rows,]

#Photo.R$Fragment.ID.full<-Photo.R$ID
#Photo.R$Fragment.ID<-NULL
Photo.R$Light_Dark<-NULL

Photo.R<-left_join(Photo.R, Sample.Info)
#Convert sample volume to mL
Photo.R$Vol.L <- Photo.R$Volume/1000 #calculate volume

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Photo.R$umol.sec <- Photo.R$umol.L.sec*Photo.R$Vol.L

#Account for blank rate by temerature
#convert character columns to factors
Photo.R <- Photo.R %>%
  mutate_if(sapply(., is.character), as.factor)
# make the blank column a factor
Photo.R$Blank<-as.factor(Photo.R$Blank)

photo.blnk <- aggregate(umol.sec ~ Species*Temp.Cat*Light_Dark*Blank, data=Photo.R, mean)
# pull out only the blanks
#photo.blnk<-photo.blnk[photo.blnk$Species=='BK',]
photo.blnk<-photo.blnk[photo.blnk$Blank==1,]

# remove the species column and join with the full data set
#photo.blnk$Species<-NULL
# remove the blank column
photo.blnk$Blank<-NULL

colnames(photo.blnk)[4]<-'blank.rate' # rename the blank rate 
# join the blank data with the rest of the data
Photo.R<-left_join(Photo.R, photo.blnk)

# subtract the blanks######################
Photo.R$umol.sec.corr<-Photo.R$umol.sec-Photo.R$blank.rate

#### Normalize to organic biomass (ash free dry weight)#####

#Calculate net P and R
Photo.R$umol.cm2.hr <- (Photo.R$umol.sec.corr*3600)/Photo.R$AFDW #mmol cm-2 hr-1

#Photo.R<-Photo.R[complete.cases(Photo.R),] # remove NAs and blanks
Photo.R<-Photo.R[Photo.R$Blank==0,]

#make respiration positive
#Photo.R$umol.cm2.hr[Photo.R$PR=='Respiration']<-abs(Photo.R$umol.cm2.hr[Photo.R$PR=='Respiration'])
Photo.R$umol.cm2.hr<- -Photo.R$umol.cm2.hr

# log the rates
Photo.R$Rate.ln<-log(Photo.R$umol.cm2.hr+0.1)
#remove empty rows
#Photo.R<-Photo.R[-which(is.na(Photo.R$ID)),]

#ggplot(Photo.R, aes(x=Temp.C, y=umol.cm2.hr,  group=c(Species), col = Organism.ID))+
#geom_point(aes(shape=Species), position = position_dodge(width = 0.2), size=4)+
#ylim(0,1.5)+  
#facet_wrap(~ Species, labeller = labeller(.multi_line = FALSE))


# pull our the site names
Photo.R<-separate(data = Photo.R, col = Organism.ID, sep= "_", into = c('sp','Location','number'), remove = FALSE)  %>%
  select(-one_of('sp','number'))

# remove the two rows with the rates over 80 that are clearly outliers
Photo.R<-Photo.R[-which(Photo.R$umol.cm2.hr>80),]
Photo.R<-Photo.R[-which(Photo.R$umol.cm2.hr>40 & Photo.R$Location=='Ibbet'),]
Photo.R<-Photo.R[-which(Photo.R$umol.cm2.hr>50 & Photo.R$Location=='Bart'),]

# remove data for Bart, doug, and ibbet after the urchins died (i.e. rates flatlined)
Photo.R<-Photo.R[-which(Photo.R$Temp.Cat>36 & Photo.R$Location=='Bart'),]
Photo.R<-Photo.R[-which(Photo.R$Temp.Cat>38 & Photo.R$Location=='Doug'),]
Photo.R<-Photo.R[-which(Photo.R$Temp.Cat>40 & Photo.R$Location=='Ibbet'),]

write.csv(Photo.R, 'GalapagosRates.csv', row.names = FALSE) # export all the uptake rates

#take the means for the plot
PhotoMeans<- Photo.R %>%
  group_by(Location, Temp.Cat)%>%
  dplyr::summarise(rates.mean = mean(umol.cm2.hr), se = sd(umol.cm2.hr)/sqrt(n()))


# plot the raw data with the means on top
ggplot()+
  theme_bw()+  
  geom_point(data=Photo.R, aes(x=Temp.Cat, y=umol.cm2.hr, group=c(Location), col = Location, alpha = 0.05), position = position_dodge(width = 0.2), size=4)+
  geom_point(data=PhotoMeans, aes(x=Temp.Cat, y=rates.mean, group=c(Location)), position = position_dodge(width = 0.2), size=4)+
  geom_line(data = PhotoMeans,  aes(x=Temp.Cat, y=rates.mean, group=c(Location)), position = position_dodge(width = 0.2), size=4)+
  geom_errorbar(data = PhotoMeans, aes(x = Temp.Cat, ymin=rates.mean-se, ymax=rates.mean+se, group=c(Location)), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~ Location, labeller = labeller(.multi_line = FALSE), scales = 'free_y')+
  ggsave('RespirationRates.png')
