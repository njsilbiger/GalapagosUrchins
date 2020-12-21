##TPC curves, maps, and temperature analysis for Urchin Galapagos paper

rm(list=ls())

##Install packages
# load packages
library(nls.multstart)
library(broom)
library(purrr)
library(tidyverse)
library(nlstools)
library(nls2)
library(grid)
library(gridExtra)
library(cowplot)
library(lubridate)
library(directlabels)
library(rgdal)
library(rgeos)
library(ggthemes)
library(ggsn)
library(sp)
library(ggrepel)
library(raster)
library(rgdal)
library(patchwork)
#load data

photo.data <- read.csv("GalapagosRates.csv")
photo.data$X <- NULL
View(photo.data)
glimpse(photo.data)

#decide color scheme for the plots
#cols<-c("#99817b", "#F2C3A7", "#FEF3E1", "#C489B9")
cols<-c("#073e3e", "#c35119", "#f896b0", "#e4e0ca")

 # remove the NAs from the data
#photo.data<-photo.data[-which(is.na(photo.data$umol.cm2.hr)),]

# remove the three organisms that are wrong (had too much messy respiration files)
remove<-c('Egala_Bart_1','Egala_Ibbet_1','Egala_Ibbet_3','Egala_Botel_2', 'Egala_Corm_10')

bad.ID<-which(photo.data$Organism.ID %in% remove)
photo.data<-photo.data[-bad.ID,]

mydata <- photo.data

mydata$log.rate <- log(mydata$umol.cm2.hr)  #logging and adding 0.1 because a log of zero does not exist

# convert temp to K
mydata$K<-mydata$Temp.C + 273.15

# define the Sharpe-Schoolfield equation
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp))) #units are eV/K, electrovolts/Kelvin
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}


# fit over each set of groupings

#droplevels(mydata$Organism.ID)

mydata$Location<-as.character(mydata$Location)
mydata$Organism.ID<-as.character(mydata$Organism.ID)

fits <- mydata %>%
  group_by(Organism.ID, Location) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~ nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 26),
                                                data = .x,
                                                iter = 1000,
                                                start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                                                start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                                supp_errors = 'Y',
                                                na.action = na.omit,
                                                lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))



#broom, models over and over again and purr for lists
#make dplyr code and export the slope, intercept and R2 values
#get r2, extract predit and pull out and join lists, shows o vs predcited and r2
PredictedEgala_Bart_1 <- predict(fits$fit[[1]])
ObservedEgala_Bart_1 <- mydata$umol.cm2.hr[mydata$Organism.ID == "Egala_Bart_10"]

po <- lm(PredictedEgala_Bart_1 ~ ObservedEgala_Bart_1)
summary(po)
plot(ObservedEgala_Bart_1,PredictedEgala_Bart_1)
abline(po)
legend("topleft", bty="n", legend=paste("r2 =", format(summary(po)$adj.r.squared, digits=4)))



# look at a single fit
summary(fits$fit[[1]])

# look at output object
#select(fits, Organism.ID, data, fit)  

# get summary info
info <- fits %>%
  unnest_legacy(fit %>% map(glance))

# get params
params <- fits %>%
  unnest_legacy(fit %>% map(tidy))

# get confidence intervals
CI <- fits %>% 
  unnest_legacy(fit %>% map(~ confint2(.x) %>%
                       data.frame() #%>%
                  #     rename(., conf.low = X2.5.., conf.high = X97.5..)
                     )) %>%
  group_by(., Organism.ID) %>%
  mutate(., term = c('lnc', 'E', 'Eh', 'Th')) %>%
  ungroup()

colnames(CI)[3:4]<-c("conf.low", "conf.high") # rename columns

# merge parameters and CI estimates
params <- merge(params, CI, by = intersect(names(params), names(CI)))

# get predictions
preds <- fits %>%
  unnest_legacy(fit %>% map(augment))

#select(info, fragment.ID, logLik, AIC, BIC, deviance, df.residual)

# new data frame of predictions, do this to set a sequence to make a smooth curve with your prediction points
new_preds <- mydata %>%  do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 150), stringsAsFactors = FALSE)) #setting a specific sequence so you can have a smooth curve

# max and min for each curve
max_min <- mydata %>% group_by(Organism.ID) %>%
  dplyr::summarise(., min_K = min(K), max_K = max(K)) %>%
  ungroup()

# create new predictions
preds2 <- fits %>%
  unnest_legacy(fit %>% map(augment, newdata = new_preds)) %>%
  merge(., max_min, by = "Organism.ID") %>%
  group_by(., Organism.ID) %>%
  filter(., K > unique(min_K) & K < unique(max_K)) %>%
  dplyr::rename(., ln.rate = .fitted) %>%
  ungroup()

#want to do ggplot where we look at fragments individually
#reorder the sites
mydata$Location<- factor(mydata$Location, levels=c('Bart','Ibbet','Corm','Doug','Espi','Botel'))
preds2$Location<- factor(preds2$Location, levels=c('Bart','Ibbet','Corm','Doug','Espi','Botel'))

# rename the labels

ggplot() +
  geom_point(aes(K - 273.15, log.rate, col = Location), size = 2, mydata) +
  geom_line(aes(K - 273.15, ln.rate, col = Location, group = Organism.ID), alpha = 0.5, preds2) +
  facet_wrap(~ Organism.ID, labeller = labeller(.multi_line = FALSE), scales = 'free_y') +
  #scale_colour_manual(values = c('green4', 'blue', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab(expression(paste("Respiration Rates (Log " *mu* "mol O"[2], "  "*g^-1 , "  "*hr^-1*")"), sep = " ") )+
  xlab('Temperature (ºC)') +
  theme(legend.position = c(0.91, 0.85))+
 # scale_color_manual(values = cols)+
  labs(color = "Rate Type")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = 'none')

# plot all values in TPCs and remove Doug and Ibbet
mydata<-mydata %>%
  filter(Location != 'Ibbet', Location != 'Doug') %>%
  droplevels() %>%
  mutate(Location = factor(Location,levels = c("Botel","Corm", "Espi", "Bart")),
         Location = factor(Location, labels = c(" La Botella"," Punta Cormorant" ," Punta Espinosa" ," Bartolomé" )))

preds2<-preds2%>%
  filter(Location != 'Ibbet', Location != 'Doug')%>%
  droplevels() %>%
  mutate(Location = factor(Location,levels = c("Botel","Corm", "Espi", "Bart")),
    Location = factor(Location, labels = c(" La Botella"," Punta Cormorant" ," Punta Espinosa" ," Bartolomé" )))


ggplot() +
  geom_point(aes(K - 273.15, log.rate, col = Location), size = 2, mydata) +
  geom_line(aes(K - 273.15, ln.rate, col = Location, group = Organism.ID), alpha = 0.5, preds2) +
  facet_wrap(~ Location, labeller = labeller(.multi_line = FALSE)) +
  scale_color_manual(values = cols)+
    #scale_color_manual(values = c("#323695", "#abd9e9", "#f4a582", "#b2182b"))+
  #scale_colour_manual(values = c('green4', 'blue', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab(expression(paste("Respiration Rates (Log " *mu* "mol O"[2], "  "*g^-1 , "  "*hr^-1*")"), sep = " ") )+
  xlab('Temperature (ºC)') +
  guides(color = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(filename = "Output/MSPlots/TPCcurves.pdf", device = "pdf", width = 6, height = 6)

# function for calculating Topt

get_topt <- function(E, Th, Eh){
  return((Eh*Th)/(Eh + (8.62e-05 *Th*log((Eh/E) - 1))))
}


# calc topts for all 
Topt_data <- params %>%
  dplyr::select(Organism.ID, term, estimate,Location) %>%
  spread(term, estimate) %>%
  mutate(Topt = get_topt(E, Th, Eh)) %>%
  group_by(., Location, Organism.ID)

#get temerature back in celcius not K
Topt_data$Topt <- Topt_data$Topt - 273.15 

#anova function
Topt.mod <- lm(Topt~Location, data=Topt_data)

#check for normality, use normality plots

qqnorm(resid(Topt.mod))
qqline(resid(Topt.mod))

#check heteroscisity with boxplots

boxplot(resid(Topt.mod)~Topt_data$Location)

#high R and low show inconsistent variances, may need to do weighted regression in the future

anova(Topt.mod)
summary(Topt.mod)
TukeyHSD(aov(Topt.mod))

# plot all the TPCs
ggplot() +
  geom_point(aes(K - 273.15, log.rate, col = Location), size = 2, mydata) +
  geom_line(aes(K - 273.15, ln.rate, col = Location, group = Organism.ID), alpha = 0.5, preds2) +
  facet_wrap(~ Location, labeller = labeller(.multi_line = FALSE)) +
  scale_color_manual(values = cols)+
  #scale_color_manual(values = c("#323695", "#abd9e9", "#f4a582", "#b2182b"))+
  #scale_colour_manual(values = c('green4', 'blue', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab(expression(paste("Respiration Rates (Log " *mu* "mol O"[2], "  "*g^-1 , "  "*hr^-1*")"), sep = " ") )+
  xlab('Temperature (ºC)') +
  guides(color = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(filename = "Output/MSPlots/TPCcurves.pdf", device = "pdf", width = 6, height = 6)

data.summary<-Topt_data %>%
  group_by(Location) %>% #tells to group by these two factors
  dplyr::select(-c(Eh, Th)) %>% # only keep the E, lnc, and Topt metrics
  dplyr::summarise_if(is.numeric, list(~mean(.), ~var(.)))
  #dplyr::summarise_if(is.numeric, list(~mean(.), ~var(.), ~sd(.)/sqrt(n())))

# colnames(data.summary)[12:16]<-c("E_SE","Eh_SE","lnc_SE","Th_SE","Topt_SE")
#  dplyr::summarise(mean=mean(Topt), se=sd(Topt)/sqrt(n()), var = var(Topt)) #calculates mean and s.e.
data.summary

#reorder the sites
data.summary$Location<- factor(data.summary$Location, levels=c('Bart','Ibbet','Corm','Doug','Espi','Botel'))


# bring in the temperature data
#tempdata.insi<-read.csv('HoboLoggersFiles/logger_gal_daily_means.csv')
#tempdata<-read.csv('MUR_sst_gal.csv')
tempdata.insi<-read.csv('HoboLoggersFiles/CSV files/AllTempData.csv')

# # make a date
 tempdata.insi$t<-mdy_hm(tempdata.insi$t)
# # add column for day
 tempdata.insi$day<-date(tempdata.insi$t)

  #Make a plot with the raw data for each site for the last two months
 rawplot<-tempdata.insi %>%
   mutate(Location = factor(Location, levels = c("Botel","Corm","Espi","Bart")), # we only had in situ temperature data for these 4 sites so they were the only ones included in the analysis
     LocationNice = recode(Location, 
                                Botel = "La Botella",
                                Corm = "Punta Cormorant",
                                Espi = "Punta Espinosa",
                                Bart = "Bartolomé")) %>%
   filter(Location!="")%>% # remove the empty site
   filter(t>max(t)-months(2))%>% # only use the last 30 days
   ggplot(aes(x = t, y = temp, group = Location, color = LocationNice))+
   geom_line()+
   xlab('Date')+
   ylab(expression("Temperature " ( degree*C)))+
 #  geom_dl(aes(label = LocationNice), method =  list("last.bumpup", cex = 1))+ # add labels at end of plot
   scale_color_manual(values = cols)+
  # scale_color_manual(values = c("#323695", "#abd9e9", "#f4a582", "#b2182b"))+
   guides(color = FALSE) +
   theme_bw()
  # labs(colour = "Location")+ 
  # ggsave("Output/MSPlots/TempRaw2.pdf", width = 8, height = 5)
   #geom_dl(aes(label = LocationNice), method =  list("last.bumpup", cex = 1)) # add labels at end of plot
   
 bplot<-tempdata.insi %>%
   mutate(Location = factor(Location, levels = c("Botel","Corm","Espi","Bart")),
          LocationNice = recode(Location, 
                                Botel = "La Botella",
                                Corm = "Punta Cormorant",
                                Espi = "Punta Espinosa",
                                Bart = "Bartolomé")) %>%
   filter(Location!="")%>% # remove the empty site
   filter(t>max(t)-months(2))%>% # only use the last 30 days
   ggplot(aes(x = LocationNice, y = temp, fill = LocationNice))+
   geom_boxplot()+
  # ylab(expression("Temperature " ( degree*C)))+
   xlab("")+
   ylab("")+
   #  geom_dl(aes(label = LocationNice), method =  list("last.bumpup", cex = 1))+ # add labels at end of plot
 #  scale_fill_manual(values = c("#323695", "#abd9e9", "#f4a582", "#b2182b"))+
   scale_fill_manual(values = cols)+
    guides(fill = FALSE) +
   theme_bw()+
   theme(axis.text.y = element_blank())
   
 rawplot+bplot +plot_annotation(tag_levels = "A")+
   ggsave("Output/MSPlots/combinedtemp.pdf", height = 5, width = 10)
   
# plot the daily maximim over time
 Dailymax<-tempdata.insi %>%
   filter(Location!="")%>% # remove the empty site
   filter(t>max(t)-months(2))%>% # only use the last 30 days
   group_by(Location, day) %>%
   dplyr::summarise(max = max(temp))
 
 #reorder sites 
 Dailymax$Location<-factor(Dailymax$Location, levels=c('Bart','Corm','Espi','Botel'), labels = c(" Bartolomé", " Punta Cormorant" ," Punta Espinosa" ," La Botella" )) 
 
 
 Dailymax%>%
   ggplot(aes(x = day, y = max, group = Location, color = Location))+
   geom_line(lwd = 1)+
   xlab('Date')+
   ylab(expression("Daily Maximum Temperature   "(degree*C)))+
   scale_color_manual(values = cols)+
   #scale_color_brewer(palette="Set2")+
   guides(color = FALSE) +
   theme_bw() +
   scale_x_date(date_labels = "%b %d", breaks = "2 weeks", limits = c(as.Date("2018-06-15"),as.Date("2018-08-15")+days(20)))+ # expand the x-axis by 5 days
   geom_dl(aes(label = Location), method =  list("last.bumpup", cex = 1))+ # add labels at end of plot
   ggsave(filename = "Output/MSPlots/TempTimeSeries.pdf", width = 7, height = 4)
 
 # take the average daily max for analysis
 Dailymax.mean<-Dailymax%>%
   group_by(Location)%>%
   summarise(meandailymax=mean(max))
 
# Dailymax$Location<-as.factor(c("Bart","Corm","Espi", "Botel"))
 
 # 90th percentile temperature and other summaries from the raw data
tempsummary<-tempdata.insi %>%
   filter(Location!="")%>% # remove the empty site
   filter(t>max(t)-months(2))%>% # only use the last 30 days
   group_by(Location) %>%
   dplyr::summarise(Q90 = quantile(temp, probs=0.95),mean.sst = mean(temp, na.rm=TRUE), 
                    max.sst = max(temp, na.rm=TRUE), var.sst = var(temp, na.rm=TRUE), 
                    range.sst = max(temp, na.rm=TRUE) - mean(temp, na.rm=TRUE),
                    min = min(temp,na.rm = TRUE))
 # left_join(.,Dailymax.mean)
 
# join with the thermal optimum data
data.summary.all<-left_join(data.summary, tempsummary)%>%
  #select(-c("mean.sst","var.sst","range.sst")) %>% # only do mean teperature for now
  gather("Parameter", "Metric",-c(Location, Q90, mean.sst, max.sst, var.sst, range.sst, min)) %>%
  filter(Location != 'Ibbet', Location != 'Doug') %>%
  separate(col = Parameter, into = c("Parameter", "Stat"),sep = "_") # split the stat and the parameter name

#change facot for the stat so that mean and var are capitalized
data.summary.all$Stat<-as.factor(data.summary.all$Stat)
levels(data.summary.all$Stat)<-c("Mean","Variance")
  

ggplot(data.summary.all, aes(x = mean.sst, y = Metric), group = Parameter)+
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
  #geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
 # facet_wrap(~Parameter, scales = 'free', ncol = 2)+
  facet_grid(Parameter~ Stat, scale = "free")+
  theme_bw()+
  ggtitle('Two months of temperature data')

ggplot(data.summary.all, aes(x = range.sst, y = Metric), group = Parameter)+
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
  #geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  facet_wrap(Parameter~ Stat, scale = "free")+
  ggtitle('Two months of temperature data')


ggplot(data.summary.all, aes(x = Q90, y = Metric), group = Parameter)+
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
  #geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  facet_wrap(Parameter~ Stat, scale = "free")+
  ggtitle('Two months of temperature data')


## Make plots with regression lines
ggplot(data.summary.all, aes(x = Q90, y = Metric, label = Location, col = Location), group = Parameter)+
  geom_point(position="dodge", size=2) +
  scale_color_manual(values = cols)+
  theme(legend.text=element_text(size=rel(1))) + #makes legend elements larger
  xlab(expression("95th Percentile Temperature  " (degree*C)))+
  geom_smooth(method = "lm", se=FALSE, col = 'grey')+
  theme_bw()+
  facet_wrap(~Parameter+Stat, scale = "free_y", ncol = 2,
             strip.position = "left",
             labeller = as_labeller(c(E = "E", lnc = "b(Tc)", Topt = "Topt", Mean = "", Variance = "", sep = ""), multi_line = FALSE ) 
           )  +
  ylab(NULL) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none")+
  ggsave(filename = "Output/MSPlots/MetricsVsTemp.pdf", device = "pdf", width = 6, height = 6, useDingbats = FALSE)
  #ggtitle('Two months of temperature data (in situ)')


# run stats for each pair and make datatable
results<-data.summary.all %>%
  nest(-c(Parameter, Stat)) %>% 
  mutate(fit = map(data, ~ lm(Metric~Q90, data = .)), 
         results = map(fit,glance))%>%
  unnest_legacy(results)%>%
  dplyr::select(-c(data,fit))

#effect sizes
data.summary.all %>%
  group_by(Parameter, Stat) %>%
  do(allfits = tidy(lm(Metric ~Q90, data = .)))%>%
  unnest(allfits) 

# look at AICs of all metrics
data.summary.all %>%
  pivot_longer(cols = c("Q90":"min"), names_to = "tempparams", values_to = "values") %>%
  group_by(Parameter, Stat, tempparams) %>%
  do(allfits = glance(lm(Metric ~values, data = .)))%>%
  unnest(allfits) %>%
  View()

# wriate a data table 
write.csv(x = results,file = 'Output/MSPlots/lmresults.csv')
         

## Make a plot with population TPC curves with bootstrapped confidence internvals
# run nls.multstart on each curve of the original data ####
fit_many <- mydata %>%
  group_by(Location) %>%
  nest() %>%
  mutate(., fit = purrr::map(data, ~nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                                                  data = .x,
                                                  iter = 500,
                                                  start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                                                  start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                                  supp_errors = 'Y',
                                                  na.action = na.omit,
                                                  lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))


# run bootstrap over many curves ####
boot_many <- mydata %>%
  group_by(Location) %>%
  # create 200 bootstrap replicates per curve
  do(., boot = modelr::bootstrap(., n = 200, id = 'boot_num')) %>%
  # unnest to show bootstrap number, .id
  unnest_legacy() %>%
  # regroup to include the boot_num
  group_by(., Location, boot_num) %>%
  # run the model using map()
  mutate(fit = map(strap, ~nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 20),
                                         data = data.frame(.),
                                         iter = 50,
                                         start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                                         start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                         lower = c(lnc=-10, E=0, Eh=0, Th=0),
                                         supp_errors = 'Y')
  ))

# new data frame for smooth predictions
new_preds <- mydata %>%
  do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 150), stringsAsFactors = FALSE))

# get max and min for each curve
max_min <- mydata %>%
  group_by(Location) %>%
  summarise(., min_K = min(K), max_K = max(K)) %>%
  ungroup()

# create smoother predictions for unbootstrapped models
preds_many <- fit_many %>%
  unnest_legacy(fit %>% map(augment, newdata = new_preds))

# create smoother predictions for bootstrapped replicates
preds_many <- boot_many %>%
  unnest_legacy(fit %>% map(augment, newdata = new_preds)) %>%
  ungroup() %>%
  # group by each value of K and get quantiles
  group_by(., Location, K) %>%
  summarise(lwr_CI = quantile(.fitted, 0.025),
            upr_CI = quantile(.fitted, 0.975)) %>%
  ungroup() %>%
  merge(., preds_many, by = c('K', 'Location')) %>%
  # merge with max_min to delete predictions outside of the max and min temperatures of each curve
  merge(., max_min, by = c('Location')) %>%
  group_by(., Location) %>%
  filter(., K >= unique(min_K) & K <= unique(max_K)) %>%
  rename(., log.rate = .fitted) %>%
  ungroup()  

# plot predictions 
ggplot(mydata, aes(K - 273.15, log.rate, group = Location)) +
  geom_point(alpha = 0.5, size = 0.5, aes(col = Location)) +
  geom_line(data = preds_many) +
  geom_ribbon(aes(ymin = lwr_CI, ymax = upr_CI, fill = Location), data = preds_many, alpha = .5) + 
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  #scale_color_brewer(palette="Set2")+
  #scale_fill_brewer(palette="Set2") +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  ylab(expression(paste("Respiration Rates (Log " *mu* "mol O"[2], "  "*g^-1 , "  "*hr^-1*")"), sep = " ") )+
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 14, base_family = 'Helvetica') +
  theme(legend.position = c(0.75, 0.15), legend.text=element_text(size=rel(0.8))) +
  ggsave(filename = 'Output/MSPlots/bootstrappedTPC.pdf', device = 'pdf', width = 6, height = 6)



## Make a map of the sites
lats<-c(-0.279722, -1.220667,-0.270361, -1.291444	)	
lons<-c(-90.544861, 	-90.422611, 	-91.435833, -90.496583)
sitetype<-c(" Bartolomé", " Punta Cormorant" ," Punta Espinosa" ," La Botella")
#colors<-c( "#b2182b", "#abd9e9", "#f4a582","#323695")
colors<-cols

pts = data.frame(lats, lons, sitetype,colors)

# Create SpatialPointsDataFrame from this data, assign WGS 84 projection

spdf <- SpatialPointsDataFrame(coords = data.frame(lons, lats), data = data.frame(sitetype),
                               proj4string = CRS("+init=epsg:4326"))

### DOWNLOAD GALAPAGOS DATA -----

URL <- "https://osm2.cartodb.com/api/v2/sql?filename=public.galapagos_islands&q=select+*+from+public.galapagos_islands&format=geojson&bounds=&api_key="
fil <- "gal.json"
if (!file.exists(fil)) download.file(URL, fil)
gal <- readOGR(fil)
#gal <- gSimplify(gUnaryUnion(spTransform(gal, CRS("+init=epsg:31983")), id=NULL), tol=0.001)

# Match projections between spdf and gal

spdf<-spTransform(spdf, CRS("+init=epsg:31983"))

# Verify data lines up

plot(gal); plot(spdf, add=T)

### GGPLOT MAP -----

# Step 1: Create dfs of both survey sites and gal

sites <- data.frame(spdf)
#gal_map <- fortify(gal)

# Step 2: Change format of sitetype field from character to factor for color mapping

spdf$sitetype<-factor(spdf$sitetype, levels = c(" Bartolomè", " Punta Cormorant" ," Punta Espinosa" ," La Botella"))

sites$sitetype<-factor(sites$sitetype, levels = c(" Bartolomè", " Punta Cormorant" ," Punta Espinosa" ," La Botella"))

# Generate ggplot

ecuador <- getData('alt', country='ECU', download = TRUE)
# convert to lat long coords
g_longlat <- spTransform(x = gal, CRSobj = crs(ecuador))
g_longlat<-gSimplify(g_longlat, tol = 0.001)
gal_map<-fortify(g_longlat)

gg<-ggplot()+
  geom_map(map=gal_map, data=gal_map,
           aes(map_id=id),
           color="black", fill="#FFFFFF", size=.5) +
  geom_point(data=pts, aes(x=lons, y=lats, col=sitetype), size = 5)+
  geom_text(data = pts, aes(x = lons, y = lats, label = sitetype), hjust = -0.1, nudge_x = .02, cex = 5)+
  #coord_equal() + 
  coord_map()+
  #theme_map()+
  theme(legend.title=element_text(size=20) , legend.text=element_text(size=16),
        legend.position="top", legend.justification = 'center',
        plot.margin=grid::unit(c(0,0,0,0),"mm")) + 
  labs(color=NULL) +
  theme(legend.position = "none")+
  scale_color_brewer(palette="Set2")+
 # ggsn::north(data=gal_map, symbol = 1, scale = 0.15, location="topright") +
  ggsn::scalebar(data=gal_map, 
                 dist=50, height=0.05, st.size = 3, 
                 location="bottomleft", dist_unit = 'km', transform = TRUE, model = "WGS84") 

# Save ggplot
## take 2 with new map lat/longs

gg<-ggplot()+
  geom_map(map=gal_map, data=gal_map, 
           aes(map_id=id),
           color="black", fill="black", size=.5) +
  geom_point(data=pts, aes(x=lons, y=lats, color = colors), size = 5)+
 # scale_color_manual(values = c("#323695", "#abd9e9", "#b2182b", "#f4a582"))+
  scale_color_manual(values = cols)+
    geom_text_repel(data = pts, aes(x = lons, y = lats, label = sitetype),cex = 5, hjust = -0.2, nudge_x = .05)+
  xlim(-93,-88)+
  ylim(-1.5, 1)+
  # coord_equal() 
  coord_map()+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.title=element_text(size=20) , legend.text=element_text(size=14),
        legend.position="top", legend.justification = 'center',
        plot.margin=grid::unit(c(0,0,0,0),"mm")) + 
  labs(color=NULL) +
  theme(legend.position = "none")+
  ggsn::scalebar(data=gal_map, 
                 dist=50, height=0.05, st.size = 3, 
                 location="bottomleft", dist_unit = 'km', transform = TRUE, model = "WGS84") 

ggsave("Output/MSPlots/Galapagos_map.pdf", gg, width=5, height=5)

