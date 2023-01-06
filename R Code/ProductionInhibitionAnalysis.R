#R code for the manuscript: Temporal variability of microbial response to crude oil exposure in the northern Gulf of Mexico

#Author: Melissa L Brock

#Code for environmental data, production inhibition, regression analysis, and figure building

####Setup R Environment####
getwd()
setwd("/path/to/directory")
getwd()

library(tidyverse)
library(caret)
library(leaps)
library(ggplot2)
library(ggpubr)
library(extrafont)
library(extrafontdb)
library(MASS)
library(car)

#Read in Supplemental Table 1
compiled.data <- read.csv("SupplementalTable1.csv", header = T, row.names = 1, sep = ",")
compiled.data$Date <- as.Date(compiled.data$Date, format='%m/%d/%Y')

####Production Inhibition Plots####
colnames(compiled.data)

#Inhibition of Bacterial Production
leu.error.plot <- ggplot() +
  theme_classic() +
  labs(x ="Date", y = "Inhibition of\nBacterial Production") + 
  geom_line(data = compiled.data, aes(x = Date, y = Positive.regression.slope.1), colour = "black", size = .25, linetype = "dashed") + 
  geom_errorbar(data = compiled.data, aes(x = Date, y = Positive.regression.slope.1, ymin = Positive.regression.slope.1 - Standard.error.of.regression.slope.1, ymax = Positive.regression.slope.1 + Standard.error.of.regression.slope.1), width= 7) + 
  geom_point(data = compiled.data, aes(x = Date, y = Positive.regression.slope.1, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") + 
  scale_y_continuous(breaks = seq(0.05, 0.10, 0.01), labels = c("0.050", "0.060", "0.070", "0.080", "0.090", "0.100")) + 
  expand_limits(y = c(0.05, 0.10)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
leu.error.plot

#Inhibition of Primary Production
bicarb.error.plot <- ggplot() +
  theme_classic() +
  labs(x ="Date", y = "Inhibition of\nPrimary Production") + 
  geom_line(data = compiled.data, aes(x = Date, y = Positive.regression.slope), colour = "black", size = .25, linetype = "dashed") + 
  geom_errorbar(data = compiled.data, aes(x = Date, y = Positive.regression.slope, ymin = Positive.regression.slope - Standard.error.of.regression.slope, ymax = Positive.regression.slope + Standard.error.of.regression.slope), width= 7) + 
  geom_point(data = compiled.data, aes(x = Date, y = Positive.regression.slope, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), legend.title = element_text(size = 14)) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") + 
  scale_y_continuous(breaks = seq(0, 0.15, 0.025)) + 
  expand_limits(y = c(0, 0.15)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
bicarb.error.plot

#####Environmental Plots####
colnames(compiled.data)

#Salinity
sal.plot <- ggplot() +
  theme_classic() +
  labs(x ="Date", y = "Salinity") + 
  geom_line(data = compiled.data, aes(x = Date, y = Salinity), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data, aes(x = Date, y = Salinity, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") + 
  scale_y_continuous(breaks = seq(25, 40, 2.5)) + 
  expand_limits(y = c(25, 40)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
sal.plot

#Temperature
temp.plot <- ggplot() +
  theme_classic() +
  labs(x ="Date", y = expression(paste("Temperature (",degree, "C)"))) + 
  geom_line(data = compiled.data, aes(x = Date, y = Temperature), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data, aes(x = Date, y = Temperature, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") + 
  scale_y_continuous(breaks = seq(10, 30, 5)) + 
  expand_limits(y = c(10, 30)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
temp.plot

#Nitrate+Nitrite
no3.no2.plot <- ggplot() +
  theme_classic() +
  labs(x ="Date", y = expression(paste("Nitrate + Nitrite (",mu, "g N/L)"))) + 
  geom_line(data = compiled.data, aes(x = Date, y = Nitrate.Nitrite), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data, aes(x = Date, y = Nitrate.Nitrite, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_y_continuous(breaks = seq(0, 30, 5)) + 
  expand_limits(y = c(0, 30)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
no3.no2.plot

#Total Kjeldahl Nitrogen (TKN)
tkn.plot <- ggplot() +
  theme_classic() +
  labs(x ="Date", y = "TKN (mg N/L)") + 
  geom_line(data = compiled.data, aes(x = Date, y = Total.Kjeldahl.Nitrogen), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data, aes(x = Date, y = Total.Kjeldahl.Nitrogen, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7")) + 
  expand_limits(y = c(0, .7)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
tkn.plot

#Orthophosphate
phos.plot <- ggplot() +
  theme_classic() +
  labs(x ="Date", y = expression(paste("Orthophosphate (",mu, "g P/L)"))) + 
  geom_line(data = compiled.data, aes(x = Date, y = Orthophosphate), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data, aes(x = Date, y = Orthophosphate, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) + 
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  scale_y_continuous(breaks = seq(10, 20, 2)) + 
  expand_limits(y = c(10, 20)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
phos.plot

####Biological Plots####
colnames(compiled.data)

#Chlorophyll a
chla.plot <- ggplot() +
  theme_classic() +
  labs(x ="Date", y = expression(paste("Chlorophyll a (",mu, "g/L)"))) + 
  geom_line(data = compiled.data, aes(x = Date, y = Chlorophyll.a), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data, aes(x = Date, y = Chlorophyll.a, shape = Season), colour = "black", size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") + 
  scale_y_continuous(breaks = seq(0, 10, 2.5)) + 
  expand_limits(y = c(0, 10)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
chla.plot

#Total Bacterial Counts
bact.counts.plot <- ggplot() +
  theme_classic() +
  labs(x ="Date", y = "Bacterial Counts (cells/mL)") + 
  geom_line(data = compiled.data, aes(x = Date, y = Bacterial.Counts), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data, aes(x = Date, y = Bacterial.Counts, shape = Season), colour = "black", size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") + 
  scale_y_continuous(breaks = seq(50000, 850000, 100000)) + 
  expand_limits(y = c(80000, 880000)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
bact.counts.plot

####Statistics for Production Inhibition####
colnames(compiled.data)
row.names(compiled.data)

###Primary Production
#Remove row with NA measurement
compiled.data.pp <- compiled.data[-24, ]

#Test of normality
shapiro.test(compiled.data.pp$Positive.regression.slope)
#Result: W = 0.94468, p-value = 0.1897

#Test of HOV
leveneTest(y = compiled.data.pp$Positive.regression.slope, group = as.factor(compiled.data.pp$Season))
#Result:     Df F value Pr(>F)
#######group  3  0.9662  0.4272
#############21              

#Build linear model
lm.phyto.inhib <- lm(compiled.data.pp$Positive.regression.slope ~ as.factor(compiled.data.pp$Season))

#Run anova
anova(lm.phyto.inhib)
#Result: p-value = 0.002274

#Post-hoc test
phyto.inhib.aov <- aov(compiled.data.pp$Positive.regression.slope ~ as.factor(compiled.data.pp$Season))
phyto.inhib.tukey <- TukeyHSD(phyto.inhib.aov, conf.level = 0.95)
phyto.inhib.tukey
#Result                diff          lwr          upr     p adj
#Spring-Fall    0.01200874 -0.02686029  0.050877780 0.8245294
#Summer-Fall    0.02452860 -0.01566741  0.064724606 0.3480851
#Winter-Fall   -0.03141126 -0.07028029  0.007457780 0.1418106
#Summer-Spring  0.01251986 -0.02441141  0.049451125 0.7812790
#Winter-Spring -0.04342000 -0.07890241 -0.007937586 0.0129522
#Winter-Summer -0.05593986 -0.09287112 -0.019008590 0.0020041



###Bacterial Production
compiled.data.bp <- compiled.data

#Test of normality
shapiro.test(compiled.data.bp$Positive.regression.slope.1)
#Result: W = 0.93905, p-value = 0.1274

#Test of HOV
leveneTest(y = compiled.data.bp$Positive.regression.slope.1, group = as.factor(compiled.data.bp$Season))
#Result:     Df F value Pr(>F)
#######group  3  2.0274 0.1394
#############22              

#Build linear model
lm.bact.inhib <- lm(compiled.data.bp$Positive.regression.slope.1 ~ as.factor(compiled.data.bp$Season))

#Run anova in base R
anova(lm.bact.inhib)
#Result: p-value = 0.0004394

#Post-hoc test
bact.inhib.aov <- aov(compiled.data.bp$Positive.regression.slope.1 ~ as.factor(compiled.data.bp$Season))
bact.inhib.tukey <- TukeyHSD(bact.inhib.aov, conf.level = 0.95)
bact.inhib.tukey
#Result                diff          lwr          upr     p adj
#Spring-Fall   -0.006455629 -0.022987246  0.010075989 0.7024544
#Summer-Fall    0.002321229 -0.014210389  0.018852846 0.9793613
#Winter-Fall   -0.023674200 -0.040205817 -0.007142583 0.0033180
#Summer-Spring  0.008776857 -0.006314376  0.023868090 0.3911389
#Winter-Spring -0.017218571 -0.032309804 -0.002127338 0.0214508
#Winter-Summer -0.025995429 -0.041086662 -0.010904196 0.0004830

####Statistics for Environmental Parameters####

###Salinity
range(compiled.data$Salinity) #27 37
lm.sal <- lm(compiled.data$Salinity ~ as.factor(compiled.data$Season))
anova(lm.sal) #p-value = 0.1779


###Temperature
range(compiled.data$Temperature) #13.6 28.2
lm.temp <- lm(compiled.data$Temperature ~ as.factor(compiled.data$Season))
anova(lm.temp) #p-value < 0.001
temp.aov <- aov(compiled.data$Temperature ~ as.factor(compiled.data$Season))
temp.tukey <- TukeyHSD(temp.aov, conf.level = 0.95)
temp.tukey
#                diff        lwr       upr     p adj
#Spring-Fall    -1.691429  -5.200185  1.817328 0.5492552
#Summer-Fall     4.051429   0.542672  7.560185 0.0197068
#Winter-Fall    -6.491429 -10.000185 -2.982672 0.0002073
#Summer-Spring   5.742857   2.539815  8.945899 0.0003026
#Winter-Spring  -4.800000  -8.003042 -1.596958 0.0021391
#Winter-Summer -10.542857 -13.745899 -7.339815 0.0000000


###Nitrate + Nitrite
range(compiled.data$Nitrate.Nitrite) #8.46 26.13
lm.No3No2 <- lm(compiled.data$Nitrate.Nitrite ~ as.factor(compiled.data$Season))
anova(lm.No3No2) #p-value = 0.257

###TKN
range(compiled.data$Total.Kjeldahl.Nitrogen) #0.257 0.630
lm.tkn <- lm(compiled.data$Total.Kjeldahl.Nitrogen ~ as.factor(compiled.data$Season))
anova(lm.tkn) #p-value = 0.2844

###Orthophosphate
range(compiled.data$Orthophosphate) #10.72 17.37
lm.ortho <- lm(compiled.data$Orthophosphate ~ as.factor(compiled.data$Season))
anova(lm.ortho) #p-value = 0.2973


####Statistics for Biological Parameters####
###Chlorophyll a
range(compiled.data$Chlorophyll.a) #0.09 9.53
lm.chla <- lm(compiled.data$Chlorophyll.a ~ as.factor(compiled.data$Season))
anova(lm.chla) #p-value = 0.008116
chla.aov <- aov(compiled.data$Chlorophyll.a ~ as.factor(compiled.data$Season))
chla.tukey <- TukeyHSD(chla.aov, conf.level = 0.95)
chla.tukey
#                    diff        lwr      upr     p adj
#Spring-Fall    0.5148571 -3.2837860 4.313500 0.9813563
#Summer-Fall   -1.5251429 -5.3237860 2.273500 0.6844147
#Winter-Fall    3.2577143 -0.5409289 7.056357 0.1103698
#Summer-Spring -2.0400000 -5.5076709 1.427671 0.3813338
#Winter-Spring  2.7428571 -0.7248138 6.210528 0.1554169
#Winter-Summer  4.7828571  1.3151862 8.250528 0.0046901


###Bacterial Counts
range(compiled.data$Bacterial.Counts) #89515 793763
lm.bcounts <- lm(compiled.data$Bacterial.Counts ~ as.factor(compiled.data$Season))
anova(lm.bcounts) #p-value = 0.002452
bcounts.aov <- aov(compiled.data$Bacterial.Counts ~ as.factor(compiled.data$Season))
bcounts.tukey <- TukeyHSD(bcounts.aov, conf.level = 0.95)
bcounts.tukey
#                    diff        lwr      upr     p adj
#Spring-Fall    180489.83  -79244.12 440223.78 0.2449709
#Summer-Fall    370993.11  111259.17 630727.06 0.0033996
#Winter-Fall     66614.83 -193119.12 326348.78 0.8911911
#Summer-Spring  190503.29  -46600.28 427606.86 0.1459911
#Winter-Spring -113875.00 -350978.57 123228.57 0.5522399
#Winter-Summer -304378.29 -541481.86 -67274.72 0.0087133

####Regression Analysis####
#Primary Production Inhibition
colnames(compiled.data)
for.mlr.bicarb <- compiled.data[-24, -c(1:7, 9:11, 18, 19)] #remove row with NA value and unnecessary columns
colnames(for.mlr.bicarb)

model1.bicarb <- lm(Positive.regression.slope ~ Temperature, data = for.mlr.bicarb)
summary(model1.bicarb) #p = 0.000164; Aj. R2 = 0.4444
plot(y = for.mlr.bicarb$Positive.regression.slope, x = for.mlr.bicarb$Temperature)

#model2.bicarb <- lm(Positive.regression.slope ~ Salinity, data = for.mlr.bicarb)
#summary(model2.bicarb)
#plot(y = for.mlr.bicarb$Positive.regression.slope, x = for.mlr.bicarb$Salinity) #not a linear relationship

#model3.bicarb <- lm(Positive.regression.slope ~ Nitrate.Nitrite, data = for.mlr.bicarb)
#summary(model3.bicarb) #p = 0.07392; Aj. R2 = 0.09454
#plot(y = for.mlr.bicarb$Positive.regression.slope, x = for.mlr.bicarb$Nitrate.Nitrite)

model4.bicarb <- lm(Positive.regression.slope ~ Total.Kjeldahl.Nitrogen, data = for.mlr.bicarb)
summary(model4.bicarb) #p = 0.001529; Adj. R2 = 0.3319
plot(y = for.mlr.bicarb$Positive.regression.slope, x = for.mlr.bicarb$Total.Kjeldahl.Nitrogen)

#model5.bicarb <- lm(Positive.regression.slope ~ Orthophosphate, data = for.mlr.bicarb)
#summary(model5.bicarb) #p = 0.4875; Adj. R2 = -0.02136
#plot(y = for.mlr.bicarb$Positive.regression.slope, x = for.mlr.bicarb$Orthophosphate)

model6.bicarb <- lm(Positive.regression.slope ~ Chlorophyll.a, data = for.mlr.bicarb)
summary(model6.bicarb) #p < 0.0001; Adj. R2 = 0.7653
plot(y = for.mlr.bicarb$Positive.regression.slope, x = for.mlr.bicarb$Chlorophyll.a)



#Bacterial Production Inhibition
colnames(compiled.data)
for.mlr.leucine <- compiled.data[, -c(1:9, 11, 17)]
colnames(for.mlr.leucine)

model1.leucine <- lm(Positive.regression.slope.1 ~ Temperature, data = for.mlr.leucine)
summary(model1.leucine) #p < 0.0001; Adj. R2 = 0.6363
plot(y = for.mlr.leucine$Positive.regression.slope.1, x = for.mlr.leucine$Temperature)

#model2.leucine <- lm(Positive.regression.slope.1 ~ Salinity, data = for.mlr.leucine)
#summary(model2.leucine) #not a linear relationship

#model3.leucine <- lm(Positive.regression.slope.1 ~ Nitrate.Nitrite, data = for.mlr.leucine)
#summary(model3.leucine) #p = 0.3391; Adj. R2 = -0.00195
#plot(y = for.mlr.leucine$Positive.regression.slope.1, x = for.mlr.leucine$Nitrate.Nitrite)

model4.leucine <- lm(Positive.regression.slope.1 ~ Total.Kjeldahl.Nitrogen, data = for.mlr.leucine)
summary(model4.leucine) #p = 0.03221; Adj. R2 = 0.143
plot(y = for.mlr.leucine$Positive.regression.slope.1, x = for.mlr.leucine$Total.Kjeldahl.Nitrogen)

#model5.leucine <- lm(Positive.regression.slope.1 ~ Orthophosphate, data = for.mlr.leucine)
#summary(model5.leucine) #p = 0.6698; Adj. R2 = -0.03364
#plot(y = for.mlr.leucine$Positive.regression.slope.1, x = for.mlr.leucine$Orthophosphate)

#model6.leucine <- lm(Positive.regression.slope.1 ~ Bacterial.Counts, data = for.mlr.leucine)
#summary(model6.leucine) #p = 0.3445; Adj. R2 = -0.002814
#plot(y = for.mlr.leucine$Positive.regression.slope.1, x = for.mlr.leucine$Bacterial.Counts)

model7a.leucine <- lm(Positive.regression.slope.1 ~ Alpha.diversity_Shannon, data = for.mlr.leucine)
summary(model7a.leucine) #p = 0.02807; Adj. R2 = 0.1515
plot(y = for.mlr.leucine$Positive.regression.slope.1, x = for.mlr.leucine$Alpha.diversity_Shannon)

#model7b.leucine <- lm(Positive.regression.slope.1 ~ Alpha.diversity_Richness, data = for.mlr.leucine)
#summary(model7b.leucine) #p-value: 0.06723; Adj. R2 =  0.09665
#plot(y = for.mlr.leucine$Positive.regression.slope.1, x = for.mlr.leucine$Alpha.diversity_Richness)

model7c.leucine <- lm(Positive.regression.slope.1 ~ Alpha.diversity_evenness, data = for.mlr.leucine)
summary(model7c.leucine) #p = 0.01177; Adj. R2 = 0.2047
plot(y = for.mlr.leucine$Positive.regression.slope.1, x = for.mlr.leucine$Alpha.diversity_evenness)

####Regression Plots####
###Phytoplankton Production Inhibition - Chlorophyll a, Temperature, TKN
colnames(compiled.data.pp)

#Chlorophyll a: p < 0.0001; Adj. R2 = 0.7653
lm.chla.phyto <- ggplot() + 
  theme_classic() +
  labs(y ="Inhibition of\nPrimary Production", x = expression(paste("Chlorophyll a (",mu, "g/L)"))) + 
  annotate("text", x = 7, y = 0.12, label = expression('Adj. R'^"2" == 0.765), size = 5, hjust = 0, family = "Times New Roman") + 
  annotate("text", x = 7, y = 0.11, label = expression('p < 0.0001'), size = 5, hjust = 0, family = "Times New Roman") +
  geom_smooth(method = "lm", data = compiled.data.pp, aes(x = Chlorophyll.a, y = Positive.regression.slope), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data.pp, aes(y = Positive.regression.slope, x = Chlorophyll.a, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), legend.title = element_text(size = 14)) + 
  scale_y_continuous(breaks = seq(0, 0.125, 0.025)) + 
  scale_x_continuous(breaks = seq(0, 10, 1)) + 
  expand_limits(x = c(0, 10), y = c(0, 0.125)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
lm.chla.phyto

#Temperature: p = 0.000164; Aj. R2 = 0.4444
lm.temp.phyto <- ggplot() + 
  theme_classic() +
  labs(y = "Inhibition of\nPrimary Production", x = expression(paste("Temperature (",degree, "C)"))) + 
  annotate("text", x = 14, y = .12, label = expression('Adj. R'^"2" ==0.444), size = 5, hjust = 0, family = "Times New Roman") + 
  annotate("text", x = 14, y = .11, label = expression('p = 0.0002'), size = 5, hjust = 0, family = "Times New Roman") +
  geom_smooth(method = "lm", data = compiled.data.pp, aes(x = Temperature, y = Positive.regression.slope), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data.pp, aes(y = Positive.regression.slope, x = Temperature, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), legend.title = element_text(size = 14)) + 
  scale_x_continuous(breaks = seq(14, 30, 2)) + 
  scale_y_continuous(breaks = seq(0, 0.125, 0.025)) + 
  expand_limits(x = c(14, 29), y = c(0, 0.125)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
lm.temp.phyto

#TKN: p = 0.001529; Adj. R2 = 0.3319
lm.tkn.phyto <- ggplot() + 
  theme_classic() +
  labs(y ="Inhibition of\nPrimary Production", x = "TKN (mg/L)") + 
  annotate("text", x = 0.5, y = 0.12, label = expression('Adj. R'^"2" == 0.332), size = 5, hjust = 0, family = "Times New Roman") + 
  annotate("text", x = 0.5, y = 0.11, label = expression('p = 0.002'), size = 5, hjust = 0, family = "Times New Roman") +
  geom_smooth(method = "lm", data = compiled.data.pp, aes(x = Total.Kjeldahl.Nitrogen, y = Positive.regression.slope), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data.pp, aes(y = Positive.regression.slope, x = Total.Kjeldahl.Nitrogen, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), legend.title = element_text(size = 14)) + 
  scale_x_continuous(breaks = seq(0.2, 0.65, 0.05)) + 
  scale_y_continuous(breaks = seq(-0.025, 0.125, 0.025)) + 
  expand_limits(x = c(0.25, .65),  y = c(-0.025, 0.125)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
lm.tkn.phyto



###Bacterial Production Inhibition - Alpha-diversity (Shannon Index & Pielou's evenness), Temperature, TKN
colnames(compiled.data.bp)

#Alpha-diversity (Shannon Index): p = 0.02807; Adj. R2 = 0.1515
lm.div.bact <- ggplot() + 
  theme_classic() +
  labs(y ="Inhibition of\nBacterial Production", x = "Alpha-diversity (Shannon Index)") + 
  annotate("text", x = 5.4, y = 0.09, label = expression('Adj. R'^"2" ==0.152), size = 5, hjust = 0, family = "Times New Roman") + 
  annotate("text", x = 5.4, y = 0.085, label = expression('p = 0.028'), size = 5, hjust = 0, family = "Times New Roman") +
  geom_smooth(method = "lm", data = compiled.data.bp, aes(x = Alpha.diversity_Shannon, y = Positive.regression.slope.1), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data.bp, aes(y = Positive.regression.slope.1, x = Alpha.diversity_Shannon, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) + 
  scale_x_continuous(breaks = seq(4, 6, 0.5)) + 
  scale_y_continuous(breaks = seq(0.04, 0.1, 0.01), labels = c("0.040", "0.050", "0.060", "0.070", "0.080", "0.090", "0.100")) + 
  expand_limits(y = c(0.04, 0.1)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
lm.div.bact

#Alpha-diversity (Pielou's evenness): p = 0.01177; Adj. R2 = 0.2047
lm.even.bact <- ggplot() + 
  theme_classic() +
  labs(y ="Inhibition of\nBacterial Production", x = "Pielou's evenness") + 
  annotate("text", x = 0.85, y = 0.09, label = expression('Adj. R'^"2" ==0.205), size = 5, hjust = 0, family = "Times New Roman") + 
  annotate("text", x = 0.85, y = 0.085, label = expression('p = 0.012'), size = 5, hjust = 0, family = "Times New Roman") +
  geom_smooth(method = "lm", data = compiled.data.bp, aes(x = Alpha.diversity_evenness, y = Positive.regression.slope.1), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data.bp, aes(y = Positive.regression.slope.1, x = Alpha.diversity_evenness, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) + 
  scale_x_continuous(breaks = seq(0.75, 0.9, 0.05)) + 
  scale_y_continuous(breaks = seq(0.04, 0.1, 0.01), labels = c("0.040", "0.050", "0.060", "0.070", "0.080", "0.090", "0.100")) + 
  expand_limits(y = c(0.04, 0.1), x = c(0.75, 0.90)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
lm.even.bact

#Temperature: p < 0.0001; Adj. R2 = 0.6363
lm.temp.bact <- ggplot() + 
  theme_classic() +
  labs(y ="Inhibition of\nBacterial Production", x = expression(paste("Temperature (",degree, "C)"))) + 
  annotate("text", x = 12, y = 0.09, label = expression('Adj. R'^"2" ==0.636), size = 5, hjust = 0, family = "Times New Roman") + 
  annotate("text", x = 12, y = 0.085, label = expression('p < 0.0001'), size = 5, hjust = 0, family = "Times New Roman") +
  geom_smooth(method = "lm", data = compiled.data.bp, aes(x = Temperature, y = Positive.regression.slope.1), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data.bp, aes(y = Positive.regression.slope.1, x = Temperature, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) + 
  scale_x_continuous(breaks = seq(12, 30, 2)) + 
  scale_y_continuous(breaks = seq(0.04, 0.1, 0.01), labels = c("0.040", "0.050", "0.060", "0.070", "0.080", "0.090", "0.100")) + 
  expand_limits(x = c(12, 30), y = c(0.04, 0.1)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
lm.temp.bact

#TKN: p = 0.03221; Adj. R2 = 0.143 
lm.tkn.bact <- ggplot() + 
  theme_classic() +
  labs(y ="Inhibition of\nBacterial Production", x = "TKN (mg N/L)") + 
  annotate("text", x = 0.5, y = 0.09, label = expression('Adj. R'^"2" == 0.143), size = 5, hjust = 0, family = "Times New Roman") + 
  annotate("text", x = 0.5, y = 0.085, label = expression('p = 0.032'), size = 5, hjust = 0, family = "Times New Roman") +
  geom_smooth(method = "lm", data = compiled.data.bp, aes(x = Total.Kjeldahl.Nitrogen, y = Positive.regression.slope.1), colour = "black", size = .25, linetype = "dashed") + 
  geom_point(data = compiled.data.bp, aes(y = Positive.regression.slope.1, x = Total.Kjeldahl.Nitrogen, shape = Season, colour = Season), size = 2.5) + 
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), text = element_text(family = "Times New Roman"), legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) + 
  scale_x_continuous(breaks = seq(0.2, 0.65, 0.05)) + 
  scale_y_continuous(breaks = seq(0.04, 0.1, 0.01), labels = c("0.040", "0.050", "0.060", "0.070", "0.080", "0.090", "0.100")) + 
  expand_limits(x = c(0.2, 0.65),  y = c(0.04, 0.1)) + 
  scale_shape_manual(values=c(15, 16, 17, 8), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter")) + 
  scale_color_manual(values = c("aquamarine4", "firebrick3", "darkorange2", "blue4"), name = "Season:", breaks = c("Spring", "Summer", "Fall", "Winter"))
lm.tkn.bact

####Figures####
###Figure 1
png("Frontiers-Figure1-138.png", units="in", width=8, height=10, res=600)

ggarrange(temp.plot, sal.plot, no3.no2.plot, tkn.plot, phos.plot, nrow = 3, ncol = 2, align = "hv", labels = c("(A)", "(B)", "(C)", "(D)", "(E)"), font.label = list(size = 16, color = "black", face = "bold", family = "Times New Roman"), vjust = 0.4, common.legend = T, legend = "top")

dev.off()



###Figures 2 and 3 are in "Microbial Community Analysis"



###Figure 4
png("Frontiers-Figure4-138.png", units="in", width=10, height=7.5, res=600)

ggarrange(bicarb.error.plot, lm.chla.phyto, lm.temp.phyto, lm.tkn.phyto, ncol = 2, nrow = 2, labels = c("(A)", "(B)", "(C)", "(D)"), font.label = list(size = 16, color = "black", face = "bold", family = "Times New Roman"), common.legend = T, legend = "top")

dev.off()



###Figure 5
png("Frontiers-Figure5-138.png", units="in", width=10, height=7.5, res=600)

ggarrange(leu.error.plot, lm.even.bact, lm.temp.bact, lm.tkn.bact, labels = c("(A)", "(B)", "(C)", "(D)"), nrow = 2, ncol = 2, font.label = list(size = 16, color = "black", face = "bold", family = "Times New Roman"), common.legend = T, legend = "top")

dev.off()

####Figure Accessibility####
#Obtained RGB codes for colors used in Figures 1, 4, and 5
col2rgb(c("aquamarine4", "firebrick3", "darkorange2", "blue4"), alpha = FALSE)

#Checked compatibility at: https://color.adobe.com/create/color-accessibility
#No conflicts found. Swatches are color blind safe.