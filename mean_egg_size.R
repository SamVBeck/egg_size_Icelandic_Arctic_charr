rm(list = ls())

library(ggplot2)
library(dplyr)
library(plyr)
library(dplyr)
library(raster)
library(RColorBrewer)
library(AICcmodavg)

# load data
mydata <- read.csv("~/Documents/PhD/R_2014:15/R_postPhD/eggSize/size_data.csv")
mydata <- mydata[!apply(is.na(mydata) | mydata == "", 1, all), ] # remove empty columns at bottom
mydata$date_processed <- as.Date(mydata$date_processed, format = "%m/%d/%y")
mydata$ID <- as.factor(mydata$ID)
mydata$size <- as.numeric(as.character(mydata$size)) # will get a warning msg as lots of NAs
mydata$year <- as.factor(mydata$year)
mydata$male <- as.factor(mydata$male)
mydata$divergence <- as.factor(mydata$divergence)
mydata$pop.morph <- factor(mydata$pop.morph, levels = c("FJ","VS","VB","SV","TP","TLB","GB"))

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

mydata <- completeFun(mydata, "size")

### subset by: 1) dev.stage; 2) Population; 3) pop.morph; 4) age
PF <- subset(mydata, dev_stage=="PF",select=lake:notes)
E<- subset(mydata, dev_stage=="E",select=lake:notes)
H <- subset(mydata, dev_stage=="H",select=lake:notes)
FF <- subset(mydata, dev_stage=="FF",select=lake:notes)

FJ <- subset(mydata, pop.morph=="FJ",select=lake:notes)
VB <- subset(mydata, pop.morph=="VB",select=lake:notes)
VS <- subset(mydata, pop.morph=="VS",select=lake:notes)
SV <- subset(mydata, pop.morph=="SV",select=lake:notes)
TLB <- subset(mydata, pop.morph=="TLB",select=lake:notes)
TP <- subset(mydata, pop.morph=="TP",select=lake:notes)
GB <- subset(mydata, pop.morph=="GB",select=lake:notes)

EggSizes <- subset(mydata, dev_stage == "PF" | dev_stage == "E")
EggSizes <- droplevels(EggSizes)



## dropped TLB females 5 and 7 as sample size was too small at H stage
H$female <- as.factor(H$female)
H <- H[ ! H$female %in% 5, ]
H <- H[ ! H$female %in% 7, ]
H <- droplevels(H)


# female used for DD timings
DDfemales <- c("292","294","295","298","299","308","309","311","312","313","314","315","338","339","343","1","2","3","5","6","7","8","9","10","11","116","119","124","128","83","84","85","86","92","240","241","250","251","252","253","16.17","18","23","24","25","26","28","30","31","54","56","65","67","68","69")
mydata_DDfemale <- mydata[mydata$female %in% DDfemales,]

PF_DD <- subset(mydata_DDfemale, dev_stage=="PF",select=lake:notes)
E_DD<- subset(mydata_DDfemale, dev_stage=="E",select=lake:notes)
H_DD <- subset(mydata_DDfemale, dev_stage=="H",select=lake:notes)
FF_DD <- subset(mydata_DDfemale, dev_stage=="FF",select=lake:notes)

## dropped TLB females 5 and 7, as removed these from H stage
H_DD$female <- as.factor(H_DD$female)
H_DD <- H_DD[ ! H_DD$female %in% 5, ]
H_DD <- H_DD[ ! H_DD$female %in% 7, ]
H_DD <- droplevels(H_DD)


# get colours
mypalette <- brewer.pal(n = 12, name = 'Paired')
#insert the 6th color into the correct place on the list
mypalette <- append(mypalette, mypalette[6], 0)
mypalette <- append(mypalette, mypalette[8], 3)



######## now do some absolute egg size variation plots:
p1 <- ggplot(PF, aes(pop.morph, size, fill = pop.morph, group=female))+
  geom_boxplot()+
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  theme(axis.text.x=element_text(size = 30, colour = "black"), 
        axis.title.x=element_text(size = 30, colour = "black"),
        axis.text.y=element_text(size = 30), axis.title.y=element_text(size = 30),
        legend.position = "bottomright")+
  theme_bw()
plotPop <- p1 + labs(y = "Egg size (mm)", x = element_blank()) + 
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# now add the alpha boxplot overlay
p1 <- plotPop + geom_boxplot(aes(fill = pop.morph, group = pop.morph), alpha = 0.25)#+


### for female FL and egg size, use means (see below)
####################################################################################
##### ALL ANALYSES ON MEAN DATA
####################################################################################
#####    Female stats = femaleData only
#####    Egg size stats = EggSizes (2 observations per female: 1) mean PF egg size; and 2) mean E egg size)
#####    Dev timings = mydata
##################################################################################
rm(list = ls())
meanData <- read.csv("~/Documents/PhD/R_2014:15/R_postPhD/eggSize/female_data.csv")
myEggSize <- read.csv("~/Documents/PhD/R_2014:15/R_postPhD/eggSize/egg_size.csv")
myEggSize$year <- as.factor(myEggSize$year)
myEggSize$male <- as.factor(myEggSize$male)
myEggSize$pop.morph <- factor(myEggSize$pop.morph, levels = c("FJ","VS","VB","SV","TP","TLB","GB"))

library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(plyr)
library(lmerTest)
library(emmeans)
library(car)
library(AICcmodavg)
library(multcomp)

# get colours
mypalette <- brewer.pal(n = 12, name = 'Paired')
#insert the 6th color into the correct place on the list
mypalette <- append(mypalette, mypalette[6], 0)
mypalette <- append(mypalette, mypalette[8], 3)



### subset by: 1) dev.stage; 2) pop.morph
PF <- subset(myEggSize, dev_stage=="PF",select=lake:eggCV)
E<- subset(myEggSize, dev_stage=="E",select=lake:eggCV)
H <- subset(myEggSize, dev_stage=="H",select=lake:eggCV)
FF <- subset(myEggSize, dev_stage=="FF",select=lake:eggCV)

FJ <- subset(myEggSize, pop.morph=="FJ",select=lake:eggCV)
VB <- subset(myEggSize, pop.morph=="VB",select=lake:eggCV)
VS <- subset(myEggSize, pop.morph=="VS",select=lake:eggCV)
SV <- subset(myEggSize, pop.morph=="SV",select=lake:eggCV)
TLB <- subset(myEggSize, pop.morph=="TLB",select=lake:eggCV)
TP <- subset(myEggSize, pop.morph=="TP",select=lake:eggCV)
GB <- subset(myEggSize, pop.morph=="GB",select=lake:eggCV)

FJ <- droplevels(FJ); VB <- droplevels(VB); VS <- droplevels(VS)
SV <- droplevels(SV); TLB <- droplevels(TLB); TP <- droplevels(TP)
GB <- droplevels(GB)

# for analyses on females
femaleData <- myEggSize %>% 
  distinct(female, .keep_all = T)
femaleData <- droplevels(femaleData)



setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")


## female FL, morph and age are all confounded. Do model selection with all to determine best result
m1 <- lm(femaleFL ~ pop.morph, data = femaleData)
m2<- lm(femaleFL ~ age, data = femaleData)

m1 <- lm(age ~ pop.morph, data = femaleData)
m2 <-lm(age ~ femaleFL, data = femaleData)

models <- list(m1,m2)
model.names <- c("pop.morph","age")  # change accordingly
model.names <- c("pop.morph","femaleFL") # change accordingly
aictab(cand.set = models, modnames = model.names)

Anova(m1)
summary(m1)

# plot
p1 <- ggplot(femaleData, aes(x = pop.morph, y = femaleFL))+ 
  geom_boxplot(aes(fill = pop.morph, group = pop.morph)) +
  xlab("Morph") +
  ylab("Female FL (cm)")
p2 <- p1 + 
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 80, colour = "black"), axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80, colour = "black"), axis.title.y=element_text(size = 80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")

# plot LSM
leastsquare <- lsmeans(m1, ~pop.morph, adjust = "tukey") # compare morphs
leastsquare <- lsmeans(m1, ~age, adjust = "tukey") # compare age
pairs(leastsquare)
test(leastsquare)  # these show differences after accounting for covariates (i.e. female FL), i.e. they are comparable as they use  predictions from uniform female FL
pwpm(leastsquare, flip = T, reverse = T) 
p1 <- plot(leastsquare, comparisons = TRUE, CIs = T, adjust = "tukey")+
  theme_bw() 
plotPop <- p1 + labs(y = "Morph", x = "Age") + 
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1 <- plotPop + coord_flip()







#m1 <- lm(meanPFandE ~ femaleFL, data = femaleData)
#m2<- lm(meanPFandE ~ age, data = femaleData)
m1 <- lm(meanPFandE ~ pop.morph, data = femaleData)
#models <- list(m1, m2, m3)
#model.names <- c("FL","age","morph")
#aictab(cand.set = models, modnames = model.names)
# to get adjusted r2
#m2 <- lm(meanPFandE ~ femaleFL, data = femaleData)
#m3 <- lm(meanPFandE ~ age, data = femaleData)
#summary(m3)
Anova(m1)
summary(m1)
################ MODEL VALIDATION  ################
plot(m1, add.smooth = FALSE, which = 1) # fitted vs. residuals (homogeneity)
mod1 <- resid(m1)
hist(mod1, xlab = "Residuals", main = "") # histogram of residuals (normality)

## least square means
leastsquare <- lsmeans(m1, ~pop.morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, flip = T, reverse = T) 
p1 <- plot(leastsquare, comparisons = TRUE, CIs = T, adjust = "tukey")+
  theme_bw() 
plotPop <- p1 + labs(y = "Egg size (mm)", x = "Morph") + 
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1 <- plotPop + coord_flip()
p1










### how does egg size change with female FL 
m1 <- lm(meanPFandE ~ femaleFL, data = femaleData)
Anova(m1)

################ MODEL VALIDATION  ################
plot(m1, add.smooth = FALSE, which = 1) # fitted vs. residuals (homogeneity)
mod1 <- resid(m1)
hist(mod1, xlab = "Residuals", main = "") # histogram of residuals (normality)

## slightly skewed, transform? ## CHECK! do I need to log transform to adjust skewedness? If so, do I transform residuals or the actual variable (i.e. fl)?
femaleData$logFemaleFL <- log10(femaleData$femaleFL)


# plot how egg size changes with female FL
p1 <- ggplot(femaleData, aes(x = femaleFL, y = meanPFandE, shape = pop.morph, fill = pop.morph)) +
  geom_point(size = 19) +
  scale_shape_manual(values = c(21:25,22,21)) +
  xlab("Female FL (cm)") +
  ylab("Egg size (mm)")
p2 <- p1 + 
  geom_smooth(aes(group = pop.morph, col = pop.morph, fill = pop.morph), method = "lm", level = 0.95) + 
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 80, colour = "black"), axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80, colour = "black"), axis.title.y=element_text(size = 80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")




# plot how egg size changes with female age
p1 <- ggplot(femaleData, aes(x = age, y = meanPFandE, shape = pop.morph, fill = pop.morph)) +
  geom_point(size = 19) +
  scale_shape_manual(values = c(21:25,22,21)) +
  xlab("Age") +
  ylab("Egg size (mm)")
p2 <- p1 + 
  geom_smooth(aes(group = pop.morph, col = pop.morph, fill = pop.morph), method = "lm", level = 0.95) + 
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  theme_bw() +
  theme(axis.text.x=element_text(size = 80, colour = "black"), axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80, colour = "black"), axis.title.y=element_text(size = 80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")


## try and remove the outlier
noOutlier_TLB <- femaleData[ ! femaleData$female %in% 115, ]
noTLB <- femaleData[ ! femaleData$pop.morph %in% "TLB", ]
m1 <- lm(meanPFandE ~ pop.morph, data = noOutlier_TLB)
Anova(m1) # still significantly correlated


### what is the correlation between age and female FL
m1 <- lm(age ~ femaleFL, data = femaleData)
m1 <- lm(meanPFandE ~ age, data = noTLB)
Anova(m1)
summary(m1) 






##### 1) does female FL differ with female age and morph?
m1 <- lm(femaleFL ~ age + pop.morph, data = femaleData)
Anova(m1)
summary(m1)
################ MODEL VALIDATION  ################
plot(m1, add.smooth = FALSE, which = 1) # fitted vs. residuals (homogeneity)
mod1 <- resid(m1)
hist(mod1, xlab = "Residuals", main = "") # histogram of residuals (normality)
###################################################
leastsquare <- lsmeans(m1, ~pop.morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, flip = T, reverse = T) # for some reason it didn't work wel with previous settings!
p1 <- plot(leastsquare, comparisons = TRUE, CIs = T, adjust = "tukey")+
  theme_bw() 
plotPop <- p1 + labs(y = "Morph", x = "Female FL (cm)") + 
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1 <- plotPop + coord_flip()
p1

# The blue bars are confidence intervals for the EMMs, and the red arrows are for the comparisons among them. 
# If an arrow from one mean overlaps an arrow from another group, the difference is not significant.






### QUESTION 1: Are there differences in egg size among populations/morphs?
  # 1a) How do morphs differ in egg size? (absolute)
m1 <- lm(meanPFandE ~ pop.morph, data = femaleData)
Anova(m1)
summary(m1)
nobs(m1)
################ MODEL VALIDATION  ################
plot(m1, add.smooth = FALSE, which = 1) # fitted vs. residuals (homogeneity)
mod1 <- resid(m1)
hist(mod1, xlab = "Residuals", main = "") # histogram of residuals (normality)
###################################################
leastsquare <- lsmeans(m1, ~pop.morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, means = T, flip = T,     # args for pwpm()
     reverse = T, sort = F,                          # args for pairs()
     side = ">", delta = 0.05, adjust = "tukey")  # args for test()
p1 <- plot(leastsquare, comparisons = T, CIs = T, adjust = "tukey")+
  theme_bw() 
plotPop <- p1 + labs(y = "Morph", x = "Egg size (mm)") + 
  coord_flip() +
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plotPop








## 2) do morphs differ in egg size variability?
m1 <- lm(eggCV ~ pop.morph, data = femaleData)
m2 <- lm(eggCV ~ femaleFL, data = femaleData) 
m3 <- lm(eggCV ~ age, data = femaleData) 
models <- list(m1, m2, m3)
model.names <- c("morph","FL", "age")
aictab(cand.set = models, modnames = model.names)
Anova(m1)
summary(m3)
################ MODEL VALIDATION  ################
plot(m1, add.smooth = FALSE, which = 1) # fitted vs. residuals (homogeneity)
mod1 <- resid(m1)
hist(mod1, xlab = "Residuals", main = "") # histogram of residuals (normality)