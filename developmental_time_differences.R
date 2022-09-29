######### now testing for developmental time differences
rm(list = ls())
mydata <- read.csv("~/Documents/PhD/R_2014:15/R_postPhD/eggSize/developmental_timings.csv")
mydata <- mydata[!apply(is.na(mydata) | mydata == "", 1, all), ] # remove empty columns at bottom
mydata$year <- as.factor(mydata$year)
mydata$morph <- factor(mydata$morph, levels = c("FJ","VS","VB","TP","TLB"))
mydata$female <- as.factor(mydata$female)
mydata$age <- as.factor(mydata$age)
levels(mydata$morph)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(plyr)

str(mydata)

FJ <- subset(mydata, morph=="FJ",select=morph:eggSize)
VB <- subset(mydata, morph=="VB",select=morph:eggSize)
VS <- subset(mydata, morph=="VS",select=morph:eggSize)
TLB <- subset(mydata, morph=="TLB",select=morph:eggSize)
TP <- subset(mydata, morph=="TP",select=morph:eggSize)

setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")


myshapes <- c("FJ"=21,  "VS"=23, "VB"=22,"TP"=21,"TLB"=24)

# DD are only between FJ, VS, VB, TP and TLB, so need to re-do colours accordingly 
library(RColorBrewer)
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n = 12, name = 'Paired')
mypalette <- brewer.pal(n = 12, name = 'Paired')
#insert the 6th color into the correct place on the list
mypalette <- append(mypalette, mypalette[6], 0)


library(lmerTest)
library(multcompView)
library(emmeans)
library(car)
library(AICcmodavg)

mydata$dev_stage

mydata %>%
  group_by(pop.morph,dev_stage) %>%
  dplyr::summarise(mean = mean(DD), sd = sd(DD),  n = n()) %>%
  print(n = 100)

mydata %>%
  group_by(HDD) %>%
  dplyr::summarise(n = n())

### first plot non-egg size adjusted patterns:
m1 <- lm(HDD ~ morph, data = mydata)
m2 <- lm(HDD ~ eggSize, data = mydata) 
models <- list(m1, m2)
model.names <- c("morph","eggSize")
aictab(cand.set = models, modnames = model.names)

Anova(m1)
summary(m1)
leastsquare <- lsmeans(m1, ~morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, flip = T, reverse = T) 


#m1 <- lm(HDD ~ morph * eggSize, data = mydata) # no female FL as female means
m1 <- lm(HDD ~ morph + eggSize, data = mydata) # NS interaction, therefore dropped
Anova(m1)
summary(m1)
leastsquare <- lsmeans(m1, ~morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, flip = T, reverse = T) 
p1 <- plot(leastsquare, comparisons = TRUE, CIs = T, adjust = "tukey")+
  theme_bw() 
plotPop <- p1 + labs(y = "Morph", x = "Hatching (DD)") + 
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1 <- plotPop + coord_flip()
p1


# plot differences between morphs in hatching time
p1 <- ggplot(mydata, aes(morph, HDD, fill = morph))+
  geom_boxplot()+
  ylim(350,500) +
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  theme(axis.text.x=element_text(size = 80, colour = "black"), 
        axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80), axis.title.y=element_text(size = 80),
        legend.position = "bottomright")+
  theme_bw()
plotPop <- p1 + labs(y = "Hatching (DD)", x = element_blank()) + 
  theme(legend.position = "none", axis.title = element_text(size = 80), 
        axis.text.x = element_text(size=80, colour = "black"),
        axis.text.y = element_text(size=80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plotPop
setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")
pdf(file = "24b_HDD_box.pdf", width = 20, height = 19, family = "Helvetica") # defaults to 7 x 7 inches
plotPop
dev.off()


# plot relationship between egg size and hatching time
p1 <- ggplot(mydata, aes(x = eggSize, y = HDD)) +
  geom_point(aes(shape = morph, fill = morph), size = 19, stroke = 1, colour = "black") +
  ylim(350,500) +
  scale_shape_manual(values= c(21:23,25,22)) +
  geom_smooth(aes(group = morph, col = morph, fill = morph), method = "lm") + 
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  xlab("Egg size (mm)") +
  ylab("Hatching (DD)")+
  theme_bw() +
  theme(axis.text.x=element_text(size = 80, colour = "black"), axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80, colour = "black"), axis.title.y=element_text(size = 80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
p1
setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")
pdf(file = "24c_HDD~eggSize.pdf", width = 20, height = 19, family = "Helvetica") # defaults to 7 x 7 inches
p1
dev.off()



#####FF
### first plot non-egg size adjusted patterns:
m1 <- lm(FFDD ~ morph, data = mydata)
m2 <- lm(FFDD ~ eggSize, data = mydata) 
m3 <- lm(FFDD ~ morph*eggSize, data = mydata) 
m4 <- lm(FFDD ~ morph+eggSize, data = mydata) 
models <- list(m1, m2, m3, m4)
model.names <- c("morph","eggSize","morph*eggSize", "morph+eggSize")
aictab(cand.set = models, modnames = model.names)

## NB: INTERESTED IN MORPH ONLY ATM (see below for morph + egg size)
m1 <- lm(FFDD ~ morph, data = mydata)  
Anova(m1)
summary(m1)
leastsquare <- lsmeans(m1, ~morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, flip = T, reverse = T) 


## is time to FF influenced by morph or egg size?
m1 <- lm(FFDD ~ morph + eggSize, data = mydata) # NS interaction, therefore dropped
Anova(m1)
summary(m1)
leastsquare <- lsmeans(m1, ~morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, flip = T, reverse = T) # for some reason it didn't work wel with previous settings!
p1 <- plot(leastsquare, comparisons = TRUE, CIs = T, adjust = "tukey")+
  theme_bw() 
plotPop <- p1 + labs(y = "Morph", x = "First feeding (DD)") + 
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1 <- plotPop + coord_flip()
p1




# plot differences between morphs in FF time
p1 <- ggplot(mydata, aes(morph, FFDD, fill = morph))+
  geom_boxplot()+
  ylim(600,750)+
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  theme(axis.text.x=element_text(size = 80, colour = "black"), 
        axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80), axis.title.y=element_text(size = 80),
        legend.position = "bottomright")+
  theme_bw()
plotPop <- p1 + labs(y = "First feeding (DD)", x = element_blank()) + 
  theme(legend.position = "none", axis.title = element_text(size = 80), 
        axis.text.x = element_text(size=80, colour = "black"),
        axis.text.y = element_text(size=80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plotPop
setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")
pdf(file = "26_FFDD_box.pdf", width = 20, height = 19, family = "Helvetica") # defaults to 7 x 7 inches
plotPop
dev.off()


# plot relationship between egg size and hatching time
p1 <- ggplot(mydata, aes(x = eggSize, y = FFDD)) +
  geom_point(aes(shape = morph, fill = morph), size = 19, stroke = 1, colour = "black") +
  ylim(600,750)+
  scale_shape_manual(values= c(21:23,25,22)) +
  geom_smooth(aes(group = morph, col = morph, fill = morph), method = "lm") + 
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  xlab("Egg size (mm)") +
  ylab("First feeding (DD)")+
  theme_bw() +
  theme(axis.text.x=element_text(size = 80, colour = "black"), axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80, colour = "black"), axis.title.y=element_text(size = 80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
p1
setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")
pdf(file = "25b_FFDD~eggSize.pdf", width = 20, height = 19, family = "Helvetica") # defaults to 7 x 7 inches
p1
dev.off()

