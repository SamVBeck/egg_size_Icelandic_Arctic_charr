######### now testing for size differences
rm(list = ls())

library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(plyr)
library(car)
library(emmeans)

mydata <- read.csv("~/Documents/PhD/R_2014:15/R_postPhD/eggSize/female_data.csv")
mydata$year <- as.factor(mydata$year)
mydata$male <- as.factor(mydata$male)
mydata$pop.morph <- factor(mydata$pop.morph, levels = c("FJ","VS","VB","SV","TP","TLB","GB"))


mypalette <- brewer.pal(n = 12, name = 'Paired')
#insert the 6th color into the correct place on the list
mypalette <- append(mypalette, mypalette[6], 0)
mypalette <- append(mypalette, mypalette[8], 3)
myshapes <- c("FJ"=21,  "VS"=23, "VB"=22,"TP"=21,"TLB"=24)



### subset by: 1) dev.stage; 2) pop.morph;
PF <- subset(mydata, dev_stage=="PF",select=lake:eggCV); PF <- droplevels(PF)
E<- subset(mydata, dev_stage=="E",select=lake:eggCV); E <- droplevels(E)
H <- subset(mydata, dev_stage=="H",select=lake:eggCV); H <- droplevels(H)
FF <- subset(mydata, dev_stage=="FF",select=lake:eggCV); FF <- droplevels(FF)


FJ <- subset(mydata, pop.morph=="FJ",select=lake:eggCV)
VB <- subset(mydata, pop.morph=="VB",select=lake:eggCV)
VS <- subset(mydata, pop.morph=="VS",select=lake:eggCV)
SV <- subset(mydata, pop.morph=="SV",select=lake:eggCV)
TLB <- subset(mydata, pop.morph=="TLB",select=lake:eggCV)
TP <- subset(mydata, pop.morph=="TP",select=lake:eggCV)
GB <- subset(mydata, pop.morph=="GB",select=lake:eggCV)

EggSizes <- subset(mydata, dev_stage == "PF" | dev_stage == "E")
EggSizes <- droplevels(EggSizes)

FJ <- droplevels(FJ); VB <- droplevels(VB); VS <- droplevels(VS)
SV <- droplevels(SV); TLB <- droplevels(TLB); TP <- droplevels(TP)
GB <- droplevels(GB)


H$female <- as.factor(H$female)
H <- H[ ! H$female %in% 5, ]
H <- H[ ! H$female %in% 7, ]
H <- droplevels(H)
levels(H$female)

# analyse how size at H varied among populations 
m1 <- lm(meanHSize ~ pop.morph, data = H)
Anova(m1)
summary(m1)
leastsquare <- lsmeans(m1, ~pop.morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, flip = T, reverse = T) # for some reason it didn't work wel with previous settings!
p1 <- plot(leastsquare, comparisons = TRUE, CIs = T, adjust = "tukey")+
  theme_bw() 
plotPop <- p1 + labs(y = "Morph", x = "Size at hatching (mm)") + 
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1 <- plotPop + coord_flip()
p1




# 1) is size at hatching influenced by egg size and/or morphs?
m1 <- lm(meanHSize ~ meanPFandE + pop.morph, data = H)
Anova(m1)
summary(m1)
leastsquare <- lsmeans(m1, ~pop.morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, flip = T, reverse = T) # for some reason it didn't work wel with previous settings!
p1 <- plot(leastsquare, comparisons = TRUE, CIs = F, adjust = "tukey")+
  theme_bw() 
plotPop <- p1 + labs(y = "Morph", x = "Size at hatching (mm)") + 
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1 <- plotPop + coord_flip()
p1


# lts plot H size diffrences between morphs
p1 <- ggplot(EggSizes, aes(pop.morph, meanHSize, fill = pop.morph))+
  geom_boxplot()+
  ylim(12,18)+
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  theme(axis.text.x=element_text(size = 80, colour = "black"), 
        axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80), axis.title.y=element_text(size = 80),
        legend.position = "bottomright")+
  theme_bw()
plotPop <- p1 + labs(y = "Size at hatching (mm)", x = element_blank()) + 
  theme(legend.position = "none", axis.title = element_text(size = 80), 
        axis.text.x = element_text(size=80, colour = "black"),
        axis.text.y = element_text(size=80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plotPop
setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")
pdf(file = "26b_Hsize_box.pdf", width = 20, height = 19, family = "Helvetica") # defaults to 7 x 7 inches
plotPop
dev.off()




### lets plot how egg size correlates with H size among morphs
p1 <- ggplot(H, aes(x = meanPFandE, y = meanHSize)) +
  geom_point(aes(shape = pop.morph, fill = pop.morph), size = 19, stroke = 1, colour = "black") +
  ylim(12,18)+
  scale_shape_manual(values= c(21:25,21,22)) +
  geom_smooth(aes(group = pop.morph, col = pop.morph, fill = pop.morph), method = "lm") + 
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  xlab("Egg size (mm)") +
  ylab("Size at hatching (mm)")+
  theme_bw() +
  theme(axis.text.x=element_text(size = 80, colour = "black"), axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80, colour = "black"), axis.title.y=element_text(size = 80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
p1
setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")
pdf(file = "26b_Hsize_eggSize.pdf", width = 20, height = 19, family = "Helvetica") # defaults to 7 x 7 inches
p1
dev.off()




## look at differences in FF size between morphs without the influence of egg size
m1 <- lm(meanFFSize ~ pop.morph, data = FF)
Anova(m1)
summary(m1)
leastsquare <- lsmeans(m1, ~pop.morph, adjust = "tukey")
test(leastsquare) 
pwpm(leastsquare, flip = T, reverse = T) 
p1 <- plot(leastsquare, comparisons = TRUE, CIs = F, adjust = "tukey")+
  theme_bw() 
plotPop <- p1 + labs(y = "Morph", x = "Size at hatching (mm)") + 
  theme(legend.position = "none", axis.title = element_text(size = 40), 
        axis.text.x = element_text(size=40, colour = "black"),
        axis.text.y = element_text(size=40, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1 <- plotPop + coord_flip()
p1



# 2) is size at first feeding influenced by egg size and/or morphs?
m1 <- lm(meanFFSize ~ meanPFandE * pop.morph, data = FF)
Anova(m1)
summary(m1)
emm_s.t <- emmeans(m1, pairwise ~ pop.morph | meanPFandE) # Obtain slopes for each pop.morph
( fit1.emt <- emtrends(m1, "pop.morph", var = "meanPFandE") )
emtrends(m1, pairwise ~ pop.morph, var = "meanPFandE")
summary(emtrends(m1, ~"pop.morph", var = "meanPFandE"), infer = T)
pwpm(fit1.emt, flip = T, reverse = T) 
emmip(m1, pop.morph ~ meanPFandE, at = list(meanPFandE = range(FF$meanPFandE)))

# plot the slopes
setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")
pdf(file = "33_FFinteraction_.pdf", width = 20, height = 19, family = "Helvetica") # defaults to 7 x 7 inches
p1 <- emmip(m1, pop.morph ~ meanPFandE, at = list(meanPFandE = range(FF$meanPFandE)))+
  geom_line(aes(group = pop.morph), size = 1.5) +
  ylab("Size at first feeding (mm)") + xlab("Mean egg diameter (mm)") +
  geom_point(aes(shape = pop.morph, fill = pop.morph), size = 10, stroke = 1, colour = "black") +
  scale_shape_manual(values= c(21:25,21,22)) +
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0)) +
  #ggtitle("Mean size at FF related to mean egg size") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text.x=element_text(colour="black", size = 40),axis.title.x=element_text(size = 40, colour = "black"),
                     axis.text.y=element_text(colour="black", size = 40),axis.title.y=element_text(size = 40, colour = "black"),
                     plot.title = element_text(size = 40),
                     legend.position = "none")
p1
dev.off()

# lets plot FF size differences
p1 <- ggplot(EggSizes, aes(pop.morph, meanFFSize, fill = pop.morph))+
  geom_boxplot()+
  ylim(15,28)+
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  theme(axis.text.x=element_text(size = 80, colour = "black"), 
        axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80), axis.title.y=element_text(size = 80),
        legend.position = "bottomright")+
  theme_bw()
plotPop <- p1 + labs(y = "Size at first feeding (mm)", x = element_blank()) + 
  theme(legend.position = "none", axis.title = element_text(size = 80), 
        axis.text.x = element_text(size=80, colour = "black"),
        axis.text.y = element_text(size=80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plotPop
setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")
pdf(file = "26_FFsize_box.pdf", width = 20, height = 19, family = "Helvetica") # defaults to 7 x 7 inches
plotPop
dev.off()




### lets plot how egg size correlates with FF size among morphs
p1 <- ggplot(FF, aes(x = meanPFandE, y = meanFFSize)) +
  geom_point(aes(shape = pop.morph, fill = pop.morph), size = 19, stroke = 1, colour = "black") +
  ylim(15,28)+
  scale_shape_manual(values= c(21:25,21,22)) +
  geom_smooth(aes(group = pop.morph, col = pop.morph, fill = pop.morph), method = "lm") + 
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  xlab("Egg size (mm)") +
  ylab("Size at first feeding (mm)")+
  theme_bw() +
  theme(axis.text.x=element_text(size = 80, colour = "black"), axis.title.x=element_text(size = 80, colour = "black"),
        axis.text.y=element_text(size = 80, colour = "black"), axis.title.y=element_text(size = 80, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")
p1
setwd("~/Documents/PhD/Figures/44_POST_PHD/02_eggSize/")
pdf(file = "26b_FFsize_eggSize.pdf", width = 20, height = 19, family = "Helvetica") # defaults to 7 x 7 inches
p1
dev.off()
