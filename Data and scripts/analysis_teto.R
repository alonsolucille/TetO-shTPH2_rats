#R version 4.2.3 (2023-03-15 ucrt) -- "Shortstop Beagle"
#RStudio 2023.03.0+386 "Cherry Blossom"

rm(list = ls())

# environment ####---------------------------------------------------------------
# install.packages("renv")
# Step 1: initialize a project level R library
# renv::init()
# Step 2: save the current status of your library to a lock file
# renv::snapshot()#to run after making changes to your project, such as installing new packages
# Step 3: restore state of your project from renv.lock
renv::restore()
# libraries ####---------------------------------------------------------------

library(tidyverse)#version 2.0.0
# v ggplot2 3.4.2     v purrr   1.0.1
# v tibble  3.2.1     v dplyr   1.1.1
# v tidyr   1.3.0     v stringr 1.5.0
# v readr   2.1.4     v forcats 1.0.0
library(gridExtra)#version2.3
library(lmPerm)#version2.1.0
library(ggbeeswarm)
library(grid)

# themes ####------------------------------------------------------------------
theme.teto <-  theme (axis.title = element_text(size=13,face="plain", hjust = .5,vjust = 0),
                      plot.title = element_text(size=12, hjust = 0),
                      #plot.margin = c(2,2,2,2),#margin(2,2,2,2, unit="pt")
                      #axis.text = element_text(margin = margin(6,6,6,6, unit="pt")),
                      axis.text.y = element_text(size = 11,angle=0,colour="grey16"),#, face = 'bold', hjust=-0.25),
                      axis.text.x = element_text(size = 11,angle=0,colour="grey16"),#,face = 'bold'
                      axis.line=element_line(colour="grey16"),
                      axis.line.y=element_line(linewidth=.3,colour="grey16"),
                      axis.line.x=element_line(linewidth=.3, colour="grey16"),
                      # strip.background = element_rect(colour="grey16", fill="white",linetype = 0),
                      strip.text = element_blank(),
                      #strip.text.x = element_text(size=14),
                      #strip.text.y = element_text(size=14),
                      #legend.title = element_text(size=12),
                      legend.title = element_blank(),
                      legend.text =element_text(size=12),#, family = "sans"
                      legend.position="right",
                      axis.ticks=element_line(linewidth=.3,colour="grey16"),
                      axis.ticks.length = unit(.1, "cm"))#.15

col.teto <- c("#331068FF","#DE4968FF")
col.teto2 <- c("#331068FF","#331068FF","#331068FF","#DE4968FF","#DE4968FF","#DE4968FF")

col.dm <- c("Control-gdm" = "white","Control-int" ="white","Control-pdm" ="white","Tph2-kd-gdm" ="white","Tph2-kd-int" ="white","Tph2-kd-pdm" ="#000004FF")

shap.dm <- c("Control-gdm" = 24,"Control-int" =22,"Control-pdm" =25,"Tph2-kd-gdm" =24,"Tph2-kd-int" =22,"Tph2-kd-pdm" =25)

carre = grobTree(rectGrob(gp = gpar(fill = "white", alpha=0)))#when white space is needed in a figure

# data ####--------------------------------------------------------------------
data1 <- read.csv2("teto_all_16082021.csv")

#rename groups and factor
data1$group <- factor(data1$group,levels = c("CONTROL", "TETO-DOX"), labels = c("Control","Tph2-kd"))

data1$groupdm <- interaction(data1$group,data1$dm)
data1$groupdm <- factor(data1$groupdm,
                      levels = c("Control.GDM","Control.INT","Control.PDM", "Tph2-kd.GDM","Tph2-kd.INT","Tph2-kd.PDM"), 
                      labels = c("Control-gdm","Control-int","Control-pdm", "Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm"))
data1$rat <- factor(data1$rat)
data1$batch <- factor(data1$batch)
data1$genotype <- factor(data1$genotype)
data1$treatment <- factor(data1$treatment)
data1$dm <- factor(data1$dm,levels = c("GDM","INT","PDM","OUT"))

turnover.hyp <- read.csv2("serotonin_19052023.csv", sep=",", dec = ".")

turnover.hyp$rat <- factor(turnover.hyp$rat)
turnover.hyp$group <- factor(turnover.hyp$group, 
                             levels = c("Tph2-wt","Tph2-kd"),
                             labels = c("Control","Tph2-kd"))
turnover.hyp$groupx <- factor(turnover.hyp$groupx, levels = c("SD-dox","TETO-water","TETO-dox"))
turnover.hyp$batch <- factor(turnover.hyp$batch)
turnover.hyp$genotype <- factor(turnover.hyp$genotype)
turnover.hyp$treatment <- factor(turnover.hyp$treatment)
turnover.hyp$expe.group <- factor(turnover.hyp$expe.group)

#result serotonin####
#calculate ratios (turnover of serotonin = ratio, normalization serotonin by TRP = ratio2, normalization 5HIAA by TRP = ratio3)
turnover.hyp <- turnover.hyp %>%
  mutate(ratio = X5HIAA/X5HT, 
         ratio2 = X5HT/TRP, 
         ratio3 = X5HIAA/TRP)

#calculate means per batch for control groups
turnover.hyp.cont.mean <- turnover.hyp %>%
  filter(group == "Control") %>%
  group_by(genotype,treatment,batch) %>%
  mutate(M5HT = mean(X5HT, na.rm = TRUE),
         M5HIAA = mean(X5HIAA, na.rm = TRUE),
         MTRP = mean(TRP, na.rm = TRUE)) %>%
  ungroup()

#mean value for control groups per batch
tt <- turnover.hyp.cont.mean %>%
  filter(!is.na(M5HT)) %>%
  group_by(batch) %>%
  summarise(M5HT2 = unique(M5HT), 
            M5HIAA2 = unique(M5HIAA),
            MTRP2 = unique(MTRP)) %>%
  ungroup()

#join and fill means for Tph2-kd per batch
tt2 <- full_join(turnover.hyp,tt)

#calculate the percentages: normalized values for all individuals
turnover.hyp.norm <- tt2 %>% 
  mutate(Norm.5HT = X5HT/M5HT2*100,
         Norm.5HIAA = X5HIAA/M5HIAA2*100,
         Norm.TRP = TRP/MTRP2*100) %>%
  ungroup()

#control for outliers
turnover.hyp.norm <- turnover.hyp.norm%>%
  group_by(batch)%>%
  mutate(sd=(sd(ratio, na.rm=T)),
         mean=mean(ratio,na.rm=T))%>%
  ungroup()

#are they outliers?
turnover.hyp.norm <- turnover.hyp.norm%>%
  mutate(is.out=ratio>(mean+2*sd))
turnover.hyp.norm$is.out#rat#230 is an outlier for the serotonin turnover (R#230 is in batch 12)

#calculate the mean serotonin decrease and sd for kd individuals
turnover.hyp.norm_kd<-turnover.hyp.norm %>% filter(group=="Tph2-kd", batch!=12)#remove outlier and B12
100-mean(turnover.hyp.norm_kd$Norm.5HT, na.rm=TRUE)
sd(turnover.hyp.norm_kd$Norm.5HT, na.rm=TRUE)

#calculate the mean 5HIAA decrease and sd for kd individuals
turnover.hyp.norm_kd2<-turnover.hyp.norm %>% filter(group=="Tph2-kd", rat!=230)#remove outlier
100-mean(turnover.hyp.norm_kd2$Norm.5HIAA, na.rm=TRUE)
sd(turnover.hyp.norm_kd2$Norm.5HIAA, na.rm=TRUE)

#fuse data (behavior and serotonin)####
data2 <- full_join(data1,turnover.hyp.norm, by=c("rat", "batch","genotype","treatment","group", "groupx"))

#data filtered (data3 used fro later analysis) ####
data3 <- data2%>%filter(rat!=230) #remove the outlier from analysis

nb.data3 <- data3 %>%
  group_by(group) %>%
  summarise(Count = n()) %>%
  spread(group,Count)
nb.data3

#data without B12
datab12out<-data3%>% filter(batch!=12)

nb.datab12out <- datab12out %>%
  group_by(group) %>%
  summarise(Count = n()) %>%
  spread(group,Count)
nb.datab12out

#H1: treatment reduces 5-HT ####
p.5ht<-ggplot(data = datab12out, aes(x=group, y=X5HT)) +  
  geom_boxplot(aes(color=group), fatten=1.5, lwd=.3, outlier.size = .4)+#remove outliers in background,outlier.shape = NA
  geom_beeswarm(aes(color=group), size=2, stroke=.5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic()+
  theme.teto+theme(legend.position ="none")+
  geom_segment(aes(x=1, xend=2, y=2150, yend=2150), colour="black", linewidth=.5) +
  annotate("text", x=1.5, y=2200, label="*", colour="black", size=5) +
  scale_color_manual(values = col.teto)+
  labs(x="",y="5-HT (pg/mg)", title = "5-HT - B12 out")
p.5ht

wilcox.test(datab12out$X5HT~datab12out$group)#W = 2208.5, p-value = 1.289e-06 serotonin is decreased in Tph2-kd

p.5HIAA<-ggplot(data = data3, aes(x=group, y=X5HIAA)) +  
  geom_boxplot(aes(color=group), fatten=1.5, lwd=.3, outlier.size = 1)+#remove outliers in background,outlier.shape = NA
  geom_beeswarm(aes(color=group), size=2, stroke=.5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic()+
  theme.teto+theme(legend.position = "none")+
  geom_segment(aes(x=1, xend=2, y=2850, yend=2850), colour="black", linewidth=.5) +
  annotate("text", x=1.5, y=2900, label="*", colour="black", size=5) +
  scale_color_manual(values = col.teto) +
  labs(x="",y="5-HIAA (pg/mg)", title = "B")
p.5HIAA

wilcox.test(data3$X5HIAA~data3$group)#W = 2649.5, p-value = 9.92e-07 YES 5HIAA is decreased in Tph2-kd

p.TRP<-ggplot(data = data3, aes(x=group, y=TRP)) +  
  geom_boxplot(aes(color=group), fatten=1.5, lwd=.3, outlier.size = 1)+#remove outliers in background,outlier.shape = NA
  # geom_beeswarm(aes(color=group), size=2, stroke=.5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic()+
  theme.teto+  theme(legend.position = "none") +
  scale_color_manual(values = col.teto)+
  labs(x="",y="TRP (pg/mg)", title = "C")
p.TRP

wilcox.test(data3$TRP~data3$group)#W = 1840, p-value = 0.5922 TRP levels are unchanged in Tph2-kd

#normalized figures####
#figure 5A 5-HT norm
p.n.5ht<-ggplot(data = datab12out, aes(x=group, y=Norm.5HT)) +  
  geom_boxplot(aes(color=group), fatten=1.5, lwd=.5, outlier.size = .5)+#remove outliers in background,outlier.shape = NA
  # geom_beeswarm(aes(color=group), size=2, stroke=.5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic() +
  theme.teto+theme(legend.position ="none") +
  geom_segment(aes(x=1, xend=2, y=178, yend=178), colour="black", linewidth=.5) +
  annotate("text", x=1.5, y=183, label="*", colour="black", size=5) +
  scale_y_continuous(limits = c(0,183),breaks=seq(0,200,50))+
  scale_x_discrete(labels= c("Tph2-wt","Tph2-kd")) +
  scale_color_manual(values = col.teto) +
  labs(x="",y="5-HT % of control levels", title = "A")#5-HT (%) - B12 out
p.n.5ht

#excluded animal (R#230) belongs to B12, so removing B12 for 5HT data also removes R#230:
wilcox.test(datab12out$Norm.5HT~datab12out$group)#W = 2259, p-value = 2.523e-07

#figure 5B 5-HIAA
p.n.5HIAA<-ggplot(data = data3, aes(x=group, y=Norm.5HIAA)) +  
  geom_boxplot(aes(color=group), fatten=1.5, lwd=.5, outlier.size = .5)+#remove outliers in background,outlier.shape = NA
  # geom_beeswarm(aes(color=group), size=2, stroke=.5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic()+
  theme.teto+theme(legend.position = "none") +
  geom_segment(aes(x=1, xend=2, y=178, yend=178), colour="black", linewidth=.5) +
  annotate("text", x=1.5, y=183, label="*", colour="black", size=5) +
  scale_y_continuous(limits = c(0,183),breaks=seq(0,200,50))+
  scale_x_discrete(labels= c("Tph2-wt","Tph2-kd")) +
  scale_color_manual(values = col.teto) +
  labs(x="",y="5-HIAA % of control levels", title = "B")
p.n.5HIAA

wilcox.test(data3$Norm.5HIAA~data3$group)#W = 2825, p-value = 5.286e-09

#figure 5C TRP
p.n.TRP<-ggplot(data = data3, aes(x=group, y=Norm.TRP)) +  
  geom_boxplot(aes(color=group), fatten=1.5, lwd=.5, outlier.size = .5)+#remove outliers in background,outlier.shape = NA
  # geom_beeswarm(aes(color=group), size=2, stroke=.5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic()+
  theme.teto+  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0,183),breaks=seq(0,200,50)) +
  scale_x_discrete(labels= c("Tph2-wt","Tph2-kd")) +
  scale_color_manual(values = col.teto) +
  labs(x="",y="TRP % of control levels", title = "C")
p.n.TRP

wilcox.test(data3$Norm.TRP~data3$group)#W = 2034, p-value = 0.1141

F5.5HT.n<-grid.arrange(p.n.5ht, p.n.5HIAA, p.n.TRP, ncol=3)

#save figure 5 ####
# ggsave("F4-20240430.png", plot = F5.5HT.n, device="png", dpi=300,height = 10, width = 20,units = "cm")

#Fig S3 supp figure all data included ratio3####
p.ratio3<-ggplot(data = data3, aes(x=batch, y=(X5HIAA/TRP))) +  
  geom_boxplot(aes(color=group), fatten=1.5, lwd=.5, outlier.shape = NA)+#remove outliers in background,outlier.shape = NA
  geom_beeswarm(aes( color=group, fill=groupdm, shape=groupdm), size=2, stroke=.5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic()+
  theme.teto+
  theme(legend.position = "right",
        # axis.ticks.x = element_blank(),axis.line.x = element_blank(),
        # axis.text.x = element_blank(),axis.title.x =element_blank(),
        strip.text= element_text(color="black",size=10),
        strip.background = element_blank(),
        plot.title = element_text(size = 10))+ #
  geom_rect(aes(xmin = 0.5, ymin = 0.15, xmax =2.5, ymax = 0.71), fill="transparent",color="black", linetype = 3) +
  geom_rect(aes(xmin = 4.5, ymin = 0.15, xmax =6.5, ymax = 0.71), fill="transparent",color="black", linetype = 3) +
  geom_rect(aes(xmin = 6.5, ymin = 0.15, xmax =8.5, ymax = 0.71), fill="transparent",color="black", linetype = 3) +
  geom_rect(aes(xmin = 8.5, ymin = 0.15, xmax =10.5, ymax = 0.71), fill="transparent",color="black", linetype = 3) +
  annotate("text", x=1.5, y=0.73, label="25", colour="black", size=3.5) + 
  annotate("text", x=3, y=0.73, label="26", colour="black", size=3.5) + 
  annotate("text", x=4, y=0.73, label="29", colour="black", size=3.5) + 
  annotate("text", x=5.5, y=0.73, label="46", colour="black", size=3.5) +
  annotate("text", x=7.5, y=0.73, label="47", colour="black", size=3.5) + 
  annotate("text", x=9.5, y=0.73, label="48", colour="black", size=3.5) +
  scale_x_discrete(limits=factor(c(16,17,18,13,11,12,14,15,19,20))) +
  scale_fill_manual(values = col.dm, labels = c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm")) +
  scale_shape_manual(values = shap.dm, labels = c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm")) +
  scale_color_manual(values = col.teto, guide=F,, labels = c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm"))+
  guides(fill = guide_legend(override.aes=list(alpha=1,size=3, col=c("#331068FF","#331068FF","#331068FF","#DE4968FF","#DE4968FF","#DE4968FF"))))+
  labs(x="Batch",y="Ratio 5-HIAA/TRP", title = "")
p.ratio3

#save figure S3 ####
# ggsave("S3-20240430.png", plot = p.ratio3, device="png", dpi=300,height = 14, width = 30,units = "cm")

wilcox.test(data3$ratio3~data3$group)#W = 2564, p-value = 9.296e-06

stat5HT3<-aovp(ratio3~group+treat.dur.d+Error(batch), data = data3, perm="Exact", Ca=0.01, maxIter=1000000)
summary(stat5HT3)
# Error: batch
# Component 1 :
#   Df R Sum Sq R Mean Sq Pr(Exact)
# group        1  0.00499  0.004990    0.7642
# treat.dur.d  1  0.02118  0.021178    0.7190
# Residuals    7  0.79594  0.113706          
# 
# 
# Error: Within
# Component 1 :
#   Df R Sum Sq R Mean Sq  Iter  Pr(Prob)    
# group       1  0.28670  0.286704 1e+06 < 2.2e-16 ***
#   Residuals 107  0.40527  0.003788                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# treatment, F(1,107) = 0.286704/0.003788 = 75, p-value < 0.001, 
# duration, F(1,7) = 0.004990/0.021178 = 0.23, p-value = 0.719

#plots for each batch ####
datab11 <- data2%>% filter(batch == 11)
datab13 <- data2%>% filter(batch == 13)
datab14 <- data2%>% filter(batch == 14)
datab15 <- data2%>% filter(batch == 15)
datab16 <- data2%>% filter(batch == 16)
datab17 <- data2%>% filter(batch == 17)
datab18 <- data2%>% filter(batch == 18)
datab19 <- data2%>% filter(batch == 19)
datab20 <- data2%>% filter(batch == 20)

# change data name to plot another batch
p.5HT.b18<-ggplot(data = datab18, aes(x=group, y=(ratio))) +  
  geom_boxplot(aes(color=group), fatten=1.5, lwd=.3, outlier.size = .4)+#remove outliers in background,outlier.shape = NA
  geom_beeswarm(aes(fill = groupdm, shape=groupdm, col = group),size=1,cex = 3.5, alpha=.8,  priority = "random",  show.legend = F) +
  theme_classic() +
  theme.teto +
  theme(legend.position = "right", strip.text = element_text(size=6), strip.background = element_blank()) +
  scale_fill_manual(values = col.dm) +
  scale_shape_manual(values = shap.dm) +
  scale_color_manual(values = col.teto) +
  labs(y="ratio 5-HIAA/5-HT", title = "Individual ratio 5-HIAA/5-HT")
p.5HT.b18

#H2 - equivalence of controls ####

#each batch 
p.5ht.cont<-ggplot(data = filter(datab12out, group == "Control"), aes(x = groupx, y = X5HT)) +  
  geom_boxplot(aes(color = groupx), fatten=1.5, lwd=.3, outlier.size = .4) +#remove outliers in background,outlier.shape = NA
  geom_beeswarm(aes(color = groupx), size=1, stroke=.5, cex = 3, priority = "density", groupOnX = T) +
  theme_classic() +
  theme.teto+theme(legend.position = "right", strip.text = element_text(size=6), strip.background = element_blank()) +
  facet_grid(~batch, scale="free_x") +
  labs(x = "",y = "5-HT (pg/mg)", title = "5-HT")
p.5ht.cont

wilcox.test(X5HT~groupx, data = filter(datab12out, group == "Control"))#W = 546.5, p-value = 0.001204
wilcox.test(TRP~groupx, data = filter(data3, group == "Control"))#W = 507.5, p-value = 0.2578

stat5HT.cont<-aovp(X5HT~groupx+Error(rat), data = filter(data3, group == "Control"),
               perm = "Exact", Ca = 0.01, maxIter = 1000000)
# Error: rat
# Component 1 :
#   Df R Sum Sq R Mean Sq  Iter Pr(Prob)    
# groupx1    1  2989623   2989623 1e+06  9.1e-05 ***
#   Residuals 58 10711307    184678                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#each pairs of batch
p.5ht.cont2<-ggplot(data = filter(datab12out, group == "Control"), aes(x = groupx, y = X5HT)) +  
  geom_boxplot(aes(color = groupx), fatten=1.5, lwd=.3, outlier.size = .4)+#remove outliers in background,outlier.shape = NA
  geom_beeswarm(aes(color = groupx), size = 1, stroke=.5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic() +
  theme.teto+theme(legend.position = "right", strip.text = element_text(size = 6), strip.background = element_blank()) +
  facet_grid(~expe.group, scale = "free_x") +
  scale_shape_manual(values = c(24,22,25,24,22,25)) +
  labs(x="",y="5-HT (pg/mg)", title = "5-HT")
p.5ht.cont2

p.trp.cont<-ggplot(data = filter(data3, group=="Control"), aes(x=groupx, y=TRP)) +  
  geom_boxplot(aes(color=groupx), fatten=1.5, lwd=.3, outlier.size = .4)+#remove outliers in background,outlier.shape = NA
  geom_beeswarm(aes(color=groupx), size=1, stroke=.5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic() +
  theme.teto+theme(legend.position ="right", strip.text = element_text(size=6), strip.background = element_blank()) +
  facet_grid(~batch, scale="free_x") +
  scale_shape_manual(values = c(24,22,25,24,22,25)) +
  labs(x="",y="TRP (pg/mg)", title = "TRP")
p.trp.cont

statTRP.cont<-aovp(TRP~groupx + Error(rat), 
                   data = filter(data3, group=="Control"),
                   perm="Exact", Ca=0.01, maxIter=1000000)
summary(statTRP.cont)
# Error: rat
# Component 1 :
#   Df R Sum Sq R Mean Sq  Iter Pr(Prob)  
# groupx1    1  1904159   1904159 97493  0.09303 .
# Residuals 58 36779995    634138                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#equivalence test
library("TOSTER")

data4<-filter(datab12out, group=="Control")
t_TOST(ratio~groupx, data = data4, eqb=.2)
# Control groups are not equivalent, SD-dox control batches have variable levels of 5-HT
# Welch Two Sample t-test
# 
# The equivalence test was non-significant, t(36.58) = -0.15, p = 0.44
# The null hypothesis test was significant, t(36.58) = 3.091p < 0.01
# NHST: reject null significance hypothesis that the effect is equal to zero 
# TOST: don't reject null equivalence hypothesis
# 
# TOST Results 
#                  t    df p.value
# t-test      3.0913 36.58   0.004
# TOST Lower  6.3350 36.58 < 0.001
# TOST Upper -0.1525 36.58    0.44
# 
# Effect Sizes 
#                Estimate      SE             C.I. Conf. Level
# Raw              0.1906 0.06166 [0.0865, 0.2946]         0.9
# Hedges's g(av)   0.8506 0.38598 [0.3665, 1.3247]         0.9
# Note: SMD confidence intervals are an approximation. See vignette("SMD_calcs").

#serotonin levels of Tph2-kd each batch
p.5ht.teto<-ggplot(data = filter(datab12out, group!="Control"), aes(x=groupx, y=X5HT)) +  
  geom_boxplot(aes(color=groupx), fatten = 1.5, lwd = .3, outlier.size = .4)+#remove outliers in background,outlier.shape = NA
  geom_beeswarm(aes(color=groupx), size = 1, stroke = .5, cex = 3, priority = "density", groupOnX = T)+
  theme_classic() +
  theme.teto+theme(legend.position = "right", strip.text = element_text(size = 6), strip.background = element_blank()) +
  facet_grid(~batch, scale = "free_x") +
  scale_shape_manual(values = c(24,22,25,24,22,25)) +
  scale_y_continuous(limits = c(0,2000))+
  # scale_color_manual(values = col.teto)+
  labs(x="",y="5-HT (pg/mg)", title = "5-HT")
p.5ht.teto

#Behavioral results  ####
#Rat Gambling Task ####

#last 20 min ####
data.rgt<-data3%>%filter(dm%in%c("GDM","PDM","INT"))#also use this dataset for dm analysis of other tests

nb<-data.rgt %>%
  group_by(group,dm) %>%
  summarise(Count = n()) %>%
  spread(group,Count) %>%
  filter(dm %in% c("GDM","PDM","INT"))
nb
# dm    Control `Tph2-kd`
# <fct>   <int>     <int>
# 1 GDM        49        39
# 2 INT         8        10
# 3 PDM         3         9

#distribution of RGT scores
shapiro.test(data.rgt$RGT)#non normal data W = 0.75757, p-value = 1.122e-12
qqnorm(data.rgt$RGT);qqline(data.rgt$RGT)

p.RGT1.1<-ggplot(data = data.rgt, aes(x= group, y=RGT)) +  
  geom_hline(yintercept=70, colour="#333333", linetype="dashed", linewidth=.3) +
  geom_hline(yintercept=30, colour="#333333", linetype="dashed", linewidth=.3) +
  geom_boxplot(aes(colour=group), outlier.shape = NA, alpha=0,show.legend = F, 
               fatten=1.5,lwd=.5) +#,alpha=.6#remove outliers in background,color="black",
  geom_segment(aes(x=1, xend=2, y=105, yend=105), colour="black", linewidth=.5) +
  annotate("text", x=1.5, y=106, label="*", colour="black", size=5) +
  geom_beeswarm(aes(fill = groupdm, shape=groupdm, col = group),size=1,cex = 3.5,
                alpha=.6,  priority = "random", stroke=1, show.legend = F) +#stroke=.5,groupOnX = T,color="black",
  theme_classic() +
  theme.teto +
  theme(legend.position = "bottom", #"none"c(.9,.5)
        axis.ticks.x = element_line(colour = "transparent"), 
        #axis.line.x = element_line(colour = "transparent"),
        axis.text.x = element_text(colour = "transparent"), 
        axis.title.x = element_text(colour = "transparent"),
        plot.margin = margin(t=4,r=0, unit = "pt"),
        plot.background = element_rect(color="transparent")) +
  scale_fill_manual(values = col.dm) +
  scale_color_manual(values = col.teto) +
  scale_shape_manual(values = c(24,22,25,24,22,25)) +
  scale_y_continuous(limits = c(0,108),breaks=seq(0,100,10), 
                     labels = c(0,"",20,"",40,"",60,"",80,"",100)) +
  scale_x_discrete() +
  guides(shape = "none") +
  labs(x = "Group",y = "RGT score (%)", title = "B")
p.RGT1.1

p.RGT1.2 <- ggplot(data = data.rgt, aes(y = RGT, fill = group, colour = group)) +  
  geom_hline(yintercept = 70, colour = "#333333", linetype = "dashed", linewidth=.3)+
  geom_hline(yintercept = 30, colour = "#333333", linetype = "dashed", linewidth=.3)+
  geom_density(data = data.rgt, aes(y=RGT, fill = group, colour = group), 
               position = "identity", alpha=0.5,size=.5,trim=T) +#trim stop the distribution at 100 precisely
  theme_classic() +
  theme.teto +
  theme(legend.position = "none", #"none""bottom"c(.9,.5)c(.7,.5)
        axis.ticks.x = element_line(colour="transparent"), 
        #axis.line.x = element_line(colour="transparent"), 
        axis.text.x= element_text(colour="transparent"),
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), 
        axis.text.y= element_blank(),
        axis.title.y= element_blank(),
        plot.title = element_text(color="transparent"),
        axis.title.x= element_text(colour="transparent"),
        plot.margin = margin(t=4,l=0, unit="pt"),
        plot.background = element_rect(color="transparent")) +
  scale_fill_manual(values = col.teto, labels=c("Tph2-wt","Tph2-kd")) +
  scale_color_manual(values = col.teto, labels=c("Tph2-wt","Tph2-kd")) +
  guides(fill = "none") +
  scale_y_continuous(limits= c(0,108),breaks=seq(0,100,10),labels = c(0,"",20,"",40,"",60,"",80,"",100)) +
  labs(x="Group",y="RGT score (%)", title = "B")
p.RGT1.2

stat.RGT1<-wilcox.test(data.rgt$RGT~group,data=data.rgt)# express stat in text: stat.RGT2$statistic, stat.RGT2$p.value
# Wilcoxon rank sum test with continuity correction
# 
# data:  data.rgt$RGT by group
# W = 2112.5, p-value = 0.04518
# alternative hypothesis: true location shift is not equal to 0

#proportion of decision maker categories ####
count_GDM_INT_PDM<- matrix(c(49,39,8,10,3,9), ncol = 2,
                           dimnames=list(c("GDM","INT","PDM"),c("control","test")),byrow=TRUE)

stat.RGT2<-fisher.test(count_GDM_INT_PDM,simulate.p.value = TRUE)# express stat in text: stat.RGT1$p.value cited in the text
# Fisher's Exact Test for Count Data with simulated p-value (based on 2000 replicates)
# 
# data:  count_GDM_INT_PDM
# p-value = 0.1079
# alternative hypothesis: two.sided

#proportion of poor decision makers compared to good decision makers in control and tph2-kd ####
count_GDM_PDM<- matrix(c(49,39,3,9), ncol = 2,
                       dimnames=list(c("GDM","PDM"),c("control","test")),byrow=TRUE)

stat.RGT2.2<-fisher.test(count_GDM_PDM,alternative = "greater")
# Fisher's Exact Test for Count Data
# 
# data:  count_GDM_PDM
# p-value = 0.04473
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.032936      Inf
# sample estimates:
# odds ratio 
#   3.721382 

# kinetics of the decision maker groups ####
RGT <- data.rgt %>% gather(time,perf,RGT10:RGT60)
RGT<- RGT%>%mutate(time=substr(time, 4, 5))
pd <- position_dodge(1)
RGTgroup<- RGT%>% group_by(time,group)%>%summarise(Mean=mean(perf, na.rm=TRUE), SD=sd(perf, na.rm=TRUE), 
                                                   Median=median(perf, na.rm = TRUE), 
                                                   Quant5=quantile(perf,c(.05), na.rm=TRUE), 
                                                   Quant95=quantile(perf,c(.95), na.rm=TRUE))

RGTgroupdm<- RGT%>% group_by(time,group,dm, groupdm)%>%summarise(Mean=mean(perf, na.rm=TRUE), SD=sd(perf, na.rm=TRUE), 
                                                                 Median=median(perf, na.rm = TRUE), 
                                                                 Quant5=quantile(perf,c(.05), na.rm=TRUE), 
                                                                 Quant95=quantile(perf,c(.95), na.rm=TRUE))
RGTgroupdm$time <- as.integer(RGTgroupdm$time)
pd2 <- position_dodge(2)
p.RGT2 <- ggplot(RGTgroupdm, aes(x = time, y = Mean, color = groupdm, group = groupdm)) + #,fill=group
  geom_hline(yintercept = 50, colour = "#333333", linetype="dashed", linewidth = .3) +
  geom_errorbar(aes(ymin = Mean, ymax = Mean+SD), linewidth = .3, width = 0, position = pd2) +
  geom_line(linewidth = .5) +#
  geom_point(aes(fill = groupdm,shape = groupdm,col = groupdm), size = 1.5, stroke = 1)+#,  stroke=.5,color="black"
  theme_classic() +
  theme.teto +
  theme(legend.position = "none")+#, plot.margin = margin(1,0,1,1, unit="pt")
  scale_y_continuous(limits = c(0,108),breaks=seq(0,100,10), labels = c("0","", "20","", "40","", "60","", "80", "", "100")) +
  scale_x_continuous(breaks=c(10,20,30,40,50,60)) +
  facet_grid(~group) +
  scale_color_manual(values = col.teto2) +#, labels=c("Tph2-wt-GDM","Tph2-wt-INT","Tph2-wt-PDM", "Tph2-kd-GDM","Tph2-kd-INT","Tph2-kd-PDM")
  scale_fill_manual(values = col.dm) +
  scale_shape_manual(values = c(24,22,25,24,22,25)) +
  guides(col = guide_legend(ncol = 2)) +
  labs(x = "Time (min)",y = "Advantageous choices (%)", title = "A")
p.RGT2

#time points stats
begin <- RGT %>% filter(time == 10)
kruskal.test(begin$perf~begin$groupdm)                 
# Kruskal-Wallis rank sum test
# data:  begin$RGT by begin$groupdm -ALL
# Kruskal-Wallis chi-squared = 10.294, df = 5, p-value = 0.06731
begin2 <- RGT %>% filter(time == 10 & groupdm %in% c("Control-gdm","Control-int","Control-pdm"))
kruskal.test(begin2$perf~begin2$groupdm)
# data:  begin$RGT by begin$groupdm -CONT
# Kruskal-Wallis chi-squared = 1.9723, df = 2, p-value = 0.373
begin3 <- RGT %>% filter(time == 10 & groupdm %in% c("Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm"))
kruskal.test(begin3$perf~begin3$groupdm)
# data:  begin$RGT by begin$groupdm -KD
# Kruskal-Wallis chi-squared = 8.2269, df = 2, p-value = 0.01635

#kinetics stats
RGT1<-RGT%>%dplyr::select(group,batch,rat,perf,time,dm,groupdm)%>%ungroup()
RGT1$time<-factor(RGT1$time)

#interaction best model :)

anova.p.rgt2<-aovp(perf~groupdm*time+Error(rat), data=RGT1,
                   perm = "Exact", Ca=0.01, maxIter=1000000)# perm="" gives the F-Statistics, perm="Exact" requires calculation MSqA/MSqE
sum.aovp1<-summary(anova.p.rgt2)

#how to calculate F-value from the table: 
#MSqA/MSqE
#https://arc.lib.montana.edu/book/statistics-with-r-textbook/item/56#The%20permutation%20p-value%20for%20the%20alternative%20hypothesis%20of%20some%20(not%20of%20greater%20or%20less%20than!)%20difference%20in%20the%20true%20means%20of%20the%20groups%20will%20involve%20counting%20the%20number%20of%20permuted%20SSA*%20results%20that%20are%20larger%20than%20what%20we%20observed.

#Latency fig S1####

p.lat.1 <- ggplot(data = data.rgt, aes(x = RGT, y = latency)) + 
  #geom_smooth(aes(x=RGT, y=flexibility_index,group=group), colour="grey", size=1,method = "glm", se = FALSE)+
  geom_boxplot(aes(colour = groupdm), outlier.shape = NA,alpha = 0,show.legend = F,fatten = 1.5, lwd = .5) +
  geom_beeswarm(aes(fill = groupdm, shape = groupdm, col = groupdm), size = 2,stroke = 1,  alpha = .6) +#stroke=0.5,groupOnX = T gives jitter to x axis#color='transparent',color="black",
  theme_classic() +
  theme.teto +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(0,100,10), labels = c("0","", "20","", "40","", "60","", "80", "", "100")) +#, limits = c(-0,103) #no limits to avoid loss of dots
  scale_shape_manual(values = c(24,22,25,24,22,25),labels = c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm")) +
  scale_fill_manual(values = col.dm,labels = c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm")) +
  scale_color_manual(values = col.teto2,labels = c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm")) +
  facet_grid(~group) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size=3),nrow = 6), shape = guide_legend(nrow = 6)) +
  labs(x = "RGT score (%)",y = "Latency for reward (s)", title = "")
p.lat.1

wilcox.test(latency~group, data = data.rgt)#W = 1725, p-value = 0.9378

anova.p.lat <- aovp(latency~group+RGT+Error(rat), data = data.rgt,
                  perm = "Exact", Ca = 0.01, maxIter = 1000000)
sum.aovp2 <- summary(anova.p.lat)
# Error: rat
# Component 1 :
#   Df R Sum Sq R Mean Sq    Iter Pr(Prob)    
# group       1   0.0155    0.0155   86385   0.1038    
# RGT         1   3.9833    3.9833 1000000   <2e-16 ***
# Residuals 115  11.9144    0.1036

# treatment, F(1,115) = 0.0155/0.1036 = 0.14, p-value = 0.1038
# RGT score, F(1,115) = 3.9833/0.1036 = 38, p-value < 0.001

#save fig S1####
# ggsave("S1-20240430.png", plot = p.lat.1, device="png", dpi=300,height = 20, width = 15,units = "cm")

#Reversed-Rat Gambling Task ####
data.rev<-data.rgt%>%filter(Do.Reverse%in%c("no","yes","undecided"))
nb.rev<-data.rev%>%
  group_by(group,dm,Do.Reverse) %>%
  summarise(Count = n()) %>%
  spread(group,Count) %>%
  filter(dm %in% c("GDM","PDM","INT"))
nb.rev

REV <- data.rev %>% gather(time,perf,REV10:REV60)
REV<- REV%>%mutate(time=substr(time, 4, 5))
REVgroup<- REV%>% group_by(time, group,dm,groupdm,Do.Reverse)%>%summarise(Mean=mean(perf, na.rm=TRUE), SD=sd(perf, na.rm=TRUE))%>%
  filter(dm!="INT")

data.rev.plot<-data.rev

pREV1<-ggplot(data = data.rev.plot, aes(x = RGT, y = flexibility_index)) + 
  geom_boxplot(aes(colour = groupdm), outlier.shape = NA,alpha = 0,show.legend = F,fatten = 1.5, lwd = .5) +
  geom_beeswarm(aes(fill = groupdm, shape = groupdm, col = groupdm), size = 1,stroke = 1,  alpha = .6) +#stroke=0.5,groupOnX = T gives jitter to x axis#color='transparent',color="black",
  theme_classic() +
  theme.teto +
  theme(legend.position = "none") +#c(.8,.2)
  scale_y_continuous(limits = c(-0,101),breaks = seq(0,100,10), labels = c("0","", "20","", "40","", "60","", "80", "", "100")) +#, limits = c(-0,103) #no limits to avoid loss of dots
  scale_x_continuous(breaks = seq(0,100,10), labels = c("0","", "20","", "40","", "60","", "80", "", "100")) +#, limits = c(-0,103) #no limits to avoid loss of dots
  scale_shape_manual(values = c(24,22,25,24,22,25)) +
  scale_fill_manual(values = col.dm) +
  scale_color_manual(values = col.teto2) +
  facet_grid(~group) +
  guides(fill = guide_legend(override.aes = list(alpha = 1),nrow = 6), shape = guide_legend(nrow = 6)) +
  labs(x = "RGT score (%)",y = "Flexibility score (%)", title = "C")
pREV1

#Figure  2####
legd.teto <- legendGrob(c("Tph2-wt","Tph2-kd"), vgap= unit(.4,"lines"), gp=gpar(col  = col.teto, lwd=2,lineend=0,cex=.9))
legd.dm <- legendGrob(c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm"), 
                      pch = c(24,22,25,24,22,25), vgap = unit(.4,"lines"),gp = gpar(fill = col.dm, col = col.teto2,lwd = 2,cex = .9))

#save figure 2
legd = arrangeGrob(legd.teto,legd.dm, heights  = c(.1,.4))
line1 = arrangeGrob(grobs = list(p.RGT2, p.RGT1.1,p.RGT1.2), widths = c(1.8,.9,.4))
line1.3 =arrangeGrob(grobs = list(pREV1,legd), widths  = c(1.8,1.3))

f2<-grid.arrange(line1,line1.3,  nrow = 2)
# ggsave("F2-20240430.png", plot = f2, device="png", dpi=300,height = 17, width = 17,units = "cm")


#Probability discounting task ####

data.pdt <- filter(data3, normAUC > 0)
data.pdt.rgt <- data.pdt %>% filter(dm!= "OUT")#remove OUT in order to add RGT factor

#data.pdt<-filter(data, normAUC>0, batch%in%c(11,12,13,15,16,17,19))
nb.pdt <- data.pdt %>%
  group_by(group,dm) %>%
  summarise(Count = n()) %>%
  spread(group,Count) %>%
  filter(dm %in% c("GDM","PDM","INT"))
nb.pdt

PDT <- data.pdt %>% gather(proba,perf,PDT100:PDT9)
PDT <- PDT %>% mutate(proba=substr(proba, 4, 6))
PDT$proba <- factor(PDT$proba, levels=c("100","33","20","14","9"))

kinpdt <- PDT %>% group_by(proba, group) %>% summarise(Mean = mean(perf, na.rm = TRUE), 
                                                       SD = sd(perf, na.rm = TRUE), 
                                                       Median = median(perf, na.rm = TRUE), 
                                                       Quant5 = quantile(perf,c(.05), na.rm = TRUE), 
                                                       Quant95 = quantile(perf,c(.95), na.rm = TRUE))

kinpdtdm <- PDT %>% filter(dm!="OUT")%>%group_by(proba, groupdm,group)%>%summarise(Mean=mean(perf, na.rm=TRUE), 
                                                                                   SD=sd(perf, na.rm=TRUE), 
                                                                                   Median=median(perf, na.rm = TRUE), 
                                                                                   Quant5=quantile(perf,c(.05), na.rm=TRUE), 
                                                                                   Quant95=quantile(perf,c(.95), na.rm=TRUE))

#kinetics####
pd3 <- position_dodge(.15) # move them 1 to the left and right
p.pdt<-ggplot(kinpdtdm, aes(x=proba, y=Mean, color=groupdm, shape=groupdm, fill=groupdm)) +
  geom_hline(yintercept=50, colour="#333333", linetype="dashed", size=.3)+
  geom_errorbar(aes(ymin=Mean, ymax=Mean+SD), size=.3, width=0, position=pd3) +
  geom_line(aes(group=groupdm),size=.5) +
  geom_point(aes(col=groupdm), size=1.5,stroke=1)+#, shape=21, stroke=.5,color="black"
  theme_classic()+
  theme.teto+
  theme(legend.position ="none", plot.margin = margin(r = 0, t = 4, l = 4, b = 4, unit = "pt"))+# c(.8,.9)
  #theme(axis.title.x = element_text(size = 18,angle=0,face = 'plain',colour="black"))+
  scale_y_continuous(breaks=seq(0,100,10), limits = c(0,100), labels = c("0","", "20","", "40","", "60","", "80", "", "100"))+
  scale_x_discrete()+
  facet_grid(~group)+
  scale_color_manual(values = col.teto2) +
  scale_shape_manual(values = c(24,22,25,24,22,25))+
  scale_fill_manual(values = col.dm)+
  guides(fill = guide_legend(nrow=3), shape= guide_legend(nrow=3), col=guide_legend(nrow=3))+
  # labs(x="Probability (%)",y="Large reward choices (%)", title = "D")
  labs(x="Probability (%)",y="Large reward choices (%)", title = "A")
p.pdt

#AUC####
p.pdt.2<-ggplot(data = data.pdt.rgt, aes(x=groupdm, y=normAUC)) + 
  geom_boxplot(aes(colour=groupdm), outlier.shape = NA,alpha=0,show.legend = F,fatten=1.5, lwd=.5)+
  geom_beeswarm(aes(fill=groupdm, shape=groupdm, col=groupdm), size=1,stroke=.8, cex=3.5,  alpha=.6)+#stroke=0.5,groupOnX = T gives jitter to x axis#color='transparent',color="black",
  theme_classic()+
  theme.teto+
  theme(legend.position = "none", 
        plot.title = element_text(size = 11),
        axis.text.x = element_text(color = "transparent"),
        axis.ticks.x = element_line(color = "transparent"), 
        axis.title.x = element_text(color = "transparent"), 
        #axis.line.x = element_line(color = "transparent"),
        plot.margin = margin(r = 4, t = 4, l = 0, b = 4, unit = "pt"),)+
  scale_shape_manual(values = c(24,22,25,24,22,25))+
  scale_fill_manual(values = col.dm)+
  scale_color_manual(values = col.teto2)+
  labs(x="",y="", title="AUC")
p.pdt.2

lay.pdt<-rbind(c(1,1,1,2,2),c(1,1,1,3,3),c(1,1,1,3,3),c(1,1,1,3,3))
p.pdt.3<-grid.arrange(p.pdt,carre,p.pdt.2, layout_matrix = lay.pdt)

#anova on the performances
PDT1<-PDT%>%dplyr::select(group,groupdm,batch,rat,perf,proba,dm,RGT)%>%filter(dm!="OUT")

anova.p.pdt<-aovp(perf~proba*groupdm+Error(rat), data=PDT1,
                  perm = "Exact", Ca=0.01, maxIter=1000000)
sum.aovp3.1<-summary(anova.p.pdt)

anova.p.pdt<-aovp(perf~proba*group*RGT+Error(rat), data=PDT1,
                  perm = "Exact", Ca=0.01, maxIter=1000000)
sum.aovp3.12<-summary(anova.p.pdt)
#anova on the performances without Control-pdm (n=1)
PDT2<-PDT%>%dplyr::select(group,groupdm,batch,rat,perf,proba,dm,RGT)%>%filter(dm!="OUT", groupdm!="Control-pdm")

anova.p.pdt2<-aovp(perf~proba*group*RGT+Error(rat), data=PDT2,
                   perm = "Exact", Ca=0.01, maxIter=1000000)
sum.aovp3.2<-summary(anova.p.pdt2)

#anova on the performances without Control
PDT3<-PDT%>%dplyr::select(group,groupdm,batch,rat,perf,proba,dm,RGT)%>%filter(dm!="OUT", group!="Tph2-kd")

anova.p.pdt3<-aovp(perf~proba*groupdm+Error(rat), data=PDT3,
                   perm = "Exact", Ca=0.01, maxIter=1000000)
sum.aovp3.3<-summary(anova.p.pdt3)

#anova without PDM at all -> no difference between groups
PDT4<-PDT%>%dplyr::select(group,groupdm,batch,rat,perf,proba,dm,RGT)%>%filter(dm!="PDM")
anova.p.pdt4<-aovp(perf~proba*group+Error(rat), data=PDT4,
                   perm = "Exact", Ca=0.01, maxIter=1000000)
sum.aovp3.4<-summary(anova.p.pdt4)

#kruskal test on the AUC
stat.PDT1<-kruskal.test(data=data.pdt.rgt, normAUC~groupdm)#Kruskal-Wallis chi-squared = 13.02, df = 5, p-value = 0.0232

#kruskal test on the AUC without Control-pdm confirms the previous test
data.pdt.rgt.nocont.pdm<-data.pdt.rgt%>%filter(groupdm!="Control-pdm")
stat.PDT1.2<-kruskal.test(data=data.pdt.rgt.nocont.pdm, normAUC~groupdm)#Kruskal-Wallis chi-squared = 12.315, df = 4, p-value = 0.01516

#wilcox test on the AUC without PDM at all -> no difference
data.pdt.rgt.nopdm<-data.pdt.rgt%>%filter(dm!="PDM")
stat.PDT1.3<-wilcox.test(data=data.pdt.rgt.nopdm, normAUC~group)#W = 644, p-value = 0.06638

#Social recognition task ####

# data3<-data3%>%mutate(E1.norm = E1/hab*100,
#                       E2.norm = E2/hab*100,
#                       E3.norm = E3/hab*100)
# data3<-data3%>%mutate(STM.norm = E1.norm/E3.norm)#normalize data srt to hab

data.srt<-data3%>%filter(hab>0)#rat%in%c(207:266, 303:326)
nb.srt<-data.srt%>%
  group_by(group,dm) %>%
  summarise(Count = n()) %>%
  spread(group,Count) %>%
  filter(dm %in% c("GDM","PDM","INT"))
nb.srt
data.srt.rgt<-data.srt%>%filter(dm!="OUT")

SR<-data.srt%>%gather(time,perf,hab:E3)
SR$time <- factor(SR$time, levels=c("hab","E1","E2","E3"))#, labels = c(1,2,3,4)as.numeric(as.character(

SR1<-SR%>%dplyr::select(group,batch,rat,perf,time,dm,RGT, groupdm)
SR2<-SR1%>%filter(dm!="OUT")#remove OUT in order to add RGT factor in stats

#kinetics####
kinSRdm<- SR%>% filter(dm!="OUT")%>%group_by(time,groupdm,group)%>%summarise(Mean=mean(perf, na.rm=TRUE), SD=sd(perf, na.rm=TRUE), Median=median(perf, na.rm = TRUE), Quant5=quantile(perf,c(.05), na.rm=TRUE), Quant95=quantile(perf,c(.95), na.rm=TRUE))
pd4<-position_dodge(.15)
p.srt<-ggplot(kinSRdm, aes(x=time, y=Mean, color=groupdm, shape=groupdm, fill=groupdm)) +
  geom_errorbar(aes(ymin=Mean, ymax=Mean+SD), size=.3, width=0, position=pd4) +
  geom_line(aes(group=groupdm),size=.5) +
  geom_point(aes(col=groupdm), size=1.5,stroke=1)+#, shape=21, stroke=.5color="black"
  theme_classic()+
  theme.teto+
  theme(legend.position ="none",
        plot.margin = margin(r = 0, t = 4, l = 4, b = 4, unit = "pt"))+
  facet_grid(~group)+
  scale_color_manual(values = col.teto2) +
  scale_shape_manual(values = c(24,22,25,24,22,25))+
  scale_fill_manual(values = col.dm)+
  # labs(y="Interaction time (s)", title = "E", x="Encounters")
  labs(y="Interaction time (s)", title = "B", x="Encounters")
p.srt

#STM####
p.srt.2<-ggplot(data = data.srt.rgt, aes(x=groupdm, y=STM)) + 
  geom_hline(yintercept = 1, colour="#333333", linetype="dashed", linewidth=.3)+
  geom_boxplot(aes(colour=groupdm), outlier.shape = NA,alpha=0,show.legend = F,fatten=1.5, lwd=.5)+
  geom_beeswarm(aes(fill=groupdm, shape=groupdm, col=groupdm), size=1,stroke=.8, cex=3.5,  alpha=.6)+#stroke=0.5,groupOnX = T gives jitter to x axis#color='transparent',color="black",
  theme_classic()+
  theme.teto+
  theme(legend.position = "none", 
        plot.title = element_text(size = 11),
        axis.text.x = element_text(color = "transparent"),
        axis.ticks.x = element_line(color = "transparent"), 
        axis.title.x = element_text(color = "transparent"), 
        #axis.line.x = element_line(color = "transparent"),
        plot.margin = margin(r = 4, t = 4, l = 0, b = 4, unit = "pt"))+
  # scale_x_continuous(breaks=seq(0,100,10), labels = c("0","", "20","", "40","", "60","", "80", "", "100"))+#, limits = c(-0,103) #no limits to avoid loss of dots
  scale_shape_manual(values = c(24,22,25,24,22,25))+
  scale_fill_manual(values = col.dm)+
  scale_color_manual(values = col.teto2)+
  # facet_grid(~group, scale = "free_x")+
  # guides(fill = guide_legend(override.aes=list(alpha=1),nrow=6), shape= guide_legend(nrow=6))+
  labs(x="",y="", title = "STM")
p.srt.2

lay.srt<-rbind(c(1,1,1,1,2,2,2),c(1,1,1,1,3,3,3),c(1,1,1,1,3,3,3),c(1,1,1,1,3,3,3))
p.srt.3<-grid.arrange(p.srt,carre,p.srt.2, layout_matrix = lay.srt)

#anova on the performances 
anova.p.srt<-aovp(perf~time*group*RGT+Error(rat), data=SR2,
                  perm = "Exact", Ca=0.01, maxIter=1000000)
sum.aovp4.1<-summary(anova.p.srt)

# Error: rat
# Component 1 :
#   Df R Sum Sq R Mean Sq   Iter Pr(Prob)  
# RGT        1        8       7.6  19920   0.3342  
# group      1     2360    2359.8  33602   0.2294  
# RGT:group  1     8042    8042.3 151279   0.0620 .
# Residuals 78   177547    2276.2                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# Error: Within
# Component 1 :
#   Df R Sum Sq R Mean Sq    Iter Pr(Prob)    
# time             3   537934    179311 1000000  < 2e-16 ***
#   RGT:time         3     1015       338   12529  0.70325    
# time:group       3      711       237   17691  0.59997    
# RGT:time:group   3     8037      2679  780570  0.01802 *  
#   Residuals      234   207603       887                     

#difference of STM ratio from 1
data.srt.rgt.cont<-data.srt%>%filter(dm!="OUT",group=="Control")
data.srt.rgt.kd<-data.srt%>%filter(dm!="OUT",group=="Tph2-kd")

srt.gdmkd<-data.srt.rgt.kd%>%filter(groupdm=="Tph2-kd-gdm")
srt.intkd<-data.srt.rgt.kd%>%filter(groupdm=="Tph2-kd-int")
srt.pdmkd<-data.srt.rgt.kd%>%filter(groupdm=="Tph2-kd-pdm")
wilcox.test(srt.gdmkd$STM, mu=1)#V = 315.5, p-value = 4.063e-05
wilcox.test(srt.intkd$STM, mu=1)#V = 45, p-value = 0.003906
wilcox.test(srt.pdmkd$STM, mu=1)#V = 9, p-value = 0.7865

srt.gdmcont<-data.srt.rgt.cont%>%filter(groupdm=="Control-gdm")
srt.intcont<-data.srt.rgt.cont%>%filter(groupdm=="Control-int")
srt.pdmcont<-data.srt.rgt.cont%>%filter(groupdm=="Control-pdm")
wilcox.test(srt.gdmcont$STM, mu=1)#V = 570, p-value = 3.311e-06
wilcox.test(srt.intcont$STM, mu=1)#V = 28, p-value = 0.01563
#wilcox.test(srt.pdmcont$STM, mu=1)#do not run n=1 here


#Figure 3 ####

legd.dm2 <- legendGrob(c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm"), 
                       pch = c(24,22,25,24,22,25), vgap = unit(.4,"lines"),gp = gpar(fill = col.dm, col = col.teto2,lwd = 2,cex = .7))
carre=grobTree(rectGrob(gp = gpar(fill = "white", alpha=0)))

#save figure 3
line2.3 =arrangeGrob(grobs = list(carre,legd.dm,carre,carre), heights  = c(.5,1,1,1))
line2 = arrangeGrob(grobs = list(p.pdt.3, p.srt.3), heights = c(1,1))
lay3 = rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,3,3,3,2))
f3<-grid.arrange(line2,carre, line2.3, nrow = 2, layout_matrix = lay3)
#ggsave("F3-20240430.png", plot = f3, device="png", dpi=300,height = 14, width = 17,units = "cm")


#Odor discrimination fig S2####

data.od <- data3 %>% filter(odor.preference>0)

odor1<-ggplot(data = data.od, aes(x=group, y=odor.preference)) +  
  geom_hline(yintercept=.5, colour="#333333", linetype="dashed", size=.3)+
  geom_boxplot(aes(col=group), outlier.shape = NA,alpha=0,show.legend = F,fatten=1.5, lwd=.5)+
  geom_beeswarm(aes(fill=groupdm, shape=groupdm, col=group), size=2, stroke=1, alpha=.8, cex = 3.5)+#stroke=0.5,groupOnX = T gives jitter to x axis#color='transparent',color="black",
  # geom_boxplot(aes(color= group),show.legend = F, lwd=0.3)+#outlier.shape = NA, remove outliers in background 
  # geom_segment(x=1, xend=2, y=103,yend=103, colour="black", size=1)+
  # geom_beeswarm(aes(fill=groupdm, shape=groupdm), size=3,alpha=.8, color="black", cex = 3.5, priority = "random",  show.legend = T)+#stroke=.5, groupOnX = T,
  theme_classic()+
  theme.teto+
  theme(legend.position = "bottom")+  
  scale_shape_manual(values = c(24,22,25,24,22,25),labels = c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm"))+
  scale_fill_manual(values = col.dm, labels = c("Tph2-wt-gdm","Tph2-wt-int","Tph2-wt-pdm","Tph2-kd-gdm","Tph2-kd-int","Tph2-kd-pdm"))+
  scale_color_manual(values = col.teto)+
  scale_x_discrete(labels= c("Tph2-wt","Tph2-kd")) +
  # guides(fill = guide_legend(override.aes=list(alpha=1),nrow=6), shape= guide_legend(nrow=6))+
  #guides(fill=T,shape=T)+
  # guides(shape = guide_legend(override.aes = list(alpha=1, fill=c("#331068FF","#331068FF","#331068FF","#DE4968FF", "#DE4968FF","#DE4968FF")),nrow=3), fill=F)+#color manual to avoid mix...
  guides(fill = guide_legend(override.aes=list(alpha=1,size=3, col=c("#331068FF","#331068FF","#331068FF","#DE4968FF","#DE4968FF","#DE4968FF")),nrow=3),
         shape = guide_legend(nrow=3),
         col = "none")+
  # guides(fill = guide_legend(override.aes=list(alpha=1, col=c("#331068FF","#331068FF","#331068FF","#DE4968FF","#DE4968FF","#DE4968FF")),nrow=6),
  #        shape = guide_legend(nrow=6),
  #        col = "none")+
  # scale_y_continuous(limits = c(0,1))+
  labs(y="Preference for social odor", title = "")
odor1

wilcox.test(odor.preference~group, data=data.od)#W = 417, p-value = 0.7923

anova.p.odor<-aovp(odor.preference~group*RGT+Error(rat), data=data.od,
                   perm = "Exact", Ca=0.01, maxIter=1000000)
summary(anova.p.odor)

#Save fig S2####
#ggsave("S2_20240430.png", plot = odor1, device="png", dpi=600,height = 13, width = 9,units = "cm")

#Network Analysis####
library(Hmisc)
library(corpcor)
library(qgraph)
library(NetworkComparisonTest)
library(igraph)

#select data
data.network<-data3%>%dplyr::select("rat", "batch","genotype","treatment","group","groupdm", "groupx","dm","pdm_rest","RGT", "latency","normAUC","Soc.pref","STM")%>%
  #flexibility_index not included
  filter(dm!="OUT")#remove dm="OUT"

#rename variables
colnames(data.network)=c("rat", "batch","genotype","treatment","group","groupdm", "groupx","dm","pdm_rest","RGT", "Lat","PDT","SP","STM")

#remove lines with NAs
data.network.2 = data.network[complete.cases(data.network),]

#filter groups
data.net_control<-data.network.2%>%filter(group=="Control")
data.net_tetodox<-data.network.2%>%filter(group=="Tph2-kd")

data.net_tetodoxpdm<-data.network.2%>%filter(groupdm=="Tph2-kd-pdm")
data.net_controlnopdm<-data.net_control%>%filter(groupdm!="Control-pdm")
data.net_tetodoxnopdm<-data.net_tetodox%>%filter(groupdm!="Tph2-kd-pdm")

#check normality - most variables are not normally distributed -
shapiro.test(data.net_control[,10])#RGT
shapiro.test(data.net_control[,11])#Lat yes
shapiro.test(data.net_control[,12])#PDT yes
shapiro.test(data.net_control[,13])#SP
shapiro.test(data.net_control[,14])#STM

shapiro.test(data.net_controlnopdm[,10])#RGT
shapiro.test(data.net_controlnopdm[,11])#Lat yes
shapiro.test(data.net_controlnopdm[,12])#PDT yes
shapiro.test(data.net_controlnopdm[,13])#SP
shapiro.test(data.net_controlnopdm[,14])#STM

shapiro.test(data.net_tetodox[,10])#RGT
shapiro.test(data.net_tetodox[,11])#Lat
shapiro.test(data.net_tetodox[,12])#PDT yes
shapiro.test(data.net_tetodox[,13])#SP yes
shapiro.test(data.net_tetodox[,14])#STM

shapiro.test(data.net_tetodoxnopdm[,10])#RGT
shapiro.test(data.net_tetodoxnopdm[,11])#Lat yes
shapiro.test(data.net_tetodoxnopdm[,12])#PDT yes
shapiro.test(data.net_tetodoxnopdm[,13])#SP
shapiro.test(data.net_tetodoxnopdm[,14])#STM

shapiro.test(data.net_tetodoxpdm[,10])#RGT 
shapiro.test(data.net_tetodoxpdm[,11])#Lat yes
shapiro.test(data.net_tetodoxpdm[,12])#PDT yes
shapiro.test(data.net_tetodoxpdm[,13])#SP yes
shapiro.test(data.net_tetodoxpdm[,14])#STM 

#check correlations Spearman
#Pearson's correlation are the classical interaction measures in network analysis
#but our data are not normal so we have to use Spearman's correlations

cont<-as.matrix(data.net_controlnopdm[,10:14])#
R.cont<-rcorr(cont,type=c("spearman"))

net.cont<-qgraph(R.cont[["r"]],layout = "circular", threshold = 0.28, 
                 minimum = 0, maximum = .7, 
                 labels = c("RGT","Lat","PDT","SP","STM"),node.width = 1.5, edge.labels = TRUE, 
                 edge.label.color = "black", edge.label.bg = "white", edge.label.cex=3, 
                 sampleSize = "maximum",title = "A", title.cex = 3, 
                 bg = "white", filename = "")

teto<-as.matrix(data.net_tetodoxnopdm[,10:14])#
R.kd<-rcorr(teto,type=c("spearman"))

net.kd<-qgraph(R.kd[["r"]],layout="circular",threshold=0.28, 
               minimum=0,maximum=.7, 
               labels =c("RGT","Lat","PDT","SP","STM"), node.width = 1.5, edge.labels=TRUE, 
               edge.label.color="black",edge.label.bg="white", edge.label.cex=3,
               sampleSize="maximum",title="B",title.cex=3, 
               bg="white",filename="")

tetopdm<-as.matrix(data.net_tetodoxpdm[,10:14])#
R.kd.pdm<-rcorr(tetopdm,type=c("spearman"))

net.pdm<-qgraph(R.kd.pdm[["r"]],layout="circular",threshold=0.28, 
                minimum=0,maximum=.7, 
                labels =c("RGT","Lat","PDT","SP","STM"), node.width = 1.5, edge.labels=TRUE, 
                edge.label.color="black",edge.label.bg="white", edge.label.cex=3, 
                sampleSize="maximum",title="C",title.cex=3, 
                bg="white",filename="")

#Saving the images
#cont-NO-pdm
png(filename = "net.tph2-wt_nopdm.png", width = 480, height = 480)
plot(net.cont)
dev.off()
#teto-NO-pdm
png(filename = "net.tph2-kd_nopdm.png", width = 480, height = 480)
plot(net.kd)
dev.off()
#teto-dox-pdm
png(filename = "net.tph2-kd_pdm.png", width = 480, height = 480)
plot(net.pdm)
dev.off()

#Combining image files into one image
library(gridExtra)
library(grid)
library(png)

img1 <- readPNG("./net.tph2-wt_nopdm.png")
img2 <- readPNG("./net.tph2-kd_nopdm.png")
img3 <- readPNG("./net.tph2-kd_pdm.png")

img1 =rasterGrob(img1)
img2 =rasterGrob(img2)
img3 =rasterGrob(img3)

#Figure 4####
f4.net<-grid.arrange(img1,img2,img3, ncol=3)

#saving figure 4
#ggsave("plot.net.spearman_30042024.png", plot = f4.net, device="png", dpi=600,height = 4, width = 12,units = "cm")
