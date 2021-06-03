
#R code accompanying Nätscher et al. "Morphological response accompanying 
#size reduction of belemnites during an Early Jurassic hyperthermal event 
#modulated by life history"


# load packages -----------------------------------------------------------

library(geomorph)
library(nlme)
library(vegan)
library(grDevices)
library(effsize)
library(car)
library(AICcmodavg)
library(raster)
library(dplyr)
library(corrplot)
library(gdata)
library(nlme)
library(RVAideMemoire)
library(rcompanion)

pal.t3.trans <- c("gray90", "gray65", "gray24")
pal.t8 <- c("#ffe600", "#95eb00", "#00b530", "#00bda7", "#0083c9", "#0003c9", "#6a00b5", "#a8006d")
pal.t8.trans <- c("#fff27a", "#b7e070", "#6cb880", "#79bab3", "#63a6c9", "#5b5dc2", "#8f5db3", "#b06395")

# GPA ---------------------------------------------------------------------

proto.gpa <- gpagen(proto, curves = curves)
plot(proto.gpa)

plot3d(proto.gpa$consensus, xlim=c(-0.5,1), ylim=c(-0.5,1), zlim=c(-0.5,1), size=10, box=F, axes=F, xlab="", ylab="", zlab="")


# PCA ---------------------------------------------------------------------

proto.pca <- plotTangentSpace(proto.gpa$coords, label=TRUE, verbose=TRUE, groups = BelDat$spec.fac)
proto.pca$pc.summary
plot3d(proto.pca$pc.scores[,1:3],type="s",pch=20,col=BelDat$spec.fac, size=2,xlab="PC1",
       ylab="PC2",zlab="PC3",xlim= c(-0.15,0.23), ylim= c(-0.15,0.23), zlim= c(-0.15,0.23))


# plot PC Axes ------------------------------------------------------------

#PC1.min 
plotRefToTarget(proto.pca$pc.shapes$PC1min, proto.gpa$consensus,
                                 method="points", links=links, gridPars=designy)
#PC1.max 
plotRefToTarget(proto.pca$pc.shapes$PC1max, proto.gpa$consensus,
                                 method="points", links= links, gridPars=designy)

#PC2.min
plotRefToTarget(proto.pca$pc.shapes$PC2min, proto.gpa$consensus,
                                 method="points", links=links, gridPars=designy)
#PC2.max
plotRefToTarget(proto.pca$pc.shapes$PC2max, proto.gpa$consensus,
                                 method="points", links=links, gridPars=designy)



# SIZE ~ MORPHOLOGY -------------------------------------------------------

# * Geometric mean vs. Centroid Size --------------------------------------
cor(BelDat$GM, proto.gpa$Csize, use = "pairwise.complete.obs")

lmod <- lm(BelDat$GM ~ proto.gpa$Csize)
summary(lmod)

# * Centroid size vs PC1 --------------------------------------------------
cor(proto.gpa$Csize, BelDat$PC1, use = "pairwise.complete.obs")

lin2.mod <- lm(proto.gpa$Csize ~ BelDat$PC1)
summary(lin2.mod)



# EFFECT SIZES (hedge's g) ------------------------------------------------

#PC1  
#all belemnites
cohen.d(BelDat$PC1[sub1],BelDat$PC1[sub2],hedges.correction=TRUE)
cohen.d(BelDat$PC1[sub2],BelDat$PC1[sub3],hedges.correction=TRUE)
cohen.d(BelDat$PC1[sub3],BelDat$PC1[sub4],hedges.correction=TRUE)
cohen.d(BelDat$PC1[sub4],BelDat$PC1[sub5],hedges.correction=TRUE)

#P bisulcata
cohen.d(BelDat$PC1[sub1.bis],BelDat$PC1[sub2.bis],hedges.correction=TRUE)
cohen.d(BelDat$PC1[sub2.bis],BelDat$PC1[sub3.bis],hedges.correction=TRUE)
cohen.d(BelDat$PC1[sub3.bis],BelDat$PC1[sub4.bis],hedges.correction=TRUE)

#C longiforma
cohen.d(BelDat$PC1[sub2.long],BelDat$PC1[sub3.long],hedges.correction=TRUE)
cohen.d(BelDat$PC1[sub3.long],BelDat$PC1[sub4.long],hedges.correction=TRUE)


#for csize
#all belemnites
cohen.d(BelDat$proto.gpa.Csize[sub1],BelDat$proto.gpa.Csize[sub2],hedges.correction=TRUE)
cohen.d(BelDat$proto.gpa.Csize[sub2],BelDat$proto.gpa.Csize[sub3],hedges.correction=TRUE)
cohen.d(BelDat$proto.gpa.Csize[sub3],BelDat$proto.gpa.Csize[sub4],hedges.correction=TRUE)
cohen.d(BelDat$proto.gpa.Csize[sub4],BelDat$proto.gpa.Csize[sub5],hedges.correction=TRUE)

#P bisulcata
cohen.d(BelDat$proto.gpa.Csize[sub1.bis],BelDat$proto.gpa.Csize[sub2.bis],hedges.correction=TRUE)
cohen.d(BelDat$proto.gpa.Csize[sub2.bis],BelDat$proto.gpa.Csize[sub3.bis],hedges.correction=TRUE)
cohen.d(BelDat$proto.gpa.Csize[sub3.bis],BelDat$proto.gpa.Csize[sub4.bis],hedges.correction=TRUE)

#C longiforma
cohen.d(BelDat$proto.gpa.Csize[sub2.long],BelDat$proto.gpa.Csize[sub3.long],hedges.correction=TRUE)
cohen.d(BelDat$proto.gpa.Csize[sub3.long],BelDat$proto.gpa.Csize[sub4.long],hedges.correction=TRUE)



# SPECIES COMPOSITION (subsampling method) --------------------------------

boot <- list()

#number of specimens in source pool
boot[[1]] <- nrow(BelDat)
#number of specimens of each species
boot[[2]] <- c(nrow(BelDat[BelDat$species == "bisulcata",]), nrow(BelDat[BelDat$species == "longiformis",]), 
               nrow(BelDat[BelDat$species == "Hastitidae sp. indet.",]), nrow(BelDat[BelDat$species == "Parapassaloteuthis sp. 1",]), 
               nrow(BelDat[BelDat$species == "milleri",]), nrow(BelDat[BelDat$species == "Bairstowius sp. A",]), 
               nrow(BelDat[BelDat$species == "Passaloteuthis sp. juv",]), nrow(BelDat[BelDat$species == "Acrocoelites sp. ",]))
#relative overall abundance of each species
boot[[3]] <- boot[[2]] / boot[[1]]
#names of species
boot[[4]] <- c("P. bisulcata", "C. longiforma", "Hastitidae sp.", "Parapassaloteuthis sp.", "P. milleri", 
               "Bairstowius sp.", "Passaloteuthis sp. juv", "Acrocoelites sp.")


boot[[5]] <- as.data.frame(matrix(nrow=29, ncol=500)) #29 = mean sample size of the 5 ammonite subzones

for (i in 1:500){
        boot[[5]][,i] <- sample(boot[[4]], size = 29, replace = TRUE, prob = boot[[3]])
}

#relative abundances for all of these bootstrapped communities
boot[[6]] <- as.data.frame(matrix(nrow = 8, ncol = 500))

for (i in 1:500){
        boot[[6]][1,i] <- length(boot[[5]][,i][boot[[5]][,i] == "P. bisulcata"]) / 29
        boot[[6]][2,i] <- length(boot[[5]][,i][boot[[5]][,i] == "C. longiforma"]) / 29
        boot[[6]][3,i] <- length(boot[[5]][,i][boot[[5]][,i] == "Hastitidae sp."]) / 29
        boot[[6]][4,i] <- length(boot[[5]][,i][boot[[5]][,i] == "Parapassaloteuthis sp."]) / 29
        boot[[6]][5,i] <- length(boot[[5]][,i][boot[[5]][,i] == "P. milleri"]) / 29
        boot[[6]][6,i] <- length(boot[[5]][,i][boot[[5]][,i] == "Bairstowius sp."]) / 29
        boot[[6]][7,i] <- length(boot[[5]][,i][boot[[5]][,i] == "Passaloteuthis sp. juv"]) / 29
        boot[[6]][8,i] <- length(boot[[5]][,i][boot[[5]][,i] == "Acrocoelites sp."]) / 29
        
}

boot[[7]] <- as.data.frame(matrix(nrow = 8, ncol = 6))
colnames(boot[[7]]) <- c("mean", "sd", "lowerCI", "upperCI", "min", "max")

for (i in 1:8){
        
        boot[[7]][i,1] <- mean(as.vector(as.numeric(boot[[6]][i,])))
        boot[[7]][i,2] <- sd(as.vector(as.numeric(boot[[6]][i,])))
        boot[[7]][i,3] <- mean(as.vector(as.numeric(boot[[6]][i,]))) - qnorm(0.975)*sd(as.vector(as.numeric(boot[[6]][i,])))/sqrt(500)
        boot[[7]][i,4] <- mean(as.vector(as.numeric(boot[[6]][i,]))) + qnorm(0.975)*sd(as.vector(as.numeric(boot[[6]][i,])))/sqrt(500)
        boot[[7]][i,5] <- min(as.vector(as.numeric(boot[[6]][i,])))
        boot[[7]][i,6] <- max(as.vector(as.numeric(boot[[6]][i,])))
}


# WITHIN-SPECIES ----------------------------------------------------------

#testing for change in size
#P.bisulcata
wilcox.test(BelDat$proto.gpa.Csize[sub1.bis], BelDat$proto.gpa.Csize[sub2.bis])
wilcox.test(BelDat$proto.gpa.Csize[sub2.bis], BelDat$proto.gpa.Csize[sub3.bis])
wilcox.test(BelDat$proto.gpa.Csize[sub3.bis], BelDat$proto.gpa.Csize[sub4.bis])

#C.longiforma
wilcox.test(BelDat$proto.gpa.Csize[sub2.long], BelDat$proto.gpa.Csize[sub3.long])
wilcox.test(BelDat$proto.gpa.Csize[sub3.long], BelDat$proto.gpa.Csize[sub4.long])

#assemblage
wilcox.test(BelDat$proto.gpa.Csize[sub1], BelDat$proto.gpa.Csize[sub2])
wilcox.test(BelDat$proto.gpa.Csize[sub2], BelDat$proto.gpa.Csize[sub3])
wilcox.test(BelDat$proto.gpa.Csize[sub3], BelDat$proto.gpa.Csize[sub4])


#testing for change in PC1
#P.bisulcata
wilcox.test(BelDat$PC1[sub1.bis], BelDat$PC1[sub2.bis])
wilcox.test(BelDat$PC1[sub2.bis], BelDat$PC1[sub3.bis])
wilcox.test(BelDat$PC1[sub3.bis], BelDat$PC1[sub4.bis])

#C.longiforma
wilcox.test(BelDat$PC1[sub2.long], BelDat$PC1[sub3.long])
wilcox.test(BelDat$PC1[sub3.long], BelDat$PC1[sub4.long])

#assemblage
wilcox.test(BelDat$PC1[sub1], BelDat$PC1[sub2])
wilcox.test(BelDat$PC1[sub2], BelDat$PC1[sub3])
wilcox.test(BelDat$PC1[sub3], BelDat$PC1[sub4])


# REGRESSIONS CSIZE/PC1 BY TEMPERATURE ------------------------------------

#C. longiforma
#by temperature
windows(w=9, h=6)
plot(proto.gpa$Csize[long], BelDat$PC1[long],   
     ylab="PC1", xlab="Centroid Size [mm]", col = BelDat$dev.fac[long], 
     pch = c(15,0)[BelDat$T.class[long]], cex=1.5)
#low = hollow, high = filled
abline(lm(BelDat$PC1[long[long %in% lowT]] ~ proto.gpa$Csize[long[long %in% lowT]]), lty = 1, lwd = 3, col="#3651c9")#punkte
abline(lm(BelDat$PC1[long[long %in% highT]] ~ proto.gpa$Csize[long[long %in% highT]]), lty = 1, lwd = 3, col="#b82e2e")#striche
legend("topright", legend=c("Low temperature", "High temperature","adult","neanic","juvenile"), 
       col=c("#3651c9", "#b82e2e", pal.t3.trans[3:1]), lty=c(1,1,0,0,0), lwd=2.5,
       pch=c(0,15,15,15,15),pt.cex=1.9,cex=1.2,bty="n")


lm.CL_lo <- lm(BelDat$PC1[long[long %in% lowT]] ~ proto.gpa$Csize[long[long %in% lowT]])
summary(lm.CL_lo)
lm.CL_hi <- lm(BelDat$PC1[long[long %in% highT]] ~ proto.gpa$Csize[long[long %in% highT]])
summary(lm.CL_hi)


#P.bisulcata
#by Temperature
windows(w=9, h=6)
plot(BelDat$proto.gpa.Csize[bis], BelDat$PC1[bis],   
     ylab="PC1", xlab="Centroid Size [mm]", col = BelDat$dev.fac[bis], 
     pch = c(17,2)[BelDat$T.class[bis]], cex=1.5)
abline(lm(BelDat$PC1[bis[bis %in% lowT]] ~ proto.gpa$Csize[bis[bis %in% lowT]]), lty = 1, lwd = 3, col="#3651c9")
abline(lm(BelDat$PC1[bis[bis %in% highT]] ~ proto.gpa$Csize[bis[bis %in% highT]]), lty = 1, lwd = 3, col="#b82e2e")
legend("topright", legend=c("Low temperature", "High temperature","adult","neanic","juvenile"), 
       col=c("#3651c9", "#b82e2e", pal.t3.trans[3:1]), lty=c(1,1,0,0,0), lwd=2.5,
       pch=c(2,17,15,15,15),pt.cex=1.9,cex=1.2,bty="n")

lm.PB_lo <- lm(BelDat$PC1[bis[bis %in% lowT]] ~ proto.gpa$Csize[bis[bis %in% lowT]])
summary(lm.PB_lo)
lm.PB_hi <- lm(BelDat$PC1[bis[bis %in% highT]] ~ proto.gpa$Csize[bis[bis %in% highT]])
summary(lm.PB_hi)


# VARIATION PARTITIONING --------------------------------------------------

varpart <- varpart(BelDat$PC1, ~ BelDat$d18O, ~ BelDat$d13C, ~ as.factor(BelDat$lit.fac), ~ BelDat$abundance..no..Of.belemnites.m2.,
                   data = BelDat, scale = T)
varpart

windows(h=7, w=6)
barplot(varpart$part$indfract$Adj.R.square[1:4],
        ylim=c(0,0.08),
        names=c("d18O", "d13C", "lithology", "abuncance"),
        col=c("dodgerblue2", "dodgerblue2", "tomato2", "tomato2"), border= F,
        ylab="Indepentent effect (R^2)")

anova(rda(BelDat$PC1 ~ BelDat$d18O + BelDat$d13C + as.factor(BelDat$lit.fac) + BelDat$abundance..no..Of.belemnites.m2.,
          data = BelDat, scale = T))
RsquareAdj(rda(BelDat$PC1 ~ BelDat$d18O + BelDat$d13C + BelDat$lit.fac + BelDat$abundance..no..Of.belemnites.m2.,
               data = BelDat, scale = T))

anova(rda(BelDat$PC1 ~ BelDat$d18O,
          data = BelDat, scale = T)) # p = 0.036 *, adjR2 = 0.02
RsquareAdj(rda(BelDat$PC1 ~ BelDat$d18O,
               data = BelDat, scale = T))

anova(rda(BelDat$PC1 ~ BelDat$d13C,
          data = BelDat, scale = T)) # p = 0.382, adjR2 = -0.0021
RsquareAdj(rda(BelDat$PC1 ~ BelDat$d13C,
               data = BelDat, scale = T))

anova(rda(BelDat$PC1 ~ as.factor(BelDat$lit.fac),
          data = BelDat, scale = T)) # p = 0.001 ***, adjR2 = 0.08
RsquareAdj(rda(BelDat$PC1 ~ as.factor(BelDat$lit.fac),
               data = BelDat, scale = T))

anova(rda(BelDat$PC1 ~ BelDat$abundance..no..Of.belemnites.m2.,
          data = BelDat, scale = T)) # p = 0.781, adjR2 = -0.00649
RsquareAdj(rda(BelDat$PC1 ~ BelDat$abundance..no..Of.belemnites.m2.,
               data = BelDat, scale = T))


# GENERALISED LEAST SQUARES MODELS -----------------------------------------

# Whole assemblage --------------------------------------------------------
# * correct for lithology -------------------------------------------------

correct_mod <- lm(MOD.dat$PC1 ~ as.factor(MOD.dat$lithology))
MOD.dat$Residuals <- correct_mod$resid #add residuals to dataframe

#correlation between the variables
corr <- cor(MOD.dat[, c(2:4)], method = "spearman",
            use = "pairwise.complete.obs")

windows(h=7,w=7)
corrplot(corr, method = "number", type="upper", order="hclust", tl.col = "black")

# Assemblage level
lm.PC1 <- lm(Residuals ~ d13C + d18O + Hg.TOC,
             data = MOD.dat, na.action = "na.omit")
acf(residuals(lm.PC1)) #test for auto-correlation of residuals
durbinWatsonTest(lm.PC1, max.lag=10) #lag at 1 and 2 significant


gls.c.1 <- gls(Residuals ~ 1, correlation = corARMA(p=1), data = MOD.dat,  method="ML", na.action = "na.omit")
summary(gls.c.1)
gls.c.2 <- gls(Residuals ~ d18O + d13C + Hg.TOC, correlation = corARMA(p=1), data = MOD.dat,  method="ML", na.action = "na.omit")
summary(gls.c.2)
gls.c.3 <- gls(Residuals ~ d18O + d13C, correlation = corARMA(p=1), data = MOD.dat,  method="ML", na.action = "na.omit")
summary(gls.c.3)
gls.c.4 <- gls(Residuals ~ d18O + Hg.TOC, correlation = corARMA(p=1), data = MOD.dat,  method="ML", na.action = "na.omit")
summary(gls.c.4)
gls.c.5 <- gls(Residuals ~ d13C + Hg.TOC, correlation = corARMA(p=1), data = MOD.dat,  method="ML", na.action = "na.omit")
summary(gls.c.5)
gls.c.6 <- gls(Residuals ~ d13C, correlation = corARMA(p=1), data = MOD.dat,  method="ML", na.action = "na.omit")
summary(gls.c.6)
gls.c.7 <- gls(Residuals ~ d18O, correlation = corARMA(p=1), data = MOD.dat,  method="ML", na.action = "na.omit")
summary(gls.c.7)
gls.c.8 <- gls(Residuals ~ Hg.TOC, correlation = corARMA(p=1), data = MOD.dat,  method="ML", na.action = "na.omit")
summary(gls.c.8)

AICcscores <- c(AICc(gls.c.1), AICc(gls.c.2), AICc(gls.c.3), AICc(gls.c.4), AICc(gls.c.5), 
                AICc(gls.c.6), AICc(gls.c.7), AICc(gls.c.8))
AICcscores #Display the AICc scores


# P. bisulcata ------------------------------------------------------------
# * correct for lithology -------------------------------------------------

#MOD.dat_bis <- MOD.dat[bis,]

#check if the residuals are autocorrelated
lm.PC1_bis <- lm(Residuals ~ d13C + d18O + Hg.TOC,
                 data = MOD.dat_bis, na.action = "na.omit")
acf(residuals(lm.PC1_bis)) #test for auto-correlation of residuals
durbinWatsonTest(lm.PC1_bis, max.lag=10) #lag at 1 and 2 significant


gls.c.1_bis <- gls(Residuals ~ 1, correlation = corARMA(p=1), data = MOD.dat_bis,  method="ML", na.action = "na.omit")
summary(gls.c.1_bis)
gls.c.2_bis <- gls(Residuals ~ d18O + d13C + Hg.TOC, correlation = corARMA(p=1), data = MOD.dat_bis,  method="ML", na.action = "na.omit")
summary(gls.c.2_bis)
gls.c.3_bis <- gls(Residuals ~ d18O + d13C, correlation = corARMA(p=1), data = MOD.dat_bis,  method="ML", na.action = "na.omit")
summary(gls.c.3_bis)
gls.c.4_bis <- gls(Residuals ~ d18O + Hg.TOC, correlation = corARMA(p=1), data = MOD.dat_bis,  method="ML", na.action = "na.omit")
summary(gls.c.4_bis)
gls.c.5_bis <- gls(Residuals ~ d13C + Hg.TOC, correlation = corARMA(p=1), data = MOD.dat_bis,  method="ML", na.action = "na.omit")
summary(gls.c.5_bis)
gls.c.6_bis <- gls(Residuals ~ d13C, correlation = corARMA(p=1), data = MOD.dat_bis,  method="ML", na.action = "na.omit")
summary(gls.c.6_bis)
gls.c.7_bis <- gls(Residuals ~ d18O, correlation = corARMA(p=1), data = MOD.dat_bis,  method="ML", na.action = "na.omit")
summary(gls.c.7_bis)
gls.c.8_bis <- gls(Residuals ~ Hg.TOC, correlation = corARMA(p=1), data = MOD.dat_bis,  method="ML", na.action = "na.omit")
summary(gls.c.8_bis)

AICcscores_bis <- c(AICc(gls.c.1_bis), AICc(gls.c.2_bis), AICc(gls.c.3_bis), AICc(gls.c.4_bis), AICc(gls.c.5_bis), 
                    AICc(gls.c.6_bis), AICc(gls.c.7_bis), AICc(gls.c.8_bis))
AICcscores_bis #Display the AICc scores

# Deriving a global p-value using ANOVA for the model no. 5, the best one for the assemblage scale:
anova(gls.c.1_bis, gls.c.7_bis) 


# C. longiforma -----------------------------------------------------------
# * correct for lithology -------------------------------------------------

#MOD.dat_long <- MOD.dat[long,]

#check if the residuals are autocorrelated
lm.PC1_long <- lm(Residuals ~ d13C + d18O + Hg.TOC,
                  data = MOD.dat_long, na.action = "na.omit")
windows(h=7, w=8)
acf(residuals(lm.PC1_long)) #test for auto-correlation of residuals
durbinWatsonTest(lm.PC1_long, max.lag=10) #lag at 1 and 2 significant


gls.c.1_long <- gls(Residuals ~ 1, correlation = corARMA(p=1), data = MOD.dat_long,  method="ML", na.action = "na.omit")
summary(gls.c.1_long)
gls.c.2_long <- gls(Residuals ~ d18O + d13C + Hg.TOC, correlation = corARMA(p=1), data = MOD.dat_long,  method="ML", na.action = "na.omit")
summary(gls.c.2_long)
gls.c.3_long <- gls(Residuals ~ d18O + d13C, correlation = corARMA(p=1), data = MOD.dat_long,  method="ML", na.action = "na.omit")
summary(gls.c.3_long)
gls.c.4_long <- gls(Residuals ~ d18O + Hg.TOC, correlation = corARMA(p=1), data = MOD.dat_long,  method="ML", na.action = "na.omit")
summary(gls.c.4_long)
gls.c.5_long <- gls(Residuals ~ d13C + Hg.TOC, correlation = corARMA(p=1), data = MOD.dat_long,  method="ML", na.action = "na.omit")
summary(gls.c.5_long)
gls.c.6_long <- gls(Residuals ~ d13C, correlation = corARMA(p=1), data = MOD.dat_long,  method="ML", na.action = "na.omit")
summary(gls.c.6_long)
gls.c.7_long <- gls(Residuals ~ d18O, correlation = corARMA(p=1), data = MOD.dat_long,  method="ML", na.action = "na.omit")
summary(gls.c.7_long)
gls.c.8_long <- gls(Residuals ~ Hg.TOC, correlation = corARMA(p=1), data = MOD.dat_long,  method="ML", na.action = "na.omit")
summary(gls.c.8_long)

AICcscores_long <- c(AICc(gls.c.1_long), AICc(gls.c.2_long), AICc(gls.c.3_long), AICc(gls.c.4_long), AICc(gls.c.5_long), 
                     AICc(gls.c.6_long), AICc(gls.c.7_long), AICc(gls.c.8_long))
AICcscores_long #Display the AICc scores

# Deriving a global p-value using ANOVA for the model no. 5, the best one for the assemblage scale:
anova(gls.c.1_long, gls.c.2_long) 



###########################################################################################################################################

###########################################################################################################################################


# Revision ----------------------------------------------------------------

#Added Code to address small sample size

#1) Pooling together neanics and adults to assess changes in PC1 between subzones with mood's median test

#first, quantify effect size
#Hedges' g P. bisulcata
#P bisulcata
cohen.d(P.bisulcata$PC1[intersect(juv.bi, bis_1)], P.bisulcata$PC1[intersect(juv.bi, bis_2)],hedges.correction=TRUE)
cohen.d(P.bisulcata$PC1[intersect(juv.bi, bis_2)], P.bisulcata$PC1[intersect(juv.bi, bis_3)],hedges.correction=TRUE)
cohen.d(P.bisulcata$PC1[intersect(juv.bi, bis_3)], P.bisulcata$PC1[intersect(juv.bi, bis_4)],hedges.correction=TRUE)

#P bisulcata
cohen.d(P.bisulcata$PC1[intersect(nojuv.bi, bis_1)], P.bisulcata$PC1[intersect(nojuv.bi, bis_2)],hedges.correction=TRUE)
cohen.d(P.bisulcata$PC1[intersect(nojuv.bi, bis_2)], P.bisulcata$PC1[intersect(nojuv.bi, bis_3)],hedges.correction=TRUE)
cohen.d(P.bisulcata$PC1[intersect(nojuv.bi, bis_3)], P.bisulcata$PC1[intersect(nojuv.bi, bis_4)],hedges.correction=TRUE)


# Hedges' g C longiforma juveniles
cohen.d(C.longiforma$PC1[intersect(juv.lo, lon_2)], C.longiforma$PC1[intersect(juv.lo, lon_3)],hedges.correction=TRUE)
cohen.d(C.longiforma$PC1[intersect(juv.lo , lon_3)], C.longiforma$PC1[intersect(juv.lo , lon_4)],hedges.correction=TRUE)

# Hedges' g C longiforma no juveniles
cohen.d(C.longiforma$PC1[intersect(nojuv.lo, lon_2)], C.longiforma$PC1[intersect(nojuv.lo, lon_3)],hedges.correction=TRUE)
cohen.d(C.longiforma$PC1[intersect(nojuv.lo , lon_3)], C.longiforma$PC1[intersect(nojuv.lo , lon_4)],hedges.correction=TRUE)


##moods median

#make dataframes for life stage groups
P.bisulcata_juv <- P.bisulcata[juv.bi,]
P.bisulcata_nea <- P.bisulcata[nea.bi,]
P.bisulcata_adu <- P.bisulcata[adu.bi,]

P.bisulcata_ju.ne <- P.bisulcata[c(juv.bi, nea.bi),]
P.bisulcata_ne.ad <- P.bisulcata[c(nea.bi, adu.bi),]

## moods median test
mood.medtest(PC1 ~ as.factor(pool.fac), data  = P.bisulcata_juv, exact = FALSE)
mood.medtest(PC1 ~ as.factor(pool.fac), data  = P.bisulcata_ne.ad, exact = FALSE)

pairwiseMedianTest(PC1 ~ as.factor(pool.fac), data  = P.bisulcata_juv, exact = FALSE)
pairwiseMedianTest(PC1 ~ as.factor(pool.fac), data  = P.bisulcata_ne.ad, exact = FALSE)



#make dataframes for life stage groups
C.longiforma_juv <- C.longiforma[juv.lo,]
C.longiforma_nea <- C.longiforma[nea.lo,]
C.longiforma_adu <- C.longiforma[adu.lo,]

C.longiforma_ju.ne <- C.longiforma[c(juv.lo, nea.lo),]
C.longiforma_ne.ad <- C.longiforma[c(nea.lo, adu.lo),]

## moods median test
mood.medtest(PC1 ~ as.factor(pool.fac), data  = C.longiforma_juv, exact = FALSE)
mood.medtest(PC1 ~ as.factor(pool.fac), data  = C.longiforma_ne.ad, exact = FALSE)

pairwiseMedianTest(PC1 ~ as.factor(pool.fac), data  = C.longiforma_juv, exact = FALSE)
pairwiseMedianTest(PC1 ~ as.factor(pool.fac), data  = C.longiforma_ne.ad, exact = FALSE)




#2) Use developmental stage as an ordinal predictor variable and show that that model
#   is more parsimonious than the null model in C. longiforma

library(vegan)
library(nlme)

BelDat$ord.LS <- rep(NA,nrow(BelDat))
BelDat$ord.LS <- as.numeric(BelDat$dev.fac)

#P. bisulcata
lm.dev_PB <- lm(BelDat$PC1[bis] ~ BelDat$ord.LS[bis])
lm_PB <- lm(BelDat$PC1[bis] ~ 1)
plot(lm.dev_PB)
summary(lm.dev_PB)
summary(lm_PB)
anova(lm.dev_PB, lm_PB)
AICc(lm.dev_PB)
AICc(lm_PB)

#C. longiforma
lm.dev_CL <- lm(BelDat$PC1[long] ~ BelDat$ord.LS[long])
lm_CL <- lm(BelDat$PC1[long] ~ 1)
plot(lm.dev_CL)
summary(lm.dev_CL)
summary(lm_CL)
anova(lm.dev_CL, lm_CL)
AICc(lm.dev_CL)
AICc(lm_CL)

#Full assemblage
lm.dev_ass <- lm(BelDat$PC1 ~ BelDat$ord.LS)
lm_ass <- lm(BelDat$PC1 ~ 1)
plot(lm.dev_ass)
summary(lm.dev_ass)
summary(lm_ass)
anova(lm.dev_ass, lm_ass)
AICc(lm.dev_ass)
AICc(lm_ass)



############################################################################################################

# Varpart with new proxy data ---------------------------------------------

varpart_new <- varpart(BelDat$PC1[1:142], ~ BelDat$d18O_new[1:142], ~ BelDat$d13C_new[1:142], 
                       ~ BelDat$d11B_new[1:142], ~ as.factor(BelDat$lit.fac)[1:142], 
                       data = BelDat, scale = T)
varpart_new

windows(w=6, h=5)
plot(varpart_new, bg= c( "steelblue1"   , "steelblue"  , "steelblue3", "lightsalmon"),
     Xnames=c("d18O", "d13C", "d11B", "lithology"),
     id.size=1)

windows(h=7, w=6)
barplot(varpart_new$part$indfract$Adj.R.square[1:4],
        ylim=c(-0.02,0.02),
        names=c("d18O", "d13C", "d11B", "lithology"),
        col=c("dodgerblue2", "dodgerblue2", "dodgerblue2", "tomato2"), border= F,
        ylab="Independent effect (R^2)")
abline(h=0)

anova(rda(BelDat$PC1[1:142], ~ BelDat$d18O_new[1:142], ~ BelDat$d13C_new[1:142], ~ BelDat$d11B_new[1:142], 
          ~ as.factor(BelDat$lit.fac)[1:142], data = BelDat, scale = T))
RsquareAdj(rda(BelDat$PC1[1:142], ~ BelDat$d18O_new[1:142], ~ BelDat$d13C_new[1:142], ~ BelDat$d11B_new[1:142], 
               ~ as.factor(BelDat$lit.fac)[1:142],data = BelDat, scale = T))

anova(rda(BelDat$PC1 ~ BelDat$d18O_new,
          data = BelDat, scale = T)) # p = 0.002 **, adjR2 = 0.05018239
RsquareAdj(rda(BelDat$PC1 ~ BelDat$d18O_new,
               data = BelDat, scale = T))

anova(rda(BelDat$PC1 ~ BelDat$d13C_new,
          data = BelDat, scale = T)) # p = 0.48, adjR2 = -0.0033669
RsquareAdj(rda(BelDat$PC1 ~ BelDat$d13C_new,
               data = BelDat, scale = T))

anova(rda(BelDat$PC1 ~ as.factor(BelDat$lithology),
          data = BelDat, scale = T)) # p = 0.001 ***, adjR2 = 0.08563982
RsquareAdj(rda(BelDat$PC1 ~ as.factor(BelDat$lithology),
               data = BelDat, scale = T))

anova(rda(BelDat$PC1[1:142] ~ BelDat$d11B_new[1:142],
          data = BelDat, scale = T)) # 0.028 *, adjR2 = 0.02451042
RsquareAdj(rda(BelDat$PC1[1:142] ~ BelDat$d11B_new[1:142],
               data = BelDat, scale = T))

####################################################


# GLS models with new proxy data ------------------------------------------------


#correct for the impact of lithology
correct_mod.new <- lm(MOD.dat$PC1 ~ as.factor(MOD.dat$lithology))
MOD.dat.new$Residuals <- correct_mod.new$resid #add residuals to dataframe

#correlation between the variables
corr.new <- cor(MOD.dat.new[, c(3:5)], method = "spearman",
                use = "pairwise.complete.obs")

windows(h=7,w=7)
corrplot(corr.new, method = "number", type="upper", order="hclust", tl.col = "black")

################################################################################
# Assemblage level
lm.PC1.new <- lm(Residuals ~ d13C_new + d18O_new + d11B_new,
                 data = MOD.dat.new, na.action = "na.omit")
acf(residuals(lm.PC1.new)) #test for auto-correlation of residuals
durbinWatsonTest(lm.PC1.new, max.lag=10) 

################################################################################
gls.new.1 <- gls(Residuals ~ 1, correlation = corARMA(p=1), data = MOD.dat.new,  method="ML", na.action = "na.omit")
summary(gls.new.1)
gls.new.2 <- gls(Residuals ~ d13C_new, correlation = corARMA(p=1), data = MOD.dat.new,  method="ML", na.action = "na.omit")
summary(gls.new.2)
gls.new.3 <- gls(Residuals ~ d18O_new, correlation = corARMA(p=1), data = MOD.dat.new,  method="ML", na.action = "na.omit")
summary(gls.new.3)
gls.new.4 <- gls(Residuals ~ d11B_new, correlation = corARMA(p=1), data = MOD.dat.new,  method="ML", na.action = "na.omit")
summary(gls.new.4)
gls.new.5 <- gls(Residuals ~ d18O_new + d13C_new, correlation = corARMA(p=1), data = MOD.dat.new,  method="ML", na.action = "na.omit")
summary(gls.new.5)
gls.new.6 <- gls(Residuals ~ d13C_new + d11B_new, correlation = corARMA(p=1), data = MOD.dat.new,  method="ML", na.action = "na.omit")
summary(gls.new.6)
gls.new.7 <- gls(Residuals ~ d18O_new + d11B_new, correlation = corARMA(p=1), data = MOD.dat.new,  method="ML", na.action = "na.omit")
summary(gls.new.7)
gls.new.8 <- gls(Residuals ~ d18O_new + d13C_new + d11B_new, correlation = corARMA(p=1), data = MOD.dat.new,  method="ML", na.action = "na.omit")
summary(gls.new.8)

AICcscores.new <- c(AICc(gls.new.1), AICc(gls.new.2), AICc(gls.new.3), AICc(gls.new.4), AICc(gls.new.5), 
                    AICc(gls.new.6), AICc(gls.new.7), AICc(gls.new.8))
AICcscores.new #Display the AICc scores



######DEV
gls.new.1_DEV <- gls(Residuals ~ d18O_new + d13C_new + d11B_new + ord.LS, correlation = corARMA(p=1), data = MOD.dat.new,  method="ML", na.action = "na.omit")
summary(gls.new.1_DEV)
AICc(gls.new.1_DEV)
AICc(gls.new.1)

anova(gls.new.1, gls.new.1_DEV) 



################################################################################
#P.bisulcata#########################

#MOD.dat_bis.new <- MOD.dat.new[bis,]

#check if the residuals are autocorrelated
lm.PC1_bis.new <- lm(Residuals ~ d13C_new + d18O_new + d11B_new,
                     data = MOD.dat_bis.new, na.action = "na.omit")
acf(residuals(lm.PC1_bis.new)) #test for auto-correlation of residuals
durbinWatsonTest(lm.PC1_bis.new, max.lag=10) 

################################################################################
gls.new.1_bis <- gls(Residuals ~ 1, correlation = corARMA(p=1), data = MOD.dat_bis.new,  method="ML", na.action = "na.omit")
summary(gls.new.1_bis)
gls.new.2_bis <- gls(Residuals ~ d13C_new, correlation = corARMA(p=1), data = MOD.dat_bis.new,  method="ML", na.action = "na.omit")
summary(gls.new.2_bis)
gls.new.3_bis <- gls(Residuals ~ d18O_new, correlation = corARMA(p=1), data = MOD.dat_bis.new,  method="ML", na.action = "na.omit")
summary(gls.new.3_bis)
gls.new.4_bis <- gls(Residuals ~ d11B_new, correlation = corARMA(p=1), data = MOD.dat_bis.new,  method="ML", na.action = "na.omit")
summary(gls.new.4_bis)
gls.new.5_bis <- gls(Residuals ~ d18O_new + d13C_new, correlation = corARMA(p=1), data = MOD.dat_bis.new,  method="ML", na.action = "na.omit")
summary(gls.new.5_bis)
gls.new.6_bis <- gls(Residuals ~ d13C_new + d11B_new, correlation = corARMA(p=1), data = MOD.dat_bis.new,  method="ML", na.action = "na.omit")
summary(gls.new.6_bis)
gls.new.7_bis <- gls(Residuals ~ d18O_new + d11B_new, correlation = corARMA(p=1), data = MOD.dat_bis.new,  method="ML", na.action = "na.omit")
summary(gls.new.7_bis)
gls.new.8_bis <- gls(Residuals ~ d18O_new + d13C_new + d11B_new, correlation = corARMA(p=1), data = MOD.dat_bis.new,  method="ML", na.action = "na.omit")
summary(gls.new.8_bis)

AICcscores_bis.new <- c(AICc(gls.new.1_bis), AICc(gls.new.2_bis), AICc(gls.new.3_bis), AICc(gls.new.4_bis), AICc(gls.new.5_bis), 
                        AICc(gls.new.6_bis), AICc(gls.new.7_bis), AICc(gls.new.8_bis))
AICcscores_bis.new #Display the AICc scores

# Deriving a global p-value using ANOVA for the model no. 5, the best one for the assemblage scale:
anova(gls.new.1_bis, gls.new.4_bis) 


################################################################################
#C. longiforma
#MOD.dat_long.new <- MOD.dat.new[long,]

#check if the residuals are autocorrelated
lm.PC1_long.new <- lm(Residuals ~ d13C_new + d18O_new + d11B_new,
                      data = MOD.dat_long.new, na.action = "na.omit")
windows(h=7, w=8)
acf(residuals(lm.PC1_long.new)) #test for auto-correlation of residuals
durbinWatsonTest(lm.PC1_long.new, max.lag=10)

################################################################################
gls.new.1_long <- gls(Residuals ~ 1, correlation = corARMA(p=1), data = MOD.dat_long.new,  method="ML", na.action = "na.omit")
summary(gls.new.1_long)
gls.new.2_long <- gls(Residuals ~ d13C_new, correlation = corARMA(p=1), data = MOD.dat_long.new,  method="ML", na.action = "na.omit")
summary(gls.new.2_long)
gls.new.3_long <- gls(Residuals ~ d18O_new, correlation = corARMA(p=1), data = MOD.dat_long.new,  method="ML", na.action = "na.omit")
summary(gls.new.3_long)
gls.new.4_long <- gls(Residuals ~ d11B_new, correlation = corARMA(p=1), data = MOD.dat_long.new,  method="ML", na.action = "na.omit")
summary(gls.new.4_long)
gls.new.5_long <- gls(Residuals ~ d18O_new + d13C_new, correlation = corARMA(p=1), data = MOD.dat_long.new,  method="ML", na.action = "na.omit")
summary(gls.new.5_long)
gls.new.6_long <- gls(Residuals ~ d13C_new + d11B_new, correlation = corARMA(p=1), data = MOD.dat_long.new,  method="ML", na.action = "na.omit")
summary(gls.new.6_long)
gls.new.7_long <- gls(Residuals ~ d18O_new + d11B_new, correlation = corARMA(p=1), data = MOD.dat_long.new,  method="ML", na.action = "na.omit")
summary(gls.new.7_long)
gls.new.8_long <- gls(Residuals ~ d18O_new + d13C_new + d11B_new, correlation = corARMA(p=1), data = MOD.dat_long.new,  method="ML", na.action = "na.omit")
summary(gls.new.8_long)

AICcscores_long.new <- c(AICc(gls.new.1_long), AICc(gls.new.2_long), AICc(gls.new.3_long), AICc(gls.new.4_long), AICc(gls.new.5_long), 
                         AICc(gls.new.6_long), AICc(gls.new.7_long), AICc(gls.new.8_long))
AICcscores_long.new #Display the AICc scores

# Deriving a global p-value using ANOVA for the model no. 5, the best one for the assemblage scale:
anova(gls.new.1_long, gls.new.8_long) 


############################################################################

