#########################
palms <- read.csv("ProjectAfricaMegafaunaOnly/results/tdwg_final.csv")
palms <- palms %>% filter_at(vars(max95FL_palms, LAT, LONG), any_vars(!is.na(.)))

## Data handling  ========================
palms$LEVEL_3_CO <- as.factor(palms$LEVEL_3_CO)
palms$accAfrica <- factor(palms$accAfrica, levels = c("1", "0"))

# Set negative values (i.e., positive body mass changes) to zero before log-transformation:
palms$mean_abs_BSchange[palms$mean_abs_BSchange < 0] <- 0
palms$mean_perc_BSchange[palms$mean_perc_BSchange < 0] <- 0

palms$mean_abs_BSchange <- palms$mean_abs_BSchange/1000

# sqrt transform
palms$sqrt_change <- sqrt(palms$mean_abs_BSchange)

# make new dataframe with transformed variables: 
df <- palms
df$FLpalms <- sqrt(df$max95FL_palms)
df$sqrt_change
df$BBM <- sqrt(exp(df$log_max95FL_palms_BBMmean)) 
df$currBS <- sqrt(df$curr_max95BodySize)
df$MAT <- (df$MAT)^2
df$TempSeas <- sqrt(df$TempSeas)
df$AnnualPrec <- sqrt(df$AnnualPrec)
df$PrecSeas <- sqrt(df$PrecSeas)
df$CH_Mean <- sqrt(df$CH_Mean)

# scale function
scaleVars <- function(x, col){
  # Scales target columns and creates new columns with specified suffix
  x[paste0(col)] <- scale(x[col], scale = T, center = T)
  return(x)
}
scale.col <- c("FLpalms", "BBM", "sqrt_change", "currBS", "MAT", 
               "TempSeas", "AnnualPrec", "PrecSeas", "CH_Mean")
df_scl <- scaleVars(x = df, col = scale.col)
df <- df_scl

df <- na.omit(df)

lm_full <- lm(FLpalms ~ 
                  sqrt_change*factor(accAfrica) + 
                  currBS *factor(accAfrica) + 
                  MAT +
                  TempSeas + 
                  AnnualPrec +
                  PrecSeas+
                  CH_Mean, 
                data = df, na.action = na.fail)
summary(lm_full)

#Spatial structure of residuals
library(ncf); library(spdep)
par(mfrow=c(2,2))

#Correlograms with latlon = FALSE
cor.OBL<-correlog(df$LONG, df$LAT, z=df$FLpalms, na.rm=T, increment=1, resamp=1, latlon = FALSE)
cor.OBL
plot(cor.OBL$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), 
     xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2, main="raw data, increment = 1, latlon=F")
abline(h=0)  

#With latlon = TRUE and increment=1
cor.OBL_latlon<-correlog(df$LONG, df$LAT, z=df$FLpalms, na.rm=T, increment=1, resamp=1, latlon = TRUE)  #uses km because latlon = TRUE
plot(cor.OBL_latlon$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), 
     xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2, main="raw data, increment = 1, latlon=T")
abline(h=0)                  

#With latlon = TRUE and increment=1447.369 (this is the maximum distance between cells as calculated below under max(dsts_palms)), resampling 999 times
cor.OBL_1000<-correlog(df$LONG, df$LAT, z=df$FLpalms, na.rm=T, increment=1447.369, resamp=999, latlon = TRUE)  #uses km because latlon = TRUE
plot(cor.OBL_1000$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), 
     xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2, main="raw data, increment = max, latlon=T")
abline(h=0) 
summary(cor.OBL_1000)
cor.OBL_1000$mean.of.class
cor.OBL_1000$n
cor.OBL_1000$correlation      
cor.OBL_1000$p    

#Correlogram for residuals
cor.res_1000<-correlog(df$LONG, df$LAT, z=residuals(lm_full), na.rm=T, increment=1447.369, resamp=999, latlon = TRUE)

#Plot both residuals and raw data
pdf("ProjectAfricaMegafaunaOnly/figures/FS_MoransI.pdf", height = 3, width=4)
plot(cor.OBL_1000$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), 
     xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2, main="Palms: Moran's I (raw, residuals)")
abline(h=0)                  
points(cor.res_1000$correlation, pch=1, cex=1.2)
lines(cor.res_1000$correlation, lwd=1.5)
dev.off()


#Implementing a spatial model

#Make coordinate list
coords_palms<-as.matrix(cbind(df$LONG,df$LAT), na.rm=T)
coords_palms <- coords_palms[-c(123, 131),]

par(mfrow=c(1,1))
plot(coords_palms)


#Minimum distance to connect to at least one neighbor
palms_knear <- knn2nb(knearneigh(coords_palms, k=1))
summary(palms_knear)
dsts_palms<-unlist(nbdists(palms_knear, coords_palms, longlat = TRUE))
summary(dsts_palms)
max(dsts_palms)

#Neighbour defined by dnearneigh()
palms_nb<-dnearneigh(coords_palms,0,max(4139.695), longlat=T) #To match Africa dataset
par(mfrow=c(3,4), mar=c(2,4,1,1))
plot(coords_palms)
plot(palms_nb, coords_palms, pch=20, lwd=2)
summary(palms_nb)

palms_nb_100<-dnearneigh(coords_palms,0,100, longlat=T)
summary(palms_nb_100)
plot(palms_nb_100, coords_palms, pch=20, lwd=2)

#Defining the spatial weights matrix


nb1_w_palms<-nb2listw(palms_nb, glist=NULL, style="W", zero.policy=T)
summary(nb1_w_palms, zero.policy=TRUE)
plot(nb1_w_palms, coords_palms)

library(spatialreg)
#Spatial autoregressive error model


sem_error_nb1_w_palms <-spatialreg::errorsarlm(lm_full, listw=nb1_w_palms, zero.policy=TRUE, na.action = na.omit) #zero.policy=FALSE when taking min distance to connect all cells
summary(sem_error_nb1_w_palms)


# Call:spatialreg::errorsarlm(formula = lm_full, listw = nb1_w_palms,     na.action = na.omit, zero.policy = TRUE)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.50978 -0.35263  0.06609  0.33329  1.59287 
# 
# Type: error 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error z value  Pr(>|z|)
# (Intercept)                     0.502184   0.355509  1.4126 0.1577794 ns
# sqrt_change                     0.184543   0.195313  0.9449 0.3447336 ns
# factor(accAfrica)0             -0.729377   0.367509 -1.9847 0.0471832 ns
# currBS                          0.301297   0.187984  1.6028 0.1089841 ns
# MAT                             0.014258   0.084798  0.1681 0.8664742 ns
# TempSeas                       -0.418863   0.095106 -4.4042 1.062e-05 ***
# AnnualPrec                      0.273698   0.077098  3.5500 0.0003852 ***
# PrecSeas                        0.146560   0.103256  1.4194 0.1557863 ns 
# CH_Mean                        -0.090946   0.087749 -1.0364 0.2999948 ns
# sqrt_change:factor(accAfrica)0  0.030071   0.219183  0.1372 0.8908753 ns
# factor(accAfrica)0:currBS      -0.082459   0.212652 -0.3878 0.6981914 ns
# 
# Lambda: 0.75943, LR test value: 14.322, p-value: 0.00015405
# Asymptotic standard error: 0.096201
# z-value: 7.8942, p-value: 2.8866e-15
# Wald statistic: 62.319, p-value: 2.8866e-15
# 
# Log likelihood: -98.1776 for error model
# ML residual variance (sigma squared): 0.29583, (sigma: 0.5439)
# Number of observations: 118 
# Number of parameters estimated: 13 
# AIC: 222.36, (AIC for lm: 234.68)

#Testing for spatial structure in the residuals of the spatial model

#Correlogram for spatial model
cor.SEM_1000<-correlog(df$LONG, df$LAT, z=residuals(sem_error_nb1_w_palms), na.rm=T, increment=1447.369, resamp=999, latlon = TRUE) #or 1583

#Plot both residuals and raw data
plot(cor.OBL_1000$correlation, type="b", pch=16, cex=1.2, lwd=1.5, ylim=c(-0.5, 1), xlab="Distance class", ylab="Moran's I", cex.lab=1.5, las=1, cex.axis=1.2)
abline(h=0)                  
points(cor.res_1000$correlation, pch=1, cex=1.2)
lines(cor.res_1000$correlation, lwd=1.5)
points(cor.SEM_1000$correlation, pch=15, cex=1.2)
lines(cor.SEM_1000$correlation, lwd=1.5)

#Calculate Moran's I value individually

# Statistics - 1 == perfect dispersion
# Statistics = 0 == perfect randomness
# Statistics + 1 == perfect clustering

palms_nb_moran1000<-dnearneigh(coords_palms,0,1447.369, longlat=T) #or 1583
nb_moran1000<-nb2listw(palms_nb_moran1000, glist=NULL, style="W", zero.policy=FALSE) #zero.policy=FALSE

moran.test(df$FLpalms, listw=nb_moran1000, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")
moran.test(residuals(lm_full), listw=nb_moran1000, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")
moran.test(residuals(sem_error_nb1_w_palms), listw=nb_moran1000, randomisation=TRUE, zero.policy=TRUE, alternative="two.sided")


#Residual maps
library(akima)
par(mfrow=c(2,3))
image(interp(df$LONG, df$LAT, df$FLpalms), col = terrain.colors(12))
image(interp(df$LONG, df$LAT, residuals(lm_full)), col = terrain.colors(12), main="OLS model residuals")
image(interp(df$LONG, df$LAT, residuals(sem_error_nb1_w_palms)), col = terrain.colors(12), main="SEM residuals")

plot(coords_palms, col=c("blue", "red")[sign(residuals(lm_full))/2+1.5], pch=19,
     cex=abs(residuals(lm_full))/max(residuals(lm_full)), xlab="geographical x-coordinates", ylab="geographical y-coordinates")

plot(coords_palms, col=c("blue", "red")[sign(residuals(sem_error_nb1_w_palms))/2+1.5], pch=19,
     cex=abs(residuals(sem_error_nb1_w_palms))/max(residuals(sem_error_nb1_w_palms)), xlab="geographical x-coordinates", ylab="geographical y-coordinates")

