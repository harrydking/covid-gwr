# Dissertation R Script - Harry King - v5.1
# R version 4.0.5
# RStudio version 1.3.959


# --- FUNCTIONS --- #

cat("\014") # clear the console 


rm(list = ls()) # clear the workspace


save.image("diss_workspace.RData") # save workspace image


# --- SETUP --- #

# set working directory
setwd("~/G3615 Dissertation/R WORKING DIR/")
getwd()

# load required packages
library(sf) # for spatial data
library(tidyverse) # for data unrangling
library(GWmodel) # to undertake GWR
library(tmap) # for mapping
library(car) # for stats / regression tools
library(MASS) # for robust regression / stepAIC
library(RColorBrewer) # to manage colours
library(ggplot2) # for better plotting
library(ggspatial) # spatial data framework for ggplot2
library(spdep) # for neighbour investigation and Moran's I
library(tictoc) # for measuring code running time


# --- IMPORT DATA ---#

# data has been combined previously using a join in R

# combine data with spatial framework
df = read_csv("working_set6_noWALES.csv") # import dataset

#shape.msoa <- st_read("england_msoa_2011_clipped.shp", quiet = T) %>% `st_crs<-`(27700) # full england clipped
shape.msoa <- st_read("england_msoa_2011_gen_clipped.shp", quiet = T) %>% `st_crs<-`(27700) # generalised and clipped - caused some issues

dim(shape.msoa)
class(shape.msoa)

shape.msoa.2 = st_make_valid(shape.msoa) # make shapefile valid

msoa.sf = left_join(shape.msoa.2, df, by = c("code" = "code"))
dim(msoa.sf)
class(msoa.sf)

# write out combined output
st_write(msoa.sf, "combined_output.shp", delete_layer = T)



# ---



# --- MAPPING INPUT VARIABLES --- #

# setting boundary box for 9 regions of England
ne.bb = st_bbox(msoa.sf %>% filter(region %in% c("North East")))
nw.bb = st_bbox(msoa.sf %>% filter(region %in% c("North West")))
yorks.bb = st_bbox(msoa.sf %>% filter(region %in% c("Yorkshire and The Humber")))
wmid.bb = st_bbox(msoa.sf %>% filter(region %in% c("West Midlands")))
emid.bb = st_bbox(msoa.sf %>% filter(region %in% c("East Midlands")))
eeng.bb = st_bbox(msoa.sf %>% filter(region %in% c("East of England")))
se.bb = st_bbox(msoa.sf %>% filter(region %in% c("South East")))
ldn.bb = st_bbox(msoa.sf %>% filter(region %in% c("London")))
sw.bb = st_bbox(msoa.sf %>% filter(region %in% c("South West")))

msoa.sf$variable = msoa.sf$covMort # setting variable for mapping below

# mapping the whole of England
tm_shape(msoa.sf)+
  tm_polygons("variable", title = "Legend", palette = "Purples",
              style = "kmeans", legend.hist = F)+
  tm_layout(title = "COVID-19 mortality rate",
            frame = F, legend.outside = T,
            legend.hist.width = 1,
            legend.format = list(digits = 1),
            legend.outside.position = c("left", "top"),
            legend.text.size = 1,
            legend.title.size = 1) #+
#tm_scale_bar(position = c("left", "top")) +
#tm_compass(position = c(0.15, 0.8)) #+
#tm_shape(regions) + tm_borders(col = "black", lwd = 1)


# mapping London on its own - too small to see in above map
tm_shape(msoa.sf, bbox = ldn.bb)+
  tm_polygons("variable", title = "Legend", palette = "Purples",
              style = "kmeans", legend.hist = T)+
  tm_layout(title = "London",
            frame = F, legend.outside = T,
            legend.hist.width = 1,
            legend.format = list(digits = 1),
            legend.outside.position = c("left", "top"),
            legend.text.size = 0.7,
            legend.title.size = 1) +
  tm_scale_bar(position = c("left", "top")) +
  tm_compass(position = c(0.15, 0.8)) #+
#tm_shape(regions) + tm_borders(col = "black", lwd = 1)



# ---



# --- EDA AND DATA TRANSFORMATION --- #

# SUMMARY STATISTICS

msoa.sf
names(msoa.sf)

summary(msoa.sf)

# set variable for summary below
summ.v = msoa.sf$covMort

mean(summ.v) # mean
sd(summ.v) # standard deviation
var(summ.v) # variance
IQR(summ.v) # interquartile range


# test for skewness
skewness(msoa.sf$covMort)

# working out optimal transformations for data
testmode = msoa.sf$covMort

testmode = sqrt(testmode)
testmode = log(testmode+0.0001)
testmode = (testmode)^(1/3)

skewness(testmode)

ans1 = ((0.97+0.38+0.71+1.80+2.02+2.14+1.10+2.53+1.23+0.22+0.79+0.31+1.14+1.24+2.34+0.83) / 16) # calculating av skewness before
ans2 = ((0.033 + (-0.23)+ (-0.03)+ 1.55+0.36+0.35+(-0.1)+(-0.13)+0.38+(-0.08)+0.11+(-0.23)+0.17+0.25+0.23+0.18) / 16) # av after


# VARIABLE TRANSFORMATIONS #

# transform predictor
#msoa.sf$covMort = sqrt(msoa.sf$covMort)

# demographic characteristics
msoa.sf$popnO70 = sqrt(msoa.sf$popnO70)
msoa.sf$popnO80 = sqrt(msoa.sf$popnO80)
msoa.sf$popnM = log(msoa.sf$popnM+0.0001)
msoa.sf$popnBME = log(msoa.sf$popnBME+0.0001)
msoa.sf$popnDen = log(msoa.sf$popnDen+0.0001)

# meansures of deprivation
msoa.sf$imdSco = log(msoa.sf$imdSco+0.0001)
msoa.sf$unempLT = (msoa.sf$unempLT)^(1/3)
msoa.sf$housePov = log(msoa.sf$housePov+0.0001)
msoa.sf$annInc = sqrt(msoa.sf$annInc)
msoa.sf$oldAlo = log(msoa.sf$oldAlo+0.0001)

# measures of health
msoa.sf$limitDis = sqrt(msoa.sf$limitDis)
msoa.sf$incLung = sqrt(msoa.sf$incLung)
msoa.sf$hospAlco = log(msoa.sf$hospAlco+0.0001)
msoa.sf$hospCHD = log(msoa.sf$hospCHD+0.0001)
msoa.sf$emergHosp = log(msoa.sf$emergHosp+0.0001)


# density histogram

# set variable for hist
hist.v = msoa.sf$emergHosp

hist(hist.v, prob = T, main = "Log of emergHosp", 
     xlab = "emergHosp",
     col = c("#f7729e")) # #009999 - blue (input) & #f7729e - pink (transformed)
# ylim = c(0, 0.8)) 

rug(hist.v)

lines(density(hist.v,na.rm=T),
      col='olivedrab',lwd=2)



# boxplot

# overall boxplot
boxplot(msoa.sf$covMort, main = "COVID-19 mortality per 100,000 people at MSOA level")

# location boxplot
boxplot(msoa.sf$emergHosp ~ msoa.sf$region, par(mar = c(8, 4, 1, 0.5)+ 0.001),
        las = 2, xlab = "", ylab = "Emergency hospital admissions for all causes",
        col = brewer.pal(11, "Spectral"), cex.axis=0.7)



# ---



# --- INITIAL OLS MODEL --- #

msoa.ols = st_drop_geometry(msoa.sf)

m = lm(covMort~popnO70+popnO80+popnM+popnBME+popnDen+
         imdSco+unempLT+housePov+annInc+oldAlo+
         limitDis+incLung+hospAlco+hospCHD+emergHosp,
       data = msoa.ols)

summary(m)


# --- MODEL SELECTION USING stepAIC --- #

m2 = stepAIC(m, trace = F)
names(m2)
summary(m2)

# create formula from optimised model
form1 = as.formula(covMort~popnO70+popnO80+popnBME+popnDen+
                     imdSco+unempLT+annInc+oldAlo+
                     limitDis+incLung+hospAlco+hospCHD+emergHosp)


# --- OPTIMISED OLS MODEL --- #

m3 = lm(form1, data = msoa.ols)
summary(m3)



# ---



# --- CHECKING GEOMETRY & MAPPING RESIDUALS --- # 

# checking geometries to ensure no errors when mapping outliers below

any(is.na(st_dimension(msoa.sf))) # checking for empty geometries
any(is.na(st_is_valid(msoa.sf))) # checking for corrupt geometries
any(na.omit(st_is_valid(msoa.sf)) == FALSE) # checking for invalid geometries
st_is_valid(msoa.sf, reason = TRUE) # querying the reason for invalidity (if present)


# determine studentised residuals and attach to MSOA
s.resids <- rstudent(m3)
msoa.sf$s.resids = s.resids

# map the spatial distribution of outliers across the whole of England
tm_shape(msoa.sf) +
  tm_polygons('s.resids', breaks = c(min(s.resids), -1.5, 1.5, max(s.resids)), 
              palette = c("indianred1","antiquewhite2","powderblue"),
              title = "Residuals") +
  tm_layout(legend.format= list(digits = 1), frame = F)


# mapping outliers for London
tm_shape(msoa.sf, bbox = ldn.bb) +
  tm_polygons('s.resids', breaks = c(min(s.resids), -1.5, 1.5, max(s.resids)), 
              palette = c("indianred1","antiquewhite2","powderblue"),
              title = "Residuals") +
  tm_layout(legend.format= list(digits = 1), frame = F)



# remove zero link areas
msoa.sf.2 = msoa.sf[-c(2597, 6640),]


msoa.nb = poly2nb(msoa.sf.2) 
msoa.nb

# Create a line layer showing Queen's case contiguities
msoa.net = nb2lines(msoa.nb,coords=st_geometry(st_centroid(msoa.sf.2)), as_sf = F)

# Plot the contiguity and the LSOA layer
tm_shape(msoa.sf.2) + tm_borders(col='grey') +
  tm_shape(msoa.net) + tm_lines(col='red')


# Moran's I Test
msoa.lw = nb2listw(msoa.nb)
msoa.lw

moran.test(msoa.sf.2$s.resids, msoa.lw)



2.710914 * 10^(-1) # converting Moran's I to full number


# determine range of possible I values
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/ 2) $values))
}
moran.range(msoa.lw)



# ---



# --- ROBUST LINEAR REGRESSION (RLM) --- #

m4 = rlm(form1, data = msoa.sf) # conduct RLM

summary(m4)

sum.t = coeff.comp = compareCoefs(m3,m4,print=T, pvals = T) # compare to optimised model

coef_tab <- data.frame(`Optimised OLS`=coef(m3),Robust=coef(m4))
coef_tab

write.csv((coef_tab), file = "table1.csv")


# ---



# --- GWR --- #

# convert to sp
msoa.sp = as(msoa.sf, "Spatial")
summary(msoa.sp@data)

form1 # check the formula

# create a distance matrix
DM <-gw.dist(dp.locat=coordinates(msoa.sp))

tic("gwr.bw") # to measure how long bw calculation takes - useful to know in case need to run again

# calculate bandwidth (flexible)
bw1 = bw.gwr(formula = form1, data = msoa.sp, approach="CV", kernel="bisquare",
             adaptive=T, p=2, theta=0, longlat=F, dMat = DM,
             parallel.method=F,parallel.arg=NULL)

toc() # to display elapsed time

bw1 # view the bandwidth
bw1/nrow(msoa.sp@data) 	

(bw1/nrow(msoa.sp@data))*100 # displays % of data that is highly localised

bw1.backup =  521  # in case have to run again, so don't have to recalculate


# do the GWR
gwr.m = gwr.basic(formula = form1, data = msoa.sp, bw = bw1, kernel="bisquare",
                  adaptive=T, p=2, theta=0, longlat=F, dMat = DM,F123.test=F,cv=F, W.vect=NULL,
                  parallel.method=FALSE,parallel.arg=NULL)

gwr.m # view the model

summary(gwr.m$SDF@data)

names(gwr.m$SDF@data)

summary(gwr.m$SDF@data[, 1:14])

summary(gwr.m$SDF@data[, 34:47]) # summary of t-values



# --- GWR VS OLS TABLE --- #

tab <- rbind(apply(gwr.m$SDF@data[, 1:14], 2, summary), coef(m3), coef(m4))
rownames(tab)[7] <- "Optimised OLS"
rownames(tab)[8] <- "Robust Regression"
tab <- round(tab, 3)


t(tab)

# write to table
write.csv(t(tab), file = "GWR_VS_OLS_2.csv")


# convert the sp format SDF to sf
gwr.sf = st_as_sf(gwr.m$SDF) 

gwr.sf
head(gwr.sf)

# Write out for use in GIS
st_write(gwr.sf, "gwr_output_sf.shp")




# --- COLLINEARITY TESTING --- #

# generate collinearity diagnostics

tic("collin")

collin = gwr.collin.diagno(form1,
                           data = msoa.sp,
                           bw = bw1,
                           kernel = "bisquare",
                           adaptive = T,
                           p=2, theta=0, longlat=F, dMat = DM)

toc() # 12913.8 seconds

head(collin)
summary(collin)


# histogram of VIFs
qplot(collin$SDF$popnO70_VIF,
      geom="histogram",
      main = "Histogram of popnO70 VIF", 
      xlab = "VIF", 
      fill=I("black"), 
      col=I("grey"),
      # alpha = 0.8
      #xlim=c(1,5),
) + geom_vline(xintercept = 5, linetype = "dashed", color = "blue", size = 1.5)



# histogram of VDPs

# VDP
qplot(collin$SDF$popnO70_VDP,
      main = "Histogram of popnO70 VDP",
      xlab = "VDP",
      fill=I("grey"),
      col = I("black"),
      #xlim=c(0,1),
) + geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1.5)


# local condition numbers

# Local Condition Numbers Analysis

collin[["SDF"]]@data[["local_CN"]] # local condition numbers in collin


summary(collin[["SDF"]]@data[["local_CN"]])


# this; and then do all the rest (done in terminal)
summary(collin[["SDF"]]@data[["Corr_PCemp.PCcol"]])



# ---



# --- MS-GWR --- #


# code to search for MSGWR bandwidth - not using for project, see below

# form1 # check the formula

# tic("msgwr") # start the timer

# ms.gwr <- gwr.multiscale(form1,
# data = msoa.sp, adaptive = T,
#max.iterations = 1000, kernel = "bisquare",
#bws0=rep(50, 6),
#verbose = F, predictor.centered=rep(T, 5))

# toc()


# was taking far too long, so had to abandon and use MGWR instead


# write out for use in MGWR
st_write(msoa.sf, "forMSGWR2.0.shp", delete_layer = T)


# import results from MGWR

# combine data with spatial framework
mgwr.results = read_csv("MGWR_session_2_results.csv") # import dataset

#shape.msoa.3 <- st_read("england_msoa_2011_clipped.shp", quiet = T) %>% `st_crs<-`(27700) # full england clipped
shape.msoa.3 <- st_read("england_msoa_2011_gen_clipped.shp", quiet = T) %>% `st_crs<-`(27700) # generalised and clipped - caused some issues

dim(shape.msoa.3)
class(shape.msoa.3)

shape.msoa.4 = st_make_valid(shape.msoa.3) # make shapefile valid

mgwr.res.sf = left_join(shape.msoa.4, mgwr.results, by = c("code" = "code"))
dim(mgwr.res.sf)
class(mgwr.res.sf)

# import bandwidths
mgwr.bws = read_csv("ms_gwr_bw_2.csv")

mgwr.bws.2 = data.frame(Variable = names(mgwr.bws[1]), Bandwidth = mgwr.bws[[2]])

ms.bw <- mgwr.bws.2[[2]]

ms.bw


# EXPORT TABLE

tab2 <- rbind(apply(mgwr.res.sf[, 1:14], 2, summary), coef(m3), coef(m4))
rownames(tab2)[7] <- "Optimised OLS"
rownames(tab2)[8] <- "Robust Regression"
tab2 <- round(tab2, 3)


t(tab2)

# write to table
write.csv(t(tab), file = "GWR_VS_OLS.csv")


# --- MAPPING GWR AND MS-GWR OUTPUTS --- #

# convert the sp format SDF to sf (if not done above)
gwr.sf = st_as_sf(gwr.m$SDF) 


tval1.1 = gwr.sf$Intercept_TV
signif1.1 = tval1.1 < -1.96 | tval1.1 > 1.96

tval1.2 = mgwr.res.sf$t_Intercept
signif1.2 = tval1.2 < -1.96 | tval1.2 > 1.96



# WATCH OUT FOR THE SLIGHTLY DIFFERENT VARIABLE NAMES IN SOME OF THE MS-GWR OUTPUTS
# e.g. limitDis --> limitDs


  
  tm_shape(gwr.sf) +
    tm_fill(c("Intercept"), palette = "RdBu", title = "Intercept", midpoint = 0) +
    tm_style("col_blind") +
    tm_layout(legend.position = c("left","top"), frame = F) +
    tm_facets(ncol = 1) +
    tm_shape(gwr.sf[signif1.1,]) + tm_borders("black")
  
  
  
  tm_shape(mgwr.res.sf) +
    tm_fill(c("beta_Intercept"), palette = "RdBu", midpoint = 0,
            title = paste0("BW = ", ms.bw[c(1)])) +
    tm_style("col_blind") +
    tm_layout(legend.position = c("left","top"), frame = F) +
    tm_facets(ncol = 1) +
    tm_shape(mgwr.res.sf[signif1.2,]) + tm_borders("black")
  


# end

# harry king 2021