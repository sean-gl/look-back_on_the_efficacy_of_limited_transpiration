
# This script will 
#  - fit segmented models using genotype-specific coefficients reported in Schoppach et al. 2017
#  - build summary data-sets
#  - reproduce plots S1 and S2

# Instructions
# 1) put all files into a single directory (R script, data):
#       - accumulative_t.R, climate_greeley.txt, climate_rf.txt, genotype_coef.csv()
# 2) load libraries (6 lines of code following these comments)
# 3) Run code 

# last tested by SMG on December 5, 2024, Ubuntu Linux 20.04, R 3.6.3, RStudio 2024.04.2 Build 764

# rm(list=ls())

# load libraries
library(smatr)
library(lubridate)
library(doBy)
library(viridisLite)
library(colourvalues)

df <- read.csv("genotype_coef.csv")
head(df)

# First, just quickly examining the data and results reported in Schoppach et al. 2017
# calculate leaf-specific water use *first segment* for all genos at VPD=1 kPa
df$LS_wu_seg1 <- df$slope_1 * 1 + df$y_inter; head(df)
ols1 <- sma(df$LS_wu_seg1 ~ df$yor, method="OLS"); plot(ols1); summary(ols1)
# and just for the slopes 
ols1 <- sma(df$slope_1 ~ df$yor, method="OLS"); plot(ols1); summary(ols1)
ols1 <- sma(df$slope_2 ~ df$yor, method="OLS"); plot(ols1); summary(ols1)


# calculate leaf-specific water use fxn for all values of VPD, 
# accounting for BP and slopes (1, 2)
# Although the slopes and intercepts are reported for slope 1 (first segement),
# the intercept is not reported for slope 2... so, we solve for it here

# 1) calculate the value of y (transpiration) for each BP (VPD)
# i.e., y = mx + b
df$y_bp <- df$slope_1 * df$vpd_bp + df$y_inter; head(df)   

# 2) slove for the "new" y interccept (y intercept for the second segment),
# using coefficents for the rist segment and the BP VPD 
# for example, using the genotype Wards's Prolific from Schoppach et al. 2017,
# the y value at the BP = 80.0, x = 2.06 (i.e. BP VPD), m = slope of second segment
# 80 = 1.5 * 20.6 + b; then solve for b to get the intercept ("y_int_new" for the second segment)
# gives, the new y-intercept (intercept of second segment) as 76.952
df$y_int_new <- df$y_bp - df$vpd_bp * df$slope_2; head(df)

# check plots to make sure they look okay
# create pretend df "woo"
vpd <- seq(from=0, to=3.5, length=200)
water <- seq(from=0, to=150, length=200)
woo <- data.frame(vpd, water); head(woo)
plot(woo$water ~ woo$vpd)

# plot first geno to make sure segmented fxn works
woo$wu_geno1 <- NA
row <- 1 # genotype row in df "df"

for(i in 1:nrow(woo)) {
  # i = 1
  if(woo$vpd[i] <= df$vpd_bp[row]) {
    woo$wu_geno1[i] <- df$slope_1[row] * woo$vpd[i] + df$y_inter[row]
  } else {woo$wu_geno1[i] <- df$slope_2[row] * woo$vpd[i] + df$y_int_new[row]}
}
head(woo)
plot(woo$wu_geno1 ~ woo$vpd, ylim=c(0,150), cex=0.4)
# looks good!!


###__________ calculate water use using whole-VPD rasponse function using climate data from LIRF and Rocky Ford 
## import 10 years of LIRF climate data
clim1 <- read.csv("/home/sean/sean_stuff/r_stuff/2023/limited_transpiration_opinion_ms/data/climate_data_LIRF_Rocky_Ford/10yrs_climate_GLY04")
head(clim1)
colnames(clim1) <- c('sta', 'date_time', 'temp_c', 'rh_frac', 'vp_kpa', 'sol_rad_kJ_M2_min', 'wind_m_s', 'wind_dir_azm', 'wind_dir_sd',
                     'precip_mm', 'soil_temp_5cm', 'soil_temp_15cm', 'wind_gust_m_s', 'wind_gust_time', 'wind_gust_dir_asm')
# import 10 years of Rocky Ford climate data
clim2 <- read.csv("/home/sean/sean_stuff/r_stuff/2023/limited_transpiration_opinion_ms/data/climate_data_LIRF_Rocky_Ford/10yrs_climate_RFD01")
head(clim2)
colnames(clim2) <- c('sta', 'date_time', 'temp_c', 'rh_frac', 'vp_kpa', 'sol_rad_kJ_M2_min', 'wind_m_s', 'wind_dir_azm', 'wind_dir_sd',
                     'precip_mm', 'soil_temp_5cm', 'soil_temp_15cm', 'wind_gust_m_s', 'wind_gust_time', 'wind_gust_dir_asm')


# calculate VPD for both sites
clim1$vpd_kpa <- (1-(  clim1$rh_frac)) * 0.61121 * exp((17.502 * clim1$temp_c) / (240.97 + clim1$temp_c))
clim2$vpd_kpa <- (1-(clim2$rh_frac)) * 0.61121 * exp((17.502 * clim2$temp_c) / (240.97 + clim2$temp_c))
# check data
max(clim1$vpd_kpa, na.rm=T); min(clim1$vpd_kpa, na.rm=T)  # vpd LIRF
max(clim2$vpd_kpa, na.rm=T); min(clim2$vpd_kpa, na.rm=T)  # vpd Rocky Ford
# if vpd values are less than zero, set them to zero
clim1$vpd_kpa <- ifelse(clim1$vpd_kpa < 0, 0, clim1$vpd_kpa)
clim2$vpd_kpa <- ifelse(clim2$vpd_kpa < 0, 0, clim2$vpd_kpa)
# convert to posix time
clim1$posix <- ymd_hms(clim1$date_time, tz="GMT")
clim2$posix <- ymd_hms(clim2$date_time, tz="GMT")
# combine into single climate df
clim_all <- merge(clim1[,c("vpd_kpa", "sol_rad_kJ_M2_min", "posix")], clim2[,c("vpd_kpa", "sol_rad_kJ_M2_min", "posix")], by="posix", suffixes=c('_lirf', '_rf'))
head(clim_all)
# estimate PAR from total shortwave 
crap1 <- clim_all$sol_rad_kJ_M2_min_lirf * 0.1666666 / 0.01; crap1 
crap2 <- clim_all$sol_rad_kJ_M2_min_rf * 0.1666666 / 0.01; crap2 
# now, covert W m2 to micromol m2 min
clim_all$par_lirf <- crap1 / 0.470  # https://www.apogeeinstruments.com/conversion-ppfd-to-watts/    
clim_all$par_rf <- crap2 / 0.470
# add "day", "month", and "year" columns (need to subset by day and year in loop below)
clim_all$day <- substr(clim_all$posix, 9, 10); head(clim_all$day) 
clim_all$month <- substr(clim_all$posix, 6, 7); head(clim_all$month) 
clim_all$year <- substr(clim_all$posix, 1, 4); head(clim_all$year) 

# subset to only growth months and growth hours
# get rid of year 2024 because we don't have a full season of data yet
clim_all <- subset(clim_all, year != 2024)

# subset to growth months
clim_all <- subset(clim_all, month=='05' | month=='06' | month=='07' | month=='08'); head(clim_all)
nrow(clim_all)
# reduce days of start and end months to better align with growing season
clim_all <- subset(clim_all, !(month=='05' & day<15)); nrow(clim_all)  # only keep data after the 15th of May
clim_all <- subset(clim_all, !(month=='08' & day>15)); nrow(clim_all)  # only keep data before the 20th of August

# subset to growth hours (only include data when PAR is >= 1 uM/m2/2
clim_all <- subset(clim_all, par_lirf >= 1 | par_rf >= 1); head(clim_all)

# check data for a few days
crap <- subset(clim_all, posix > as.POSIXct('2014-06-01 00:00:00', tz='GMT') &
                         posix < as.POSIXct('2014-06-06 00:00:00', tz='GMT'))
plot(crap$vpd_kpa_lirf ~ crap$posix, cex=0.2)
plot(crap$vpd_kpa_rf ~ crap$posix, cex=0.2)
plot(crap$sol_rad_kJ_M2_min_lirf ~ crap$posix, cex=0.2)
plot(crap$sol_rad_kJ_M2_min_rf ~ crap$posix, cex=0.2)
# looks good

# check all data
plot(clim_all$vpd_kpa_lirf ~ clim_all$posix, cex=0.2)
plot(clim_all$vpd_kpa_rf ~ clim_all$posix, cex=0.2)


####_________  now that we have our VPD datga from LIRF (aka Greeley) and Rocky Ford, apply these to the segmented functions in df "df"
####_________  now that we have our VPD datga from LIRF (aka Greeley) and Rocky Ford, apply these to the segmented functions in df "df"
woo2 <- c() # empty container to hold  leaf-specific water use values for all years within a single genotype
year_list <- unique(clim_all$year); year_list

# calculate leaf-specific water use for low VPD day for all genos
df_vpd <- clim_all; head(df_vpd)
colnames(df_vpd) <- c('time_seq', 'vpd_low', 'sol_rad_lirf', 'vpd_high', 'sol_rad_rf', 'par_lirf', 'par_rf', 'day', 'month', 'year')
head(df_vpd)   # note that we've changed VPD values at Greeley (LIRF) to "vpd_low" and VPD at Rocky Ford (dry site) to "VPD_high"

for(i in 1:nrow(df)) {
  # i=1
  woo1 <- c() # empty container to hold leaf-specific water use values for a single year within a single genotype
  row <- i  # genotype row in df "df"
  geno <- df$geno[row]; geno

  # subset to year
  for(j in 1:length(year_list)) {
    # j=1
    year_woo <- year_list[j]; year_woo
    year_sub <- subset(df_vpd, year==year_woo); tail(year_sub)
    
    # start hourly water use columns for low (Greeley) and high (Rocky Ford) sites... set to NA 
    year_sub$ls_wu_g1_low <- NA
    year_sub$ls_wu_g1_high <- NA

    # low VPD scenario
    for(k in 1:nrow(year_sub)) {
      # k = 1
      if(is.na(year_sub$vpd_low[k])==FALSE) {
        if(year_sub$vpd_low[k] <= df$vpd_bp[row]) {
          year_sub$ls_wu_g1_low[k] <- df$slope_1[row] * year_sub$vpd_low[k] + df$y_inter[row] # low vpd site, below BP
        } else { year_sub$ls_wu_g1_low[k] <- df$slope_2[row] * year_sub$vpd_low[k] + df$y_int_new[row] }# low vpd site, above BP
      }
    }
    # high VPD scenario
    for(k in 1:nrow(year_sub)) {
      # k = 1
      if(is.na(year_sub$vpd_high[k])==FALSE) {
        if(year_sub$vpd_high[k] <= df$vpd_bp[row]) {
          year_sub$ls_wu_g1_high[k] <- df$slope_1[row] * year_sub$vpd_high[k] + df$y_inter[row] # high vpd site, below BP
        } else { year_sub$ls_wu_g1_high[k] <- df$slope_2[row] * year_sub$vpd_high[k] + df$y_int_new[row] }# high vpd site, above BP
      }
    }
  
    tail(year_sub)
    # eliminate negative water use values (error arising from Schoppach_2017 functions (impossible y intercepts)
    year_sub$ls_wu_g1_low <- ifelse(year_sub$ls_wu_g1_low <= 0, 0, year_sub$ls_wu_g1_low)
    year_sub$ls_wu_g1_high <- ifelse(year_sub$ls_wu_g1_high <= 0, 0, year_sub$ls_wu_g1_high)
    
    # add data to woo1 
    year_sub$geno <- geno
    woo1 <- rbind(woo1, year_sub[,c('ls_wu_g1_low', 'ls_wu_g1_high', 'geno', 'time_seq', 'year')])
  }
  woo2 <- rbind(woo2, woo1); nrow(woo1); nrow(woo2)
}
  
# check dataframe woo2
unique(woo2$geno); unique(woo2$year)

# covert water use from mg/m2/s to kg/m2/hour (time step)
woo2$ls_wu_g1_low_kg_m2_h <- woo2$ls_wu_g1_low * 60 * 60 / 1000 / 1000 
woo2$ls_wu_g1_high_kg_m2_h <- woo2$ls_wu_g1_high * 60 * 60 / 1000 / 1000 
  
# create summary dataframe, giving total water use for each year within each genotype
df_sum <- summaryBy(ls_wu_g1_low_kg_m2_h + ls_wu_g1_high_kg_m2_h ~ geno + year, data=woo2, FUN=sum, na.rm=TRUE, keep.names = TRUE)
# change units... to the actual water use for the summed time period (ca 92.6 days)
colnames(df_sum) <- c('geno', 'year', 'ls_wu_g1_low_kg_m2', 'ls_wu_g1_high_kg_m2'); df_sum

# need to include year of release (yor)
crap <- merge(df[,c('geno', 'yor')], df_sum, by="geno")

# save as summarized dataframe
write.csv(crap, "schoppach_2017_summary_sean.csv", row.names = FALSE)


# now, loop through each genotype and each year and create accumulative sums
geno_list <- unique(woo2$geno); geno_list
year_list <- unique(woo2$year); year_list
foo2 <- c() # empty container to put all foo1 dfs (all years for individual genos)

for(i in 1:length(geno_list)) {
  #  i=1
  geno_sub <- geno_list[i]; geno_sub
  sub1 <- subset(woo2, geno==geno_sub); tail(sub1)
  foo1 <- c() # create a new empty container to put all years of a single geno into

  for(j in 1:length(year_list)) {
    #  j=1
    year_sub <- year_list[j]; year_sub
    sub2 <- subset(sub1, year==year_sub); head(sub2)
    
    for(k in 1:nrow(sub2)) {
      #  k=6
      if(k == 1) {
        sub2$wu_LIRF_accum_kg_m2[k] <- sub2$ls_wu_g1_low_kg_m2_h[k]
        sub2$wu_RF_accum_kg_m2[k] <- sub2$ls_wu_g1_high_kg_m2_h[k]
      } else {
        sub2$wu_LIRF_accum_kg_m2[k] <- sub2$wu_LIRF_accum_kg_m2[k-1] + sub2$ls_wu_g1_low_kg_m2_h[k]
        sub2$wu_RF_accum_kg_m2[k] <- sub2$wu_RF_accum_kg_m2[k-1] + sub2$ls_wu_g1_high_kg_m2_h[k]
      }
    }
    foo1 <- rbind(foo1, sub2); head(foo1) 
  }
  foo2 <- rbind(foo2, foo1)
}

# check data
nrow(woo2); nrow(foo2) # should be the same
unique(foo2$geno); unique(foo2$year)
head(foo2)

# need to include year of release (yor)
crap <- merge(df[,c('geno', 'yor')], foo2, by="geno")

# save as extended df
write.csv(crap, "schoppach_2017_extended_sean.csv", row.names = FALSE)



### plotting whole growing season transpiration vs year of release for LIRF
### plotting whole growing season transpiration vs year of release for LIRF
### plotting whole growing season transpiration vs year of release for LIRF

rm(list=ls())

df <- read.csv("schoppach_2017_summary_sean.csv")
head(df); nrow(df)

# create mean (across years) for each geno summary
df_sum <- summaryBy(ls_wu_g1_low_kg_m2 + ls_wu_g1_high_kg_m2 + yor ~ geno, 
                    data=df, FUN=mean, na.rm=TRUE, keep.names = TRUE)
df_sum


### plotting whole growing season transpiration vs year of release for LIRF
### plotting whole growing season transpiration vs year of release for LIRF
### plotting whole growing season transpiration vs year of release for LIRF
cols1 <- viridis(20, alpha=1, begin=0, end=1, direction=1, option="H"); cols1

jpeg(filename = "VPD_regression_schoppach.jpeg", height=4, width=10, units="in", res = 800)

par(mfrow=c(1,1),  oma=c(0, 0, 0, 0), mai=c(0, 0, 0, 0), mgp=c(2.4, 0.7, 0))

# plotting lirf (low VPD)
# placement panel 1
x1 = 0.09; x2 = 0.52
y1 = 0.18; y2 = 0.95
par(fig=c(x1=x1, x2=x2, y1=y1, y2=y2), new=FALSE)

max_y <- max(df$ls_wu_g1_low_kg_m2, df$ls_wu_g1_high_kg_m2)
min_y <- min(df$ls_wu_g1_low_kg_m2, df$ls_wu_g1_high_kg_m2)

plot(df$ls_wu_g1_low_kg_m2 ~ df$yor, ylim=c(min_y, max_y), col="lightgray",
     xlab="Genotype year of release", cex=0.2,
     ylab=expression("Leaf-specific transpiration (kg " * cm^-2 * ")"))
par(new=T)
plot(df_sum$ls_wu_g1_low_kg_m2 ~ df_sum$yor, pch=21, bg=cols1[5], type='p', ylim=c(min_y, max_y),
     yaxt='n', xaxt='n', xlab="", ylab="")

lm1 <- lm(df_sum$ls_wu_g1_low_kg_m2 ~ df_sum$yor); summary(lm1); coef1<-coef(lm1); coef1
curve(coef1[2] * x + coef1[1], add=TRUE, col=cols1[5])

# add axes
mtext("Year of release", side=1, line=2, cex=1)
mtext(expression("Leaf-specific transpiration (kg " * m^-2 * ")"), side=2, line=2, cex=1)
mtext("Greeley, Colorado (low VPD)", side=3, line=-1, cex=1)

# add r2 and p values
crap <- summary(lm1); attributes(crap)
val1 <- round(crap[[8]], 2); val1
val2 <- round(as.data.frame(crap[4])[2,4],3); val2
lab1 <- bquote(R^2~" = "~.(val1)); lab1
lab2 <- paste(expression("P = "), val2); lab2
legend('topleft', legend=c(as.expression(lab1), lab2),
       cex=0.8, bty="n", inset=c(0.2, 0.2)) 


# plotting Rocky Ford (high VPD)
# placement panel 1
x1 = 0.57; x2 = 0.98 # keep the same y values
par(fig=c(x1=x1, x2=x2, y1=y1, y2=y2), new=TRUE)

plot(df$ls_wu_g1_high_kg_m2 ~ df$yor, ylim=c(min_y, max_y), col="lightgray",
     xlab="Genotype year of release", cex=0.2,
     ylab=expression("Leaf-specific transpiration (kg " * cm^-2 * ")"),
     yaxt='n')
axis(side=2, labels=F)
par(new=T)
plot(df_sum$ls_wu_g1_high_kg_m2 ~ df_sum$yor, pch=21, bg=cols1[13], type='p', ylim=c(min_y, max_y),
     yaxt='n', xaxt='n', xlab="", ylab="")

lm1 <- lm(df_sum$ls_wu_g1_high_kg_m2 ~ df_sum$yor); summary(lm1); coef1<-coef(lm1); coef1
curve(coef1[2] * x + coef1[1], add=TRUE, col=cols1[13])

# add axes
mtext("Year of release", side=1, line=2, cex=1)
#mtext(expression("Leaf-specific transpiration (kg " * cm^-2 * ")"), side=2, line=2, cex=1)
mtext("Rocky Ford, Colorado (high VPD)", side=3, line=-1, cex=1)

# add r2 and p values
crap <- summary(lm1); attributes(crap)
val1 <- round(crap[[8]], 2); val1
val2 <- round(as.data.frame(crap[4])[2,4],3); val2
lab1 <- bquote(R^2~" = "~.(val1)); lab1
lab2 <- paste(expression("P = "), val2); lab2
legend('topleft', legend=c(as.expression(lab1), lab2),
       cex=0.8, bty="n", inset=c(0.2, 0.2)) 

dev.off()

# plot inset
# x1 = 0.1; x2 = 0.35
# y1 = 0.6; y2 = 0.92
# par(fig=c(x1=x1, x2=x2, y1=y1, y2=y2), new=TRUE)
# plot(df_sum$ls_wu_g1_low_kg_m2 ~ df_sum$yor, pch=21, bg=cols1[5], type='p', ylim=c(min_y, max_y),
#      yaxt='n', xaxt='n', xlab="", ylab="")


# stats on summary data
# plotting leaf-specific water use against year of release
# first, considering each year of climate data as independent observations (NOT VALID)
ols1 <- sma(df$ls_wu_g1_low_kg_m2 ~ df$yor); plot(ols1); summary(ols1); coef_ols1 <- coef(ols1)
ols1 <- sma(df$ls_wu_g1_high_kg_m2 ~ df$yor); plot(ols1); summary(ols1); coef_ols1 <- coef(ols1)
# now, considering the mean across all climate years as independent observations (VALID)
ols1 <- sma(df_sum$ls_wu_g1_low_kg_m2 ~ df_sum$yor); plot(ols1); summary(ols1); coef_ols1 <- coef(ols1)
ols1 <- sma(df_sum$ls_wu_g1_high_kg_m2 ~ df_sum$yor); plot(ols1); summary(ols1); coef_ols1 <- coef(ols1)




## plotting accumulative transpiration plots for each genotype across each year
## plotting accumulative transpiration plots for each genotype across each year
## plotting accumulative transpiration plots for each genotype across each year

rm(list=ls())

df <- read.csv("schoppach_2017_extended_sean.csv")

# get lists
geno_yor_lists <- summaryBy(yor~geno, data=df, FUN=mean, keep.names = TRUE); geno_yor_lists
year_list <- unique(df$year); year_list

# assign colors to genotypes
df$cols1 <- colour_values(df$yor, palette='plasma', alpha=255)

# create a plotting posix time, which will take in account all different days and months
# but have the same year... so we can plot all genos, and years on the same plot
crap1 <- substr(df$time_seq, 6, 19); head(crap1)
df$plot_posix <- paste('2014', crap1, sep="-"); head(df$plot_posix)
df$plot_posix <- ymd_hms(df$plot_posix, tz="GMT") 


#### plotting LIRF data

jpeg(filename = "2-panel_schoppach.jpeg", height=4, width=10, units="in", res = 800)

par(mfrow=c(1,1),  oma=c(0, 0, 0, 0), mai=c(0, 0, 0, 0), mgp=c(2.4, 0.7, 0))

# placement panel 1
x1 = 0.09; x2 = 0.52
y1 = 0.18; y2 = 0.95
par(fig=c(x1=x1, x2=x2, y1=y1, y2=y2), new=FALSE)

# set max and min for y axis
ymax <- max(df$wu_LIRF_accum_kg_m2, df$wu_RF_accum_kg_m2, na.rm=T)
ymin <- min(df$wu_LIRF_accum_kg_m2, df$wu_RF_accum_kg_m2, na.rm=T)

# loop through each genotype by year combination and plot data for 10 years of lirf climate data
for(i in 1:nrow(geno_yor_lists)) {
  # i=1
  sub1 <- subset(df, geno==geno_yor_lists$geno[i]); head(sub1) 
  for(j in 1:length(year_list)) {
    # j=1
    sub2 <- subset(sub1, year==year_list[j]); subset(sub2)
    if(i==1 & j==1) {
      plot(sub2$wu_LIRF_accum_kg_m2 ~ sub2$plot_posix, type='l', col=sub2$cols1, ylim=c(ymin, ymax),
           xlab="Early season growth",
           ylab=expression("Accumulative transpiration (kg " * m^-2 * ")"),
           cex.axis=0.8)
    } else {
      points(sub2$wu_LIRF_accum_kg_m2 ~ sub2$plot_posix, type='l', col=sub2$cols1)
    }
  }
}

#time_start <- as.POSIXct("2014-05-10 05:00:00", tz="GMT")
#time_end <- as.POSIXct("2014-08-21 05:00:00", tz="GMT")
#axis.POSIXct(1, at = seq(time_start, time_end, by = "day"), format = "%m/%d")

mtext("Early season growth", side=1, line=2, cex=1)
mtext(expression("Accumulative transpiration (kg " * m^-2 * ")"), side=2, line=2, cex=1)
mtext("Greeley, Colorado (low VPD)", side=3, line=-1, cex=1)

# add gradient legend
lgd_ = rep(NA, length(unique(df$yor))/2, length(unique(df$yor))); lgd_
lgd_[c(1, length(unique(df$yor))/2, length(unique(df$yor)))] = c(max(geno_yor_lists$yor), median(geno_yor_lists$yor),
                                                                 min(geno_yor_lists$yor))
legend(x = as.POSIXct("2014-05-20 05:00:00", tz="GMT"), 
       y = 430,
       legend = lgd_,
       fill = colorRampPalette(colors = c("#FADB24FF", "#9B169FFF"))(length(unique(df$yor))),
       border = NA,
       y.intersp = 0.5,
       cex = 0.7, text.font = 1)
text(x=as.POSIXct("2014-05-27 05:00:00", tz="GMT"), y=460, "Year of release", cex=0.7)


#### plotting RF data

x1 = 0.55; x2 = 0.98
par(fig=c(x1=x1, x2=x2, y1=y1, y2=y2), new=TRUE)

# set max and min for y axis
# ymax <- max(df$wu_RF_accum_kg_m2, na.rm=T)
# ymin <- min(df$wu_accum_rf, na.rm=T)

# loop through each genotype by year combination and plot data for 10 years of lirf climate data
for(i in 1:nrow(geno_yor_lists)) {
  # i=1
  sub1 <- subset(df, geno==geno_yor_lists$geno[i]); head(sub1) 
  for(j in 1:length(year_list)) {
    # j=1
    sub2 <- subset(sub1, year==year_list[j]); subset(sub2)
    if(i==1 & j==1) {
      plot(sub2$wu_RF_accum_kg_m2 ~ sub2$plot_posix, type='l', col=sub2$cols1, ylim=c(ymin, ymax),
           xlab="Early season growth",
           ylab="", yaxt='n',
           cex.axis=0.8)
    } else {
      points(sub2$wu_RF_accum_kg_m2 ~ sub2$plot_posix, type='l', col=sub2$cols1)
    }
  }
}
axis(side=2, labels=F)
mtext("Early season growth", side=1, line=2, cex=1)
mtext("Rocky Ford, Colorado (high VPD)", side=3, line=-1, cex=1)

# add gradient legend
# lgd_ = rep(NA, length(unique(df$yor))/2, length(unique(df$yor))); lgd_
# lgd_[c(1, length(unique(df$yor))/2, length(unique(df$yor)))] = c(max(geno_yor_lists$yor), median(geno_yor_lists$yor),
#                                                                 min(geno_yor_lists$yor))
#legend(x = as.POSIXct("2014-05-25 05:00:00", tz="GMT"), 
#       y = 470,
#       legend = lgd_,
#       fill = colorRampPalette(colors = c("#FADB24FF", "#9B169FFF"))(length(unique(df$yor))),
#       border = NA,
#       y.intersp = 0.5,
#       cex = 0.7, text.font = 1)
#text(x=as.POSIXct("2014-05-30 05:00:00", tz="GMT"), y=490, "Year of release", cex=0.9)

dev.off()

