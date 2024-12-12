
# rm(list=ls())

library(data.table)
library(parallel)


###______________ START 
###______________ START 
###______________ START 

# load scaling fxn
# range02 <- function(x, newMax, newMin){ ((x - min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T)) * (newMax - newMin) + newMin) }

# load WOS data (created at beginning of this script)
p_all <- readRDS(file = "/home/sean/sean_stuff/r_stuff/2023/limited_transpiration_opinion_ms/data/p_all.rds")
class(p_all)
# write.csv(p_all, "/home/sean/sean_stuff/r_stuff/2023/limited_transpiration_opinion_ms/data/p_all.csv")

# pull out items of interest from list "p_all"
# plot number of publications by year
attributes(p_all)$names
head(p_all)
# add pubs column
p_all$pubs <- 1




###___________ number of publications by year ________###
###___________ number of publications by year ________###
library(doBy)
library(minpack.lm)
df2 <- summaryBy(pubs + TC ~ PY, data=p_all, FUN=sum, na.rm=TRUE, keep.names = TRUE)
# truncate to after 1970 and before 2023 (not 2023)
my <- min(df2$PY); max(df2$PY)
unique(df2$PY)
df2 <- subset(df2, PY < 2023) # remove 2023 (not one complete year of cites at time of writing)
# plot
# first plot scatter plot and fit exponential model... then overlay this on bargraph
plot(df2$pubs ~ df2$PY)
# simple exponential growth model: y = a * exp(b*x)
# create time step x axis (years past since 1970)
df2$year_step <- df2$PY - my 

# save df2
# write.csv(df2, "LT_pubs_summary.csv", row.names = FALSE)
x2 <- df2$year_step; max(x2); min(x2)
y2 <- df2$pubs; max(y2); min(y2)
plot(y2 ~ x2)
# y3 <- 0.067 * exp(0.142 * x2); y3
# points(y3 ~ x2, cex=0.2, col="orange")
fit1 <- nlsLM(y2 ~ a * exp(b * x2), start = list(a = 0.1, b = 0.15), algorithm="port")
pars1 <- as.list(coef(fit1)); pars1
with(pars1, curve(a * exp(b * x), add=TRUE,lwd=2, col="orange")) 

# create barplot
# plot(df2$count ~ df2$year_step, type='s')
barplot(df2$pubs ~ df2$PY, xlab='Year', ylab='Publications',
        axis.lty=1, mgp=c(3,1,0), ylim=c(0, 700))
box()
# overlay the fitted exponential model
par(new=TRUE)
plot(df2$pubs ~ df2$year_step, axes=FALSE, cex=0.2, col="red", cex.axis=0.5, xlab="", ylab="", 
     xlim=c(-0.5, 38.5), ylim=c(25, 675))
with(pars1, curve(a * exp(b * x), from=0, to=39.3, add=TRUE,lwd=2, col="orange")) 
# holy crudmonkeys. This is nearly perfectly exponential


##___________ plotting number of pubs by year ___________##
##___________ plotting number of pubs by year ___________##
##___________ plotting number of pubs by year ___________##

# load scaling fxn
range02 <- function(x, newMax, newMin){ ((x - min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T)) * (newMax - newMin) + newMin) }

jpeg(filename = "/home/sean/sean_stuff/r_stuff/2023/limited_transpiration_opinion_ms/figures/pubs_by_year_v1.jpeg", 
     height=5, width=6, units="in", res = 800)
par(mfrow=c(1,1),  oma=c(1, 1, 1, 1), mai=c(.7, .7, .3, .3), mgp=c(2.4,0.7,0))

# another possible approach
# plot points and make cex equal to total citations
library(viridisLite)
cols1 <- viridis(20, alpha=1, begin=0, end=1, direction=1, option="H"); cols1
# scale symbol size and color (transparency)
df2$sym_cex <- range02(df2$TC, 4, 0.1); head(df2)
df2$sym_col <- range02(df2$TC, 255, 10); head(df2) 

# plot
library(plotrix)  
plot(df2$pubs ~ df2$PY, xlab='Year', ylab='Annual publications', cex=df2$sym_cex,
     pch=21, bg=rgb(50, 0, 255, max = 255, alpha = df2$sym_col)  )
par(new=TRUE)
plot(df2$pubs ~ df2$year_step, axes=FALSE, cex=0.0, col="red", cex.axis=0.5, xlab="", ylab="") 
with(pars1, curve(a * exp(b * x), from=0, to=49.5, add=TRUE,lwd=2, col=cols1[15])) 

# add legend title
text(17, 600, "Annual citations")
# create ciricle sizes and spacing vectors
radii_sizes <- c(0.4, 0.8, 1.2, 1.6, 2.0)*0.6; radii_sizes
spacing <- c(530, 0, 0, 0, 0) 
es <- 20
x_y_scale <- 25.5 # should be about 2.3 this will change the **relative** distance between
spacing[2] <- (spacing[1] - (radii_sizes[1]*x_y_scale + radii_sizes[2]*x_y_scale)) - es
spacing[3] <- (spacing[1] - (radii_sizes[1]*x_y_scale + radii_sizes[2]*x_y_scale*2 + radii_sizes[3]*x_y_scale)) - es*2
spacing[4] <- (spacing[1] - (radii_sizes[1]*x_y_scale + radii_sizes[2]*x_y_scale*2 + radii_sizes[3]*x_y_scale*2 +
                               radii_sizes[4]*x_y_scale)) - es*3
spacing[5] <- (spacing[1] - (radii_sizes[1]*x_y_scale + radii_sizes[2]*x_y_scale*2 + radii_sizes[3]*x_y_scale*2 +
                               radii_sizes[4]*x_y_scale*2 + radii_sizes[5]*x_y_scale)) - es*4
# draw legend
draw.circle(20, spacing[1], radii_sizes[1], border="black", col=rgb(50, 0, 255, max = 255, alpha = 10))
draw.circle(20, spacing[2], radii_sizes[2], border="black", col=rgb(50, 0, 255, max = 255, alpha = 70))
draw.circle(20, spacing[3], radii_sizes[3], border="black", col=rgb(50, 0, 255, max = 255, alpha = 130))
draw.circle(20, spacing[4], radii_sizes[4], border="black", col=rgb(50, 0, 255, max = 255, alpha = 190))
draw.circle(20, spacing[5], radii_sizes[5], border="black", col=rgb(50, 0, 255, max = 255, alpha = 250))
# add text labels
text(15, spacing[1], "400", cex=0.7)
text(15, spacing[2], "800", cex=0.7)
text(15, spacing[3], "1200", cex=0.7)
text(15, spacing[4], "1600", cex=0.7)
text(15, spacing[5], "2000", cex=0.7)
# checking calibration (radii -> cex -> cites)
# points(21, spacing[1], cex=0.8, pch=21, bg=rgb(50, 0, 255, max = 255, alpha = 120))
# points(22, spacing[2], cex=1.6, pch=21, bg=rgb(50, 0, 255, max = 255, alpha = 120))
# points(23, spacing[3], cex=2.45, pch=21, bg=rgb(50, 0, 255, max = 255, alpha = 120))
# points(24, spacing[4], cex=3.3, pch=21, bg=rgb(50, 0, 255, max = 255, alpha = 120))
# points(25, spacing[5], cex=4.0, pch=21, bg=rgb(50, 0, 255, max = 255, alpha = 120))

dev.off()




###________ network plots ________________###
###________ network plots ________________###
###________ network plots ________________###

# rm(list = ls())

library(igraph)
library(doBy)
library(plotrix)  

# load scaling fxn
range02 <- function(x, newMax, newMin){ ((x - min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T)) * (newMax - newMin) + newMin) }

# get data
df <- readRDS(file = "/home/sean/sean_stuff/r_stuff/2023/limited_transpiration_opinion_ms/data/p_all.rds")

# add pub count column (for summarized dataframes where we want the number of pubs in each category)
df$pubs_count <- 1

# "research areas"
unique(df$SC)

# try to reduce "research areas" (SC) by taking only the first time in the description
woo <- c()
crap <- strsplit(df$SC,"; "); head(crap)
for(i in 1:length(crap)) {
  # i = 5
  sub <- as.data.frame(crap[i]); sub
  sub <- sub[1,1]; sub
  woo <- rbind(woo, paste(sub[1])); woo 
}
# add woo to list "df"
df$SC_reduced_1 <- woo
unique(df$SC_reduced_1) # reduced to 34 different "research areas"

# shorten category names and reduce the number of categories
df$SC_reduced_2 <- df$SC_reduced_1
# shorten category name
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Science & Technology - Other Topics", "Science & Technology", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Life Sciences & Biomedicine - Other Topics", "Life Sciences", df$SC_reduced_2)
# reduce the number of categories
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Environmental Sciences & Ecology", "Ecology & Physiology", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Food Science & Technology", "Nutrition & Dietetics", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Biochemistry & Molecular Biology", "Biotechnology & Molecular Biology", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Biotechnology & Applied Microbiology", "Biotechnology & Molecular Biology", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Biodiversity & Conservation", "Ecology & Physiology", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Physiology", "Ecology & Physiology", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Development Studies", "Ecology & Physiology", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Remote Sensing", "Engineering", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Pharmacology & Pharmacy", "Pharmacology", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Thermodynamics", "Physics & Biophysics", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "Biophysics", "Physics & Biophysics", df$SC_reduced_2)
df$SC_reduced_2 <- ifelse(df$SC_reduced_2 == "NA", "Life Sciences", df$SC_reduced_2)
# check
unique(df$SC_reduced_2)

# first, try a reduced network where each symbol represents a single research area (color-coded) and size indicates the number of publications
df2 <- summaryBy(pubs_count ~ SC_reduced_2, data=df, FUN=sum, na.rm=TRUE)
df2


##_________ star network summary, 2023 ___________#
##_________ star network summary, 2023 ___________#
##_________ star network summary, 2023 ___________#

# plot network as star network (all articles cite articles in table 1... so they should all radiate outward from a single point)
# set colors

# make base graph first, with the correct number of vertices
st <- make_star(nrow(df2))

# get colors
library(viridisLite)
cols1 <- viridis(nrow(df2), alpha=1, begin=0, end=1, direction=1, option="H"); cols1

# assign colors to vertices
vert_length <- length(V(st)); vert_length
V(st)$vert_col <- cols1

# assign sizes to vertices symbols to denote number of publications
df2$pubs_scaled <- range02(df2$pubs_count.sum, 30, 1)
# log transforme to provide more detail
df2$pubs_log <- log10(df2$pubs_count.sum)
# rescale log values
df2$pubs_log_scaled <- range02(df2$pubs_log, 30, 1)
# assign size array to graph object "st"
V(st)$vert_size <- df2$pubs_log_scaled

# make plot
jpeg(filename = "/home/sean/sean_stuff/r_stuff/2023/limited_transpiration_opinion_ms/figures/pubs_by_subject_2023_v1.jpeg", 
     height=4, width=6.5, units="in", res = 800)
# par(mfrow=c(1,1), oma=c(0.1, 0.5, 0.1, 0.5), mai=c(0.3, 2.6, 0.3, 1.0)) #, mgp=c(2.4,0.7,0))
par(mar=c(0, 0, 0, 0))

plot(st, edge.arrow.size=0.5, vertex.color=V(st)$vert_col, vertex.size=V(st)$vert_size, 
     vertex.frame.color="black", vertex.label=NA, edge.curved=0.2) #, margin=c(0, 0, 0, 0)) 

legend(x=-1.9, y=1.00, df2$SC_reduced_2, pch=21, col="#777777", pt.bg=V(st)$vert_col, 
       pt.cex=0.6, cex=0.6, bty="n", ncol=1, xpd = TRUE)

radius_scale <- 1
# draw second legend (symbol size)
draw.circle(x=1.7, y=0.8, radius=0.01 * radius_scale, border="black", col="lightgray")
draw.circle(x=1.7, y=0.6, radius=0.05 * radius_scale, border="black", col="lightgray")
draw.circle(x=1.7, y=0.27, radius=0.14 * radius_scale, border="black", col="lightgray")

# annotate second legend
text(x=1.4, y=0.8, "5", cex=0.6)
text(x=1.4, y=0.6, "25", cex=0.6)
text(x=1.4, y=0.27, "1200", cex=0.6)

# legend titles
text(x=-1.5, y=1.05, "Research discipline", cex=0.8)
text(x=1.55, y=1.05, "Total citations", cex=0.8)

dev.off()


##_________ star network summary, 2000 ___________#
##_________ star network summary, 2000 ___________#
##_________ star network summary, 2000 ___________#

# plot network as star network (all articles cite articles in table 1... so they should all radiate outward from a single point)
# set colors

# reduce dataset to include only reports published prior to 2001
df3 <- subset(df, PY < 2001)

# first, try a reduced network where each symbol represents a single research area (color-coded) and size indicates the number of publications
df3 <- summaryBy(pubs_count ~ SC_reduced_2, data=df3, FUN=sum, na.rm=TRUE)
df3


# make base graph first, with the correct number of vertices
st <- make_star(nrow(df3))

# get colors
library(viridisLite)
cols1 <- viridis(nrow(df3), alpha=1, begin=0, end=1, direction=1, option="H"); cols1

# assign colors to vertices
vert_length <- length(V(st)); vert_length
V(st)$vert_col <- cols1

# assign sizes to vertices symbols to denote number of publications
df3$pubs_scaled <- range02(df3$pubs_count.sum, 30, 1)
# log transforme to provide more detail
df3$pubs_log <- log10(df3$pubs_count.sum)
# rescale log values
df3$pubs_log_scaled <- range02(df3$pubs_log, 30, 1)
# assign size array to graph object "st"
V(st)$vert_size <- df3$pubs_log_scaled

# make plot
jpeg(filename = "/home/sean/sean_stuff/r_stuff/2023/limited_transpiration_opinion_ms/figures/pubs_by_subject_2000_v1.jpeg", 
     height=4, width=6.5, units="in", res = 800)
# par(mfrow=c(1,1), oma=c(0.1, 0.5, 0.1, 0.5), mai=c(0.3, 2.6, 0.3, 1.0)) #, mgp=c(2.4,0.7,0))
par(mar=c(0, 0, 0, 0))

plot(st, edge.arrow.size=0.5, vertex.color=V(st)$vert_col, vertex.size=V(st)$vert_size, 
     vertex.frame.color="black", vertex.label=NA, edge.curved=0.2) #, margin=c(0, 0, 0, 0)) 

legend(x=-1.9, y=1.00, df3$SC_reduced_2, pch=21, col="#777777", pt.bg=V(st)$vert_col, 
       pt.cex=0.6, cex=0.6, bty="n", ncol=1, xpd = TRUE)

radius_scale <- 1
# draw second legend (symbol size)
draw.circle(x=1.7, y=0.8, radius=0.01 * radius_scale, border="black", col="lightgray")
draw.circle(x=1.7, y=0.6, radius=0.05 * radius_scale, border="black", col="lightgray")
draw.circle(x=1.7, y=0.27, radius=0.14 * radius_scale, border="black", col="lightgray")

# annotate second legend
text(x=1.4, y=0.8, "5", cex=0.6)
text(x=1.4, y=0.6, "9", cex=0.6)
text(x=1.4, y=0.27, "60", cex=0.6)

# legend titles
text(x=-1.5, y=1.05, "Research discipline", cex=0.8)
text(x=1.55, y=1.05, "Total citations", cex=0.8)

dev.off()

###_____END

