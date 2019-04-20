# SALMON DATA - PINK
#
# **********************************************************************************
# Date:  March 12, 2011
# Author: Oksana Chkrebtii
# References: Salmon Spawning data by Dr. Peterman
#
# **********************************************************************************

# READING IN DATA

# read in .csv data
data<-read.csv('pinkdata.csv')
# add vector of observation numbers
data[['obs']]<-seq(1,dim(data)[1])
# calculate log recruits to spawners variable and add to list
data[['Y']]<-log(data$R/data$S)
# create a new variable X and give it the values in S
data[['X']]<-data$S
# create a new variable yr and give it the values in BY
data[['yr']]<-data$BY
# create a new variable asd and give it the values in AlongShore_Distance
data[['asd']]<-data$AlongShore_Distance
# convert list to data frame
df.pink<-data.frame(data)
# remove missing observations
missing<-which(is.nan(df.pink[,'Y']))
if (length(missing)>0) {df.pink<-df.pink[-missing,]}


## to consider odd year observations do this:
df.pink<-df.pink[which(df.pink[,'yr'] %% 2 != 0),]
### to consider even year observations do this:
# df.pink<-df.pink[which(df.pink[,'yr'] %% 2 == 0),]

# function to calculate Great-Circle distance between two stocks
gcdcalc<-function(lon1,lat1,lon2,lat2)
{
  r<-6378
  gcd<-r*acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos((lon1-lon2)*pi/180))
  gcd
}


# find southernmost stock location
zerogcd<-min(df.pink[,'Latitude'], na.rm=TRUE)
smost.ind<-which(df.pink[,'Latitude']==zerogcd)[1]
# measure Great Circle Distance to all other stock locations from Southernmost stock
df.pink[,'gcd']<-gcdcalc(df.pink[smost.ind,'Longitude'],df.pink[smost.ind,'Latitude'],df.pink[,'Longitude'],df.pink[,'Latitude'])

# get column names
attributes(df.pink)$names

# plot available data
par(mfrow=c(1,1))
plot(df.pink[,'asd'],df.pink[,'yr'],type='p',pch=1,xlab='Along-Shore Distance in km',ylab='Brood Year')
# pch=20 gives filled dots

# transform Y values to the interval [0,1]
ygray<-abs(df.pink[,'Y']-rep(1,length(df.pink[,'Y']))*min(df.pink[,'Y']))/abs(max(df.pink[,'Y'])-min(df.pink[,'Y']))
# plot design points in grayscale value proportional to Y value
par(mfrow=c(1,2))
plot(df.pink[,'asd'],df.pink[,'yr'],type='p',pch=20, col=gray(rep(1,length(ygray))-ygray), xlab='Along-Shore Distance in km',ylab='Brood Year')
plot(df.pink[,'gcd'],df.pink[,'yr'],type='p',pch=20, col=gray(rep(1,length(ygray))-ygray), xlab='Great-Circle Distance in km',ylab='Brood Year')
# check to make sure that none of the dots came out white
# lines(df.pink[,'asd'],df.pink[,'BY'],type='p',pch=1)

# plot design points with  marginal densities on the sides

par(mfrow=c(1,1))
def.par <- par(no.readonly = TRUE) # save default, for resetting...

d.x <- density(df.pink[,'asd'], na.rm=TRUE)
d.y <- density(df.pink[,'yr'], na.rm=TRUE)

sort(unique(df.pink[,'asd']))
line1 <- 970 + (1234-970)/2
line2 <- 2082 +(2302-2082)/2
line3 <- 3156 + (3364-3156)/2

nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
par(mar=c(3,3,1,1))
plot(df.pink[,'asd'],df.pink[,'yr'],type='p',pch=20, col=gray(rep(1,length(ygray))-ygray), xlab='Along-Shore Distance in km',ylab='Brood Year')
par(mar=c(0,3,1,1))
plot(d.x$x, d.x$y, xlim=range(df.pink[,'asd'], na.rm=TRUE), type='l', col='gray90',axes=FALSE, ann=FALSE)
polygon(d.x$x,d.x$y,col='gray90')
box()
par(mar=c(3,0,1,1))
plot(d.y$y, d.y$x, ylim=range(df.pink[,'yr'], na.rm=TRUE), xlim=range(d.y$y, na.rm=TRUE), type='l',col='gray90',axes=FALSE, ann=FALSE) 
polygon(d.y$y,d.y$x,col='gray90')
box()

par(def.par)

