############################SALMON RUN#########################################
#########################STAT 6570 FINAL#######################################
#############################C.Sousa###########################################

# read in .csv data
data<-read.csv('pinkdata.csv')
#choose only even-year salmon
data<-data[which(data$BY%%2==0),]
#reserve 1991-1996 for testing
testdat<-data[which(data$BY>1990),]
#Use 1950-1990 data for training
data<-data[which(data$BY<=1990),]
data<-as.data.frame(data)

# plot even year salmon recruits versus spawners (in log scale)
plot(log(data$R), log(data$S), main="Recruits vs. Spawners (Even Years)")
abline(a=0,b=1, lwd=2)
abline(lm(log(data$S)~log(data$R), data = data), col="red", lwd=2) #population appears to not be stable


data$rep<-log(data$R/data$S)#reproduction rate
summary(data$rep)
plot(data$BY,data$rep, col=data$Region)
abline(h=0)
with(data[data$Region=="Alaska Peninsula",], plot(BY, rep))

#plot reproduction rate by sites. (see if there is any spatial correlation here?)
plot(data$Stock,data$rep, las=2)

