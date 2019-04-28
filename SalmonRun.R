############################SALMON RUN#########################################
#########################STAT 6570 FINAL#######################################
#############################C.Sousa###########################################

# Set working directory
# read in .csv data
data<-read.csv('pinkdata.csv')
#choose only even-year salmon
data<-data[which(data$BY%%2==0),]
#reserve 1992-1996 for testing
testdat<-data[which(data$BY>1990),]
#Use 1950-1990 data for training
data<-data[which(data$BY<=1990),]
data<-as.data.frame(data)

#Sort Stock levels by along shore distance, keep only the ones in the data
data$Stock <- factor(data$Stock)
data$Stock <- reorder(data$Stock, data$AlongShore_Distance)

#Sort Region levels by along shore distance, keep only the ones in the data
data$Region <- factor(data$Region)
data$Region <- reorder(data$Region, data$AlongShore_Distance)

#Add survival rate variable
#this year's spawners are the recruits from two years ago
#also keep the R2 "recruit minus two" figures - these will be our 
#n for the binomial trials
surv<-rep(NA,478)
R2<-rep(NA,478)
for (s in levels(data$Stock)){
  indices<-which(data$Stock==s)
  l<-length(indices)
  for (i in indices[1:(l-1)]){
    surv[i+1]<-data$S[i+1]/data$R[i]
    R2[i+1]<-data$R[i]
  }
}
data<-cbind(data,surv,R2)
remove(surv,R2)
summary(data$surv)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.01631 0.32283 0.50739 0.54621 0.78561 1.87149      35 
which(data$surv>1) #this observation had more returning spawners than it had recruits - probably measurement error
#[1] 458
#Remove this observation
data<-data[-458,]
#Some observations have 100% survival rate- perturb these to 
#avoid NaN's in logit transform
which(data$surv==1)
#[1]   4   6   8 228 230 297 298 299 326 327 328 332 333 371 372 383
#[17] 384 399 475 477
data[which(data$surv==1),"surv"]<-0.9999

#Export data for JAGS simulation (use 1972 and later)
for (s in levels(data$Region)){
  region<- ifelse(data$Region==s,1,0)
  data<-cbind(data,region)
}
colnames(data)[13:24]<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12")
#write.table(data,"salmon.txt",sep="\t",row.names=FALSE)
write.table(data[which(data$BY>1970),c(9,12,13:24)],"salmon.txt",sep="\t",row.names=FALSE)

#Add reproductive rate variable
data$rep<-log(data$R/data$S)#reproduction rate
summary(data$rep)


##################################EDA######################################


# plot even year salmon recruits versus spawners (in log scale)
plot(log(data$R), log(data$S), main="Recruits vs. Spawners (Even Years)")
abline(a=0,b=1, lwd=2)
abline(lm(log(data$S)~log(data$R), data = data), 
       col="red", lwd=2) #population appears to not be stable

#plot reproduction rate by sites. (see if there is any spatial correlation here?)
#reproduction rates appear to be rather stable
par(mar=c(10,2,2,2))
plot(data$Stock,data$rep, las=2, main = "Reproduction Rate by Stock")

#plot survival rate by sites
plot(data[-458,"Stock"],data[-458,"surv"], las=2, main = "Survival Rate by Stock")
#plot survival rate by region
plot(data[-458,"Region"],data[-458,"surv"],las=2, main = "Survival Rate by Region")
#plot survival rate by year 
#suggests possibly to only look at data since 1970 (change in fisheries mgt?)
plot(as.factor(data[-458,"BY"]), data[-458,"surv"], las=2, main = "Survival Rate by Year")

#plot survival rate by region since 1972 (still leaves us with 313 data points)
indices<-which(data$BY>1970 & data$surv<=1)
plot(data[indices,"Region"],data[indices,"surv"],las=2, main = "Survival Rate by Region")

#plot logit survival rate by region since 1972 (suspect heteroscedasticity)
#could leave out Prince William Sound
plot(data[indices,"Region"],log(data[indices,"surv"]/(1-data[indices,"surv"])),las=2, main = "Logit(Survival) by Region")



par(mar=c(5,4,4,2)+0.1)
