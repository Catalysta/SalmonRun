############################SALMON RUN#########################################
#########################STAT 6570 FINAL#######################################
#############################C.Sousa###########################################

# Set working directory
# read in .csv data
data<-read.csv('pinkdata.csv')
#choose only even-year salmon
data<-data[which(data$BY%%2==0),]

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
surv<-rep(NA,nrow(data))
R2<-rep(NA,nrow(data))
for (s in levels(data$Stock)){
  indices<-which(data$Stock==s)
  l<-length(indices)
  for (i in indices[1:(l-1)]){
    surv[i+1]<-data$S[i+1]/data$R[i]
    R2[i+1]<-data$R[i]
  }
}
data<-cbind(data,surv,R2)
remove(surv,R2,indices,l,s,i)
summary(data$surv)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.01631 0.32283 0.50739 0.54621 0.78561 1.87149      35 
which(data$surv>1) #this observation had more returning spawners than it had recruits - probably measurement error
#[1] 538
#Remove this observation
data<-data[-538,]
#Some observations have 100% survival rate- perturb these to 
#avoid NaN's in logit transform
which(data$surv==1)
#[1]   4   6   8 228 230 297 298 299 326 327 328 332 333 371 372 383
#[17] 384 399 475 477
data[which(data$surv==1),"surv"]<-0.9999


#reserve 1992-1996 for testing
testdat<-data[which(data$BY>1990),]
testdat<-as.data.frame(testdat)
#Use 1950-1990 data for training
data<-data[which(data$BY<=1990),]
data<-as.data.frame(data)

#Export data for JAGS simulation (use 1972 and later)
#Collapse sites into per Region per year survival stats
#However must account for overdispersion from different sites with
#potentially different viability
S.sum<-numeric()
R2.sum<-numeric()
reg<-numeric()
year<-numeric()
k=1
for(r in levels(data$Region)){
  years<-unique(data[which(data$Region==r),]$BY)
  for (y in years){
    S.sum[k]<-sum(data[which(data$Region==r & data$BY==y),"S"])
    R2.sum[k]<-sum(data[which(data$Region==r & data$BY==y),"R2"])
    reg[k]<-r
    year[k]<-y
    k<-(k+1)
  }
}
result<-as.data.frame(cbind(reg,year,S.sum,R2.sum))
result$year<-as.integer(as.character(result$year))
for (s in levels(result$reg)){
  region<- ifelse(result$reg==s,1,0)
  result<-cbind(result,region)
}
colnames(result)[3:16]<-c("S","R2","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12")
levels(result$reg)<-levels(data$Region)
write.table(result[which(result$year>1970),],"salmon.txt",sep="\t",row.names=FALSE)
remove(k,r,R2.sum,reg,region,s,S.sum,y,year,years)


##################################EDA######################################

par(mfrow=c(1,2))
# plot even year salmon spawners vs. recruits (in log scale)
plot(log(data$S), log(data$R), main="Spawners vs. Recruits (log scale, Even Years)")
abline(a=0,b=1, lwd=2)
abline(lm(log(data$R)~log(data$S), data = data), 
       col="red", lwd=2) 
legend("topleft", cex = 0.75, legend = c("Replacement Rate", "Estimated Rate"), 
       lty = c(1,1), col = c("black", "red"))


# plot even year salmon recruits vs. spawners (in log scale)
plot(log(data$R2), log(data$S), main="Recruits vs. Spawners (log scale, Even Years)")
abline(a=0,b=1, lwd=2)
abline(lm(log(data$S)~log(data$R2), data = data), 
       col="red", lwd=2) 
legend("topleft", cex = 0.75, legend = c("100% Survival", "Estimated Survival"), 
       lty = c(1,1), col = c("black", "red"))

#plot reproduction rate by sites. (see if there is any spatial correlation here?)
#reproduction rates appear to be rather stable
par(mfrow=c(1,1),mar=c(10,3,3,3))
#plot(data$Stock,data$rep, las=2, main = "Reproduction Rate by Stock")

#plot survival rate by sites
#plot(data[-458,"Stock"],data[-458,"surv"], las=2, main = "Survival Rate by Stock")

#plot survival rate by region
plot(data[,"Region"],data[,"surv"],las=2, main = "Survival Rate by Region")

par(mfrow=c(1,1),mar=c(5,3,3,3))
#plot survival rate by year 
#suggests possibly to only look at data since 1970 (change in fisheries mgt?)
plot(as.factor(data[,"BY"]), data[,"surv"], las=2, main = "Survival Rate by Year")

par(mfrow=c(1,1),mar=c(10,3,3,3))
#plot survival rate by region since 1972 (still leaves us with 313 data points)
indices<-which(data$BY>1970)
plot(data[indices,"Region"],data[indices,"surv"],las=2, main = "Survival Rate by Region (>1970)")

#plot logit survival rate by region since 1972 (suspect heteroscedasticity)
#could leave out Prince William Sound
plot(data[indices,"Region"],log(data[indices,"surv"]/(1-data[indices,"surv"])),las=2, main = "Logit(Survival) by Region")

par(mar=c(5,4,4,2)+0.1)


############################JAGS####################################
#Load rjags Library
library(rjags)

#set up expit function to facilitate interpretation
expit<-function(x){
  return(exp(x)/(1+exp(x)))
}

# Read in the Salmon data
salmon.data <- read.table("salmon.txt",header=T,sep='\t')
# Remove NA (initial years for which survival rates not available)
salmon.data <- salmon.data[-which(is.na(salmon.data$R2)),]
# Multiply S and R2 by 100 to give discrete counts for binomial model
salmon.data$S<-as.integer(100*salmon.data$S)
salmon.data$R2<-as.integer(100*salmon.data$R2)
# Remove records with 100% survival rates
salmon.data<-salmon.data[-which(salmon.data$S==salmon.data$R2),]

attach(salmon.data)

N = dim(salmon.data)[1]

######################################
#Model 1 -- logistic regression
######################################

# Create a data list
dataList = list("S"=S,
                "R2" = R2,
                "X1" = X1,
                "X2" = X2,
                "X3" = X3,
                "X4" = X4,
                "X5" = X5,
                "X6" = X6,
                "X7" = X7,
                "X8" = X8,
                "X9" = X9,
                "X10" = X10,
                "X11" = X11,
                "X12" = X12,
                "N" = N)

# List of parameters to be monitored  
parameters = c(
  "beta01",
  "beta02",
  "beta03",
  "beta04",
  "beta05",
  "beta06",
  "beta07",
  "beta08",
  "beta09",
  "beta10",
  "beta11",
  "beta12")

# Set initial values
initsValues = list(
  "beta01" = 0,
  "beta02" = 0,
  "beta03" = 0,
  "beta04" = 0,
  "beta05" = 0,
  "beta06" = 0,
  "beta07" = 0,
  "beta08" = 0,
  "beta09" = 0,
  "beta10" = 0,
  "beta11" = 0,
  "beta12" = 0) 

# JAGS Set-up
adaptSteps = 5000              #number of steps to "tune" the samplers
burnInSteps = 5000             #number of steps to "burn-in" the samplers
nChains = 2                    #number of chains to run
numSavedSteps = 10000          #total number of steps in chains to save
thinSteps = 1                 	#number of steps to "thin" (1=keep every step)
nIter = ceiling((numSavedSteps*thinSteps )/nChains) 	#steps per chain

# Create, initialize, and adapt the model
jagsModel1 = jags.model("model2.txt", 
                        data=dataList, 
                        inits=initsValues, 
                        n.chains=nChains, 
                        n.adapt=adaptSteps)

# Burn-in the algorithm
cat( "Burning in the MCMC chain...\n")
update(jagsModel1, n.iter=burnInSteps)

# Run MCMC algorithm
cat("Sampling final MCMC chain...\n" )
codaSamples = coda.samples(jagsModel1, 
                           variable.names=parameters, 
                           n.iter=nIter, 
                           thin=thinSteps)

# Make trace plots and density plots
par(ask=F)
plot(codaSamples)

# Calculate numerical summaries for the posterior samples
summary(codaSamples)

# Retrieve posterior samples
mcmcChain = as.matrix(codaSamples)

#Look at boxplots of the posteriors samples of the betas
par(mfrow=c(1,1), ask = F, mar=c(10,3,3,3))
boxplot(as.data.frame(cbind(mcmcChain[, "beta01"], 
                            mcmcChain[ ,"beta02"], 
                            mcmcChain[ ,"beta03"],
                            mcmcChain[ ,"beta04"],
                            mcmcChain[ ,"beta05"],
                            mcmcChain[ ,"beta06"],
                            mcmcChain[ ,"beta07"],
                            mcmcChain[ ,"beta08"],
                            mcmcChain[ ,"beta09"],
                            mcmcChain[ ,"beta10"],
                            mcmcChain[ ,"beta11"],
                            mcmcChain[ ,"beta12"])), 
        names=levels(salmon.data$reg), las=2,
        main = "Posterior Distributions of Logit(Survival)")
abline(h=0)

means<-numeric(12)
for (i in 1:12){
  means[i]<-mean(as.data.frame(mcmcChain)[,i])
}
means
expit(means)

dic1 = dic.samples(jagsModel1, nIter)
dic1

##################################################################
#Model 2 -- overdispersed logistic regression (DOES NOT CONVERGE) :(
##################################################################

# Create a data list
dataList = list("S"=S,
                "R2" = R2,
                "X1" = X1,
                "X2" = X2,
                "X3" = X3,
                "X4" = X4,
                "X5" = X5,
                "X6" = X6,
                "X7" = X7,
                "X8" = X8,
                "X9" = X9,
                "X10" = X10,
                "X11" = X11,
                "X12" = X12,
                "N" = N)

# List of parameters to be monitored  
parameters = c("beta01",
               "beta02",
               "beta03",
               "beta04",
               "beta05",
               "beta06",
               "beta07",
               "beta08",
               "beta09",
               "beta10",
               "beta11",
               "beta12", 
               "alpha",
               "sigma")

# Set initial values
initsValues = list("beta01" = 0,
                   "beta02" = 0,
                   "beta03" = 0,
                   "beta04" = 0,
                   "beta05" = 0,
                   "beta06" = 0,
                   "beta07" = 0,
                   "beta08" = 0,
                   "beta09" = 0,
                   "beta10" = 0,
                   "beta11" = 0,
                   "beta12" = 0,
                   "alpha" = rep(0,N),
                   "sigma" = 1) 

# JAGS Set-up
adaptSteps = 5000              #number of steps to "tune" the samplers
burnInSteps = 5000             #number of steps to "burn-in" the samplers
nChains = 2                    #number of chains to run
numSavedSteps = 10000           #total number of steps in chains to save
thinSteps = 1                  #number of steps to "thin" (1=keep every step)
nIter = ceiling((numSavedSteps*thinSteps )/nChains) 	#steps per chain

# Create, initialize, and adapt the model
jagsModel2 = jags.model("model2_ovr.txt", 
                        data=dataList, 
                        inits=initsValues, 
                        n.chains=nChains, 
                        n.adapt=adaptSteps)

# Burn-in the algorithm
cat( "Burning in the MCMC chain...\n")
update(jagsModel2, n.iter=burnInSteps)

# Run MCMC algorithm
cat("Sampling final MCMC chain...\n" )
codaSamples = coda.samples(jagsModel2, variable.names=parameters, n.iter=nIter, thin=thinSteps)

# Make trace plots and density plots
par(ask=F)
plot(codaSamples)

# Calculate numerical summaries for the posterior samples
summary(codaSamples)

# Retrieve posterior samples
mcmcChain_ovr = as.matrix(codaSamples)

#Look at boxplots of the posteriors samples of the betas
#Look at boxplots of the posteriors samples of the betas
par(mfrow=c(1,1), ask = F, mar=c(10,3,3,3))
boxplot(as.data.frame(cbind(mcmcChain_ovr[, "beta01"], 
                            mcmcChain_ovr[ ,"beta02"], 
                            mcmcChain_ovr[ ,"beta03"],
                            mcmcChain_ovr[ ,"beta04"],
                            mcmcChain_ovr[ ,"beta05"],
                            mcmcChain_ovr[ ,"beta06"],
                            mcmcChain_ovr[ ,"beta07"],
                            mcmcChain_ovr[ ,"beta08"],
                            mcmcChain_ovr[ ,"beta09"],
                            mcmcChain_ovr[ ,"beta10"],
                            mcmcChain_ovr[ ,"beta11"],
                            mcmcChain_ovr[ ,"beta12"])), 
        names=levels(data$Region), las=2,
        main = "Posterior Distributions of Logit(Survival)")
abline(h=0)

alpha.samples = matrix(NA, nIter*nChains, N)
for(k in 1:N) alpha.samples[,k]  = mcmcChain_ovr[,paste("alpha[",k,"]", sep="")]

par(mfrow=c(1,1), ask = F)   
boxplot(as.data.frame(alpha.samples),names=as.character(1:N),"Posterior of the Plate Random Effects")
abline(h=0)

###############################################
# DIC for model comparison
###############################################

dic1 = dic.samples(jagsModel1, nIter)
dic1


dic2 = dic.samples(jagsModel2, nIter)
dic2

diffdic(dic1, dic2)

###############################PREDICTION#############################

# We wish to predict the number of salmon that will return
# out of N = 100 salmon from Region 1.

salmon.observed = read.table("salmon.txt",header=T)
# Remove NA (initial years for which survival rates not available)
salmon.observed <- salmon.observed[-which(is.na(salmon.observed$R2)),]
# Multiply S and R2 by 100 to give discrete counts for binomial model
salmon.observed$S<-as.integer(100*salmon.observed$S)
salmon.observed$R2<-as.integer(100*salmon.observed$R2)
# Remove records with 100% survival rates
salmon.observed<-salmon.observed[-which(salmon.observed$S==salmon.observed$R2),3:16]


salmon.unobserved = c(NA,100,1,0,0,0,0,0,0,0,0,0,0,0)

salmon = rbind(salmon.observed,salmon.unobserved)
attach(salmon)  # each column now becomes its own variable

N = dim(salmon)[1]

######################################
#Model 1 -- logistic regression
######################################

# Create a data list
dataList = list("S"=S,
                "R2" = R2,
                "X1" = X1,
                "X2" = X2,
                "X3" = X3,
                "X4" = X4,
                "X5" = X5,
                "X6" = X6,
                "X7" = X7,
                "X8" = X8,
                "X9" = X9,
                "X10" = X10,
                "X11" = X11,
                "X12" = X12,
                "N" = N)

# List of parameters to be monitored  
parameters = c("beta01",
               "beta02",
               "beta03",
               "beta04",
               "beta05",
               "beta06",
               "beta07",
               "beta08",
               "beta09",
               "beta10",
               "beta11",
               "beta12", 
               "survpred")

# Set initial values
initsValues = list(
  "beta01" = 0,
  "beta02" = 0,
  "beta03" = 0,
  "beta04" = 0,
  "beta05" = 0,
  "beta06" = 0,
  "beta07" = 0,
  "beta08" = 0,
  "beta09" = 0,
  "beta10" = 0,
  "beta11" = 0,
  "beta12" = 0) 

# JAGS Set-up
adaptSteps = 5000              #number of steps to "tune" the samplers
burnInSteps = 5000             #number of steps to "burn-in" the samplers
nChains = 2                    #number of chains to run
numSavedSteps = 10000          #total number of steps in chains to save
thinSteps = 1                 	#number of steps to "thin" (1=keep every step)
nIter = ceiling((numSavedSteps*thinSteps )/nChains) 	#steps per chain

# Create, initialize, and adapt the model
# note this model file has been modified to make predictions 
jagsModel3 = jags.model("model2pred.txt", 
                        data=dataList, 
                        inits=initsValues, 
                        n.chains=nChains, 
                        n.adapt=adaptSteps)

# Burn-in the algorithm
cat( "Burning in the MCMC chain...\n")
update(jagsModel3, n.iter=burnInSteps)

# Run MCMC algorithm
cat("Sampling final MCMC chain...\n" )
codaSamples = coda.samples(jagsModel3, 
                           variable.names=parameters, 
                           n.iter=nIter, 
                           thin=thinSteps)

# Make trace plots and density plots
par(ask=F)
plot(codaSamples)

# Calculate numerical summaries for the posterior samples
summary(codaSamples)

# Retrieve posterior samples
mcmcChain = as.matrix(codaSamples)

##################################CHECK ON TESTING SET ##############################

S.sum<-numeric()
R2.sum<-numeric()
reg<-numeric()
year<-numeric()
k=1
for(r in unique(levels(testdat$Region))){
  years<-unique(testdat[which(testdat$Region==r),]$BY)
  for (y in years){
    S.sum[k]<-sum(testdat[which(testdat$Region==r & testdat$BY==y),"S"])
    R2.sum[k]<-sum(testdat[which(testdat$Region==r & testdat$BY==y),"R2"])
    reg[k]<-r
    year[k]<-y
    k<-(k+1)
  }
}
testresult<-as.data.frame(cbind(reg,year,S.sum,R2.sum))
testresult$year<-as.integer(as.character(testresult$year))
testresult$S.sum<-as.numeric(as.character(testresult$S.sum))
testresult$R2.sum<-as.numeric(as.character(testresult$R2.sum))
testresult$surv<-testresult$S.sum/testresult$R2.sum
levels(testresult$reg)<-levels(data$Region)

expmeans<-c(0.3219796,0.4524453,0.3157405,0.4622027,0.5962606,0.4840185,
          0.6164405,0.3325325,0.9012106,0.2088575,0.4819286,0.8395565)
par(mfrow=c(4,3))
i=1
for(r in unique(levels(testresult$reg))){
  plot(testresult[which(testresult$reg==r),"year"],
       testresult[which(testresult$reg==r),"surv"],
       main = paste(r), xlab = "Year",
       ylab="Survival Rate", ylim=c(0,1))
  abline(h=expmeans[i], col = "red")
  i<-(i+1)
}


