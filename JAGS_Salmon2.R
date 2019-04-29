############################################
## Bayesian Logistic Regression Example
## Seeds data
############################################

#Set working directory

#Load rjags Library
library(rjags)

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
#Priors -- implications of standard vague priors
######################################

beta = rnorm(500, 0, 1.5)							# try prior std. dev. of 2 instead
hist(exp(beta)/(1+exp(beta)))

y = rnorm(500, rep(0, 500), rep(1, 500))	# try prior std. dev. of 3 instead
hist(exp(y)/(1+exp(y)))

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
adaptSteps = 10000              #number of steps to "tune" the samplers
burnInSteps = 10000             #number of steps to "burn-in" the samplers
nChains = 2                    #number of chains to run
numSavedSteps = 20000          #total number of steps in chains to save
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
		names=levels(data$Region), las=2,
		main = "Posterior Distributions of Logit(Survival)")
abline(h=0)



###############################################
#Model 2 -- overdispersed logistic regression
###############################################

# Create a data list
dataList = list("S"=S,
                "R2" = R2,
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
parameters = c("beta02",
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
               "tau")

# Set initial values
initsValues = list("beta02" = 0,
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
                   "tau" = 1) 

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

############################################
## Bayesian Logistic Regression Example
## Seeds data - prediction
############################################

# We wish to predict the number of seeds that germinate
# out of N = 100 seeds of type 73 (x1=1) with cucumber
# extract added (x2 = 0). 


#Set working directory

#Load rjags Library
library(rjags)

seeds.observed = read.table("seeds.txt",header=T)
seeds.unobserved = c(NA,100,1,0)

seeds = rbind(seeds.observed,seeds.unobserved)
attach(seeds)  # each column now becomes its own variable

N = dim(seeds)[1]

######################################
#Priors -- implications of standard vague priors
######################################

beta = rnorm(500, 0, 1000)							# try prior std. dev. of 2 instead
hist(exp(beta)/(1+exp(beta)))

y = rnorm(500, rep(0, 500), runif(500, 0, 1000))	# try prior std. dev. of 3 instead
hist(exp(y)/(1+exp(y)))

######################################
#Model 1 -- logistic regression
######################################

# Create a data list
dataList = list("r"=r,
                "n" = n,
                "x1" = x1,
                "x2" = x2,
                "N" = N)

# List of parameters to be monitored  
parameters = c("beta0",
               "beta1",
               "beta2",
               "beta3",
               "germpred")

# Set initial values
initsValues = list("beta0" = 0,
                   "beta1" = 0,
                   "beta2" = 0,
                   "beta3" = 0) 

# JAGS Set-up
adaptSteps = 5000              #number of steps to "tune" the samplers
burnInSteps = 5000             #number of steps to "burn-in" the samplers
nChains = 2                    #number of chains to run
numSavedSteps = 10000          #total number of steps in chains to save
thinSteps = 1                 	#number of steps to "thin" (1=keep every step)
nIter = ceiling((numSavedSteps*thinSteps )/nChains) 	#steps per chain

# Create, initialize, and adapt the model
# note this model file has been modified to make predictions 
jagsModel1 = jags.model("model1-with-postpred.txt", 
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
par(mfrow=c(1,1), ask = F)
boxplot(as.data.frame(cbind(mcmcChain[, "beta0"], 
                            mcmcChain[ ,"beta1"], 
                            mcmcChain[ ,"beta2"], 
                            mcmcChain[ ,"beta3"])), 
        names=c("beta0","beta1","beta2","beta3"))
abline(h=0)





