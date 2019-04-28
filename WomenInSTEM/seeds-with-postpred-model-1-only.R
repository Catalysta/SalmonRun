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



