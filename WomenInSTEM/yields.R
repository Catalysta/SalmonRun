################################################
## Yield Example
############################################### 

# Set working directory

# Load the rjags library
library(rjags)

# Read in the Yield data
Yield.data <- read.table("yields.txt",header=T,sep='\t')

N.batches <- dim(Yield.data)[2]    #Number of batches
N.samples <- dim(Yield.data)[1]    #Number of samples of each batch

# Create a list contining the data (note: use meaningful variable names)
dataList <- list(
    'y' = as.matrix(Yield.data),
    "N.batches"=N.batches,
    "N.samples"=N.samples)

# Create a list containing the parameters of interest 
parameters <- c(
	"theta", 
	"tau2.with", 
	"tau2.btw", 
	"sigma.btw", 
	"sigma.with", 
	"mu") 

# Set initial values for the sampler
initsValues <- list(
	"theta" = 1500, 
	"sigma.with" = 1, 
	"sigma.btw" = 1, 
	"mu" = rep(1500, N.batches))

# Specify options for the JAGS sampler
adaptSteps <- 5000              #stop adapting the proposal variance after 5000 steps
burnInSteps <- 5000             #number of samples to throw out as burn-in
nChains <- 2                    #number of chains to run (n>1 for Gelman-Rubin deiagnostic)
numSavedSteps <- 5000           #total number of iteration in each chain to save
thinSteps <- 1                  #number of iterations to "thin" (1=keep every step)
nIter <- ceiling((numSavedSteps*thinSteps )/nChains) 	#iterations per chain

# Create, initialize, and adapt the model
jagsModel <- jags.model("model.txt", 
				data=dataList, 
				inits=initsValues, 
				n.chains=nChains, 
				n.adapt=adaptSteps)

# Burn-in the algorithm
cat( "Burning in the MCMC chain...\n")
update(jagsModel, n.iter=burnInSteps)

# Run MCMC algorithm
cat("Sampling final MCMC chain...\n" )
codaSamples <- coda.samples(jagsModel, 
					variable.names=parameters, 
					n.iter=nIter, 
					thin=thinSteps)

# Make trace plots and density plots
par(ask=F)
plot(codaSamples)
               
# Calculate numerical summaries for the posterior samples
summary(codaSamples)

# Retrieve posterior samples
mcmcChain <- as.matrix(codaSamples)
thetaSamples <- mcmcChain[, "theta"]
tau2.withSamples <- mcmcChain[, "tau2.with"]
tau2.btwSamples <- mcmcChain[, "tau2.btw"]
sigma.withSamples <- mcmcChain[, "sigma.with"]
sigma.btwSamples <- mcmcChain[, "sigma.btw"]
muSamples <- matrix(NA, numSavedSteps, N.batches)
for(i in 1:N.batches) muSamples[,i] <- mcmcChain[,paste("mu[",i,"]", sep="")]

# Approximate the posterior distribution of the proportion 
# of the total variance due to between batch variation
prop.varSamples <- sigma.btwSamples^2/(sigma.btwSamples^2 + sigma.withSamples^2)
par(mfrow=c(1,1), ask=F)
hist(prop.varSamples, prob=T)


