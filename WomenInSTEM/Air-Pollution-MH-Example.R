################################################
## Bayesian Regression Example
## Air Quality Data
############################################### 

# Set working directory
setwd("C:/Users/chkrebtii.1/Box Sync/Teaching - STAT6570/3_Computation/R-code")

# Load MASS library
library(MASS)

# Read in data
air.data = read.table("air.txt",header=T)
N = dim(air.data)[1]

# Tranform the response variable
y = log(air.data[,1])  #use log of SO2

#############
# Some EDA
#############

# Determine explanatory variables
X = as.matrix(air.data[,-1])

# Pairs plot
pairs(cbind(y,X))

# Remove manufacturing
X.2 = X[,-2]

# Standardize the remaining variables
X.c = (X.2-matrix(rep(apply(X.2,2,mean),N),N,dim(X.2)[2],byrow=T))/matrix(rep(apply(X.2,2,sd),N),N,dim(X.2)[2],byrow=T)
cov.names = names(air.data[,c(-1,-3)])

# Add a column of 1s for the intercepts
X.int = cbind(rep(1,N),X.c)

# Determine the number of explanatory variables
P = dim(X.int)[2]

#############################################
#Classical Regression Analysis
#############################################

lm.fit = lm(y ~ X.c)
summary(lm.fit)

plot(lm.fit$fitted.values,lm.fit$residuals, ylim=c(-1.2,1.2), xlab="Fitted Values", ylab="Residuals")
abline(h=0)
abline(h=2*summary(lm.fit)$sigma,lty=2)
abline(h=-2*summary(lm.fit)$sigma,lty=2)

#############################################
#Bayesian Regression Analysis
#w/ non-informative prior
#############################################

V.beta = chol2inv(chol(t(X.int)%*%X.int))
beta.hat = V.beta%*%t(X.int)%*%y
s2 = t(y-X.int%*%beta.hat)%*%(y-X.int%*%beta.hat)/(N-P)

# Function for sampling from p(beta|sigma2,X,y)
sample.beta = function(sigma2, V.beta, beta.hat)
					{
                    		return(mvrnorm(1,beta.hat,sigma2*V.beta))
                    }
                    
# Function for sampling from p(sigma2|X,y)
sample.sigma2 = function(s2,N,P)
					{
                       return(1/rgamma(1, (N-P)/2, ((N-P)/2)*s2))
                    }
                       
# Set number of samples
K = 1000

# Initialize variables for the posterior samples
beta.post.samples = matrix(NA, K, P)
sigma2.post.samples = rep(NA, K)

for(i in 1:K)
	{
      	sigma2.curr = sample.sigma2(s2, N, P)
      	sigma2.post.samples[i] = sigma2.curr
      	beta.curr = sample.beta(sigma2.curr, V.beta, beta.hat)
      	beta.post.samples[i,] = beta.curr
     }
      
      
# Examine posterior samples
hist(sigma2.post.samples)
boxplot(as.data.frame(beta.post.samples[,-1]),names=cov.names)
abline(h=0)

# IMPORTANT:  Why can we look at boxplots of the coefficient posterior samples here?

# Calculate posterior means
apply(beta.post.samples, 2, mean)

# Look at the correlation between the betas
pairs(beta.post.samples,labels=c("Intercept",cov.names))

# Find the posterior probability that each of the parameters is greater than 0
apply(beta.post.samples > 0, 2,mean)


