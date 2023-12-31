 Small Walkthrough for the updated drf function. The drf package on CRAN is a package for forest-based conditional distribution estimation of a possibly multivariate response. The estimated distribution is in a simple form which allows for simple and fast computation of different functionals of the conditional distributions such as, for example, conditional quantiles, conditional correlations or conditional probability statements. One can do a heterogeneity adjustment with DRF by obtaining the weighting function which describes the relevance of each training point for a given test point and which can further be used as an input to some other method.

 This poorly written update in R adds recent functionalities to drf, such as:
 - Ability to deal with missing values in covariates X
 - Uncertainty Quantification
 - Variable Importance


```R
# load functions and packages
library(drf)
library(mice)
source("drfnew_v2.R")

set.seed(10)
n<-1000
ntest<-round(n*0.1)

##Simulate Data that experiences both a mean as well as sd shift

x1 <- runif(n,-1,1)
x2 <- runif(n,-1,1)
x3 <- x1+ runif(n,-1,1)
X0 <- matrix(runif(7*n,-1,1), nrow=n, ncol=7)
Y <- as.matrix(rnorm(n,mean = 0.8*(x1 > 0), sd = 1 + 1*(x2 > 0)))
Xfull <- cbind(x1,x2, x3, X0)
colnames(Xfull)<-paste0("X", 1:10)


##Also add MAR missing values using ampute from the mice package
X<-ampute(Xfull)$amp

head(cbind(Y,X))

x<-matrix(c(0.2, 0.4, runif(8,-1,1)), nrow=1, ncol=10)
print(x)

# Fit DRF with uncertainty quantification to the data
DRF<-drfCI(X=X, Y=Y, B=50,num.trees=1000, min.node.size = 5)
DRFpred<-predictdrf(DRF, newdata=x)

## Sample from P_{Y| X=x}
Yxs<-Y[sample(1:n, size=n, replace = T, DRFpred$weights[1,])]


# Calculate quantile prediction as weighted quantiles
qx <- quantile(Yxs, probs = c(0.05,0.95))

# Calculate conditional mean estimate
mux <- mean(Yxs)

# Calculate uncertainty
alpha<-0.05
B<-length(DRFpred$weightsb)
qxb<-matrix(NaN, nrow=B, ncol=2)
muxb<-matrix(NaN, nrow=B, ncol=1)
for (b in 1:B){
Yxsb<-Y[sample(1:n, size=n, replace = T, DRFpred$weightsb[[b]][1,])]
qxb[b,] <- quantile(Yxsb, probs = c(0.05,0.95))
muxb[b] <- mean(Yxsb)
}
CI.lower.q1 <- qx[1] - qnorm(1-alpha/2)*sqrt(var(qxb[,1]))
CI.upper.q1 <- qx[1] + qnorm(1-alpha/2)*sqrt(var(qxb[,1]))

CI.lower.q2 <- qx[2] - qnorm(1-alpha/2)*sqrt(var(qxb[,2]))
CI.upper.q2 <- qx[2] + qnorm(1-alpha/2)*sqrt(var(qxb[,2]))

CI.lower.mu <- mux - qnorm(1-alpha/2)*sqrt(var(muxb))
CI.upper.mu <- mux + qnorm(1-alpha/2)*sqrt(var(muxb))



# True quantiles and conditional expectation
q1<-qnorm(0.05, mean=0.8 * (x[1] > 0), sd=(1+(x[2] > 0)))
q2<-qnorm(0.95, mean=0.8 * (x[1] > 0), sd=(1+(x[2] > 0)))
mu<-0.8 * (x[1] > 0)

# Plot
hist(Yxs, prob=T)
z<-seq(-6,7,by=0.01)
d<-dnorm(z, mean=0.8 * (x[1] > 0), sd=(1+(x[2] > 0)))
lines(z,d, col="darkred"  )
abline(v=q1,col="darkred" )
abline(v=q2, col="darkred" )
abline(v=qx[1], col="darkblue")
abline(v=qx[2], col="darkblue")
abline(v=mu, col="darkred")
abline(v=mux, col="darkblue")
abline(v=CI.lower.q1, col="darkblue", lty=2)
abline(v=CI.upper.q1, col="darkblue", lty=2)
abline(v=CI.lower.q2, col="darkblue", lty=2)
abline(v=CI.upper.q2, col="darkblue", lty=2)
abline(v=CI.lower.mu, col="darkblue", lty=2)
abline(v=CI.upper.mu, col="darkblue", lty=2)



## Variable importance for conditional Distribution Estimation


## For the conditional quantiles we use a measure that considers the whole distribution,
## i.e. the MMD based measure of DRF.
## Again missing values are no problem here
MMDVimp <- compute_drf_vimp(X=X,Y=Y)

sort(MMDVimp, decreasing = T)

```
