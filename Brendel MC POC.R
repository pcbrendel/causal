# A Regression weighting approach to misclassification bias
# by Paul Brendel

# Create expit function

expit <- function(p) {
  exp(p) / (1 + exp(p))
}

# Create data --------------------------------------------------------------------------------------

set.seed(1234)
n <- 100000

C <- rbinom(n, 1, .5)
X  <- rbinom(n, 1, expit(-.5 + .5 * C))
Y  <- rbinom(n, 1, expit(-.5 + log(2) * X + .5 * C))
Xstar <- ifelse(X == 1 & Y == 1, rbinom(n, 1, .75), 
                (ifelse(X == 1 & Y == 0, rbinom(n, 1, .65),
                        (ifelse(X == 0 & Y == 1, rbinom(n, 1, .25), rbinom(n, 1, .35))))))

#note that P(Y=1|X=1,C=c,U=u)/P(Y=1|X=0,C=c,U=u) should equal expit(log(2))
#thus odds(Y=1|X=1,C=c,U=u)/odds(Y=1|X=0,C=c,U=u) = ORyx = exp(log(2)) = 2

df <- data.frame(X, Xstar, Y, C)

# Compare biased model to bias-free model ------------------------------------------------------------

Nobias <- glm(Y ~ X + C, family = binomial(link = "logit")) #model with no bias
exp(coef(Nobias)[2])
exp(coef(Nobias)[2] + (2 * summary(Nobias)$coef[2, 2]))
exp(coef(Nobias)[2] - (2 * summary(Nobias)$coef[2, 2]))
#ORyx = 2.04 (1.99, 2.10)

MCbias <- glm(Y ~ Xstar + C, family = binomial(link = "logit")) #model with misclassification bias
exp(coef(MCbias)[2])
exp(coef(MCbias)[2] + (2 * summary(MCbias)$coef[2, 2]))
exp(coef(MCbias)[2] - (2 * summary(MCbias)$coef[2, 2]))
#ORyx = 1.27 (1.24, 1.31)

# Model P(X=1|Xstar,C,Y) and obtain regression parameters -----------------------------------------------

ModelX <- glm(X ~ Xstar + C + Y, family=binomial(link="logit"))
X0 <- coefficients(ModelX)[1]
X1 <- coefficients(ModelX)[2]
X2 <- coefficients(ModelX)[3]
X3 <- coefficients(ModelX)[4]

# Create bootstrap function --------------------------------------------------------------------------

fun <- function (cX, cXXstar, cXC, cXY) {
  set.seed(1234)
  est <- vector()
  nreps <- 10 #can vary number of bootstrap samples
  
  for(i in 1:nreps){
    bdf <- df[sample(1:nrow(df), n, replace=TRUE), ] #random samping with replacement
    
    pX <- expit(cX + cXXstar * bdf$Xstar + cXC * bdf$C + cXY * bdf$Y) #model the probability of X
    
    combined <- bdf[rep(seq_len(nrow(bdf)), 2), ] #duplicate data
    combined$Xsim <- rep(c(1, 0), each=n) #Xsim=1 in first copy, Xsim=0 in second copy
    combined$pX <- c(pX, 1 - pX) #when Xsim=1, pX is prob of X=1; when Xsim=0, pX is prob of X=0
    
    Final <- glm(Y ~ Xsim + C, family = binomial(link = "logit"), weights=pX, data=combined)
    est[i] <- coef(Final)[2]
  }
  
  out <- list(exp(median(est)), exp(quantile(est, c(.025, .975))))
  return(out)
}

# Apply functions ------------------------------------------------------------------------------------

# Known, correct bias parameters 

fun(cX = X0, cXXstar = X1, cXC = X2, cXY = X3)
# ORyx = 2.04 (2.03, 2.05)

# Incorret bias parameters

# Double each bias parameter
fun(cX = X0*2, cXXstar = X1*2, cXC = X2*2, cXY = X3*2)
# ORyx = 2.69 (2.67, 2.71)
















