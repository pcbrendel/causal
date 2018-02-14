# A Regression weighting approach to selection bias
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
S <- rbinom(n, 1, expit(-.5 + 1.5 * X + 1.5 * Y))

#note that P(Y=1|X=1,C=c,U=u)/P(Y=1|X=0,C=c,U=u) should equal expit(log(2))
#thus odds(Y=1|X=1,C=c,U=u)/odds(Y=1|X=0,C=c,U=u) = ORyx = exp(log(2)) = 2

df <- data.frame(X, Y, C, S)

# Compare biased model to bias-free model ------------------------------------------------------------

Nobias <- glm(Y ~ X + C, family = binomial(link = "logit")) #model with no bias
exp(coef(Nobias)[2])
exp(coef(Nobias)[2] + (2 * summary(Nobias)$coef[2, 2]))
exp(coef(Nobias)[2] - (2 * summary(Nobias)$coef[2, 2]))
#ORyx = 2.04 (1.99, 2.10)

set.seed(1234)
Selbias <- glm(Y ~ X + C, family = binomial(link="logit"), data = df[sample(1:nrow(df), n, replace = TRUE, S),])
exp(coef(Selbias)[2])
exp(coef(Selbias)[2] + (2 * summary(Selbias)$coef[2, 2]))
exp(coef(Selbias)[2] - (2 * summary(Selbias)$coef[2, 2]))
#ORyx = 1.32 (1.29, 1.36)

# Model P(S=1|X,Y) and obtain regression parameters ----------------------------------------------------

ModelS <- glm(S ~ X + Y, data = df, family = binomial(link = "logit"))
S0 <- coef(ModelS)[1] 
S1 <- coef(ModelS)[2] 
S2 <- coef(ModelS)[3] 

# Create bootstrap function ---------------------------------------------------------------------------

fun <- function (cS, cSX, cSY) {
  set.seed(1234)
  est <- vector()
  nreps <- 10 #can vary number of bootstrap samples
  
  for(i in 1:nreps){
    bdf <- df[sample(1:nrow(df), n, replace = TRUE, S), ] #random samping with replacement among S=1
    
    pS <- expit(cS + cSX * bdf$X + cSY * bdf$Y) #model the probability of S
    
    Final <- glm(Y ~ X + C, family = binomial(link = "logit"), weights = (1/pS), data = bdf)
    est[i] <- coef(Final)[2]
  }
  
  out <- list(exp(median(est)), exp(quantile(est, c(.025, .975))))
  return(out)
}

# Apply functions ------------------------------------------------------------------------------------

# Known, correct bias parameters 

fun(cS = S0, cSX = S1, cSY = S2)
# ORyx = 2.05 (2.01, 2.07)

# Incorret bias parameters

# Double each bias parameter
fun(cS = S0*2, cSX = S1*2, cSY = S2*2)
# ORyx = 3.91 (3.84, 3.96)





