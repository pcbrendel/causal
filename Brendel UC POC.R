# A Regression weighting approach to uncontrolled confounding
# by Paul Brendel

# Create expit function

expit <- function(p) {
  exp(p) / (1 + exp(p))
}

# Create data --------------------------------------------------------------------------------------

set.seed(1234)
n <- 100000

C <- rbinom(n, 1, .5)
U <- rbinom(n, 1, .5)
X  <- rbinom(n, 1, expit(-.5 + .5 * C + 1.5 * U))
Y  <- rbinom(n, 1, expit(-.5 + log(2) * X + .5 * C + 1.5 * U))

#note that P(Y=1|X=1,C=c,U=u)/P(Y=1|X=0,C=c,U=u) should equal expit(log(2))
#thus odds(Y=1|X=1,C=c,U=u)/odds(Y=1|X=0,C=c,U=u) = ORyx = exp(log(2)) = 2 

df <- data.frame(X, Y, C, U)

# Compare biased model to bias-free model ------------------------------------------------------------

Nobias <- glm(Y ~ X + C + U, family = binomial(link = "logit")) #model with no bias
exp(coef(Nobias)[2])
exp(coef(Nobias)[2] + (2 * summary(Nobias)$coef[2, 2]))
exp(coef(Nobias)[2] - (2 * summary(Nobias)$coef[2, 2]))
#ORyx = 2.02 (1.96, 2.09)

UCbias <- glm(Y ~ X + C, family = binomial(link = "logit")) #model with uncontrolled confounding
exp(coef(UCbias)[2])
exp(coef(UCbias)[2] + (2*summary(UCbias)$coef[2, 2]))
exp(coef(UCbias)[2] - (2*summary(UCbias)$coef[2, 2]))
#ORyx = 3.11 (3.03, 3.11)

# Model P(U=1|X,C,Y) and obtain regression parameters ------------------------------------------------

ModelU <- glm(U ~ X + C + Y, data=df, family=binomial(link = "logit"))
U0 <- coef(ModelU)[1]
U1 <- coef(ModelU)[2]
U2 <- coef(ModelU)[3]
U3 <- coef(ModelU)[4]

# Create bootstrap function --------------------------------------------------------------------------

fun <- function (cU, cUX, cUC, cUY) {
  set.seed(1234)
  est <- vector()
  nreps <- 10 #can vary number of bootstrap samples
  
  for(i in 1:nreps){
    bdf <- df[sample(1:nrow(df), n, replace=TRUE), ] #random samping with replacement
    
    pU <- expit(cU + cUX * bdf$X + cUC * bdf$C + cUY * bdf$Y) #model the probability of U
    
    combined <- bdf[rep(seq_len(nrow(bdf)), 2), ] #duplicate data
    combined$Usim <- rep(c(1, 0), each=n) #Usim=1 in first copy, Usim=0 in second copy
    combined$pU <- c(pU, 1 - pU) #when Usim=1, pU is prob of U=1; when Usim=0, pU is prob of U=0
    
    Final <- glm(Y ~ X + C + Usim, family = binomial(link = "logit"), weights=pU, data=combined)
    est[i] <- coef(Final)[2]
  }
  
  out <- list(exp(median(est)), exp(quantile(est, c(.025, .975))))
  return(out)
}

# Apply functions ------------------------------------------------------------------------------------

# Known, correct bias parameters   
fun(cU = U0, cUX = U1, cUC = U2, cUY = U3)
# ORyx = 2.02 (1.98, 2.06)

# Incorret bias parameters
fun(cU = ?, cUX = ?, cUC = ? cUY = ?)









