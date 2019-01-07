# A Regression weighting approach to uncontrolled confounding
# by Paul Brendel

# Create data --------------------------------------------------------------------------------------

set.seed(1234)
n <- 100000

C <- rbinom(n, 1, .5)
U <- rbinom(n, 1, .5)
X  <- rbinom(n, 1, plogis(-.5 + .5 * C + 1.5 * U))
Y  <- rbinom(n, 1, plogis(-.5 + log(2) * X + .5 * C + 1.5 * U))

#note that P(Y=1|X=1,C=c,U=u)/P(Y=1|X=0,C=c,U=u) should equal expit(log(2))
#thus odds(Y=1|X=1,C=c,U=u)/odds(Y=1|X=0,C=c,U=u) = ORyx = exp(log(2)) = 2 

df <- data.frame(X, Y, C, U)
rm(C, U, X, Y)

# Compare biased model to bias-free model ------------------------------------------------------------

no_bias <- glm(Y ~ X + C + U, family = binomial(link = "logit"), data = df) #model with no bias
exp(coef(no_bias)[2])
exp(coef(no_bias)[2] + summary(no_bias)$coef[2, 2] * qnorm(.025))
exp(coef(no_bias)[2] + summary(no_bias)$coef[2, 2] * qnorm(.975))
#ORyx = 2.02 (1.96, 2.09)

uc_bias <- glm(Y ~ X + C, family = binomial(link = "logit"), data = df) #model with uncontrolled confounding
exp(coef(uc_bias)[2])
exp(coef(uc_bias)[2] + summary(uc_bias)$coef[2, 2] * qnorm(.025))
exp(coef(uc_bias)[2] + summary(uc_bias)$coef[2, 2] * qnorm(.975))
#ORyx = 3.11 (3.03, 3.20)

# Model P(U=1|X,C,Y) and obtain regression parameters ------------------------------------------------

u_model <- glm(U ~ X + C + Y, family = binomial(link = "logit"), data = df)
U0 <- coef(u_model)[1]
U1 <- coef(u_model)[2]
U2 <- coef(u_model)[3]
U3 <- coef(u_model)[4]

# Create bootstrap function --------------------------------------------------------------------------

fun <- function (cU, cUX, cUC, cUY) {
  set.seed(1234)
  est <- vector()
  nreps <- 10 #can vary number of bootstrap samples
  
  for(i in 1:nreps){
    bdf <- df[sample(1:nrow(df), n, replace=TRUE), ] #random samping with replacement
    
    pU <- plogis(cU + cUX * bdf$X + cUC * bdf$C + cUY * bdf$Y) #model the probability of U
    
    combined <- bdf[rep(seq_len(nrow(bdf)), 2), ] #duplicate data
    combined$Usim <- rep(c(1, 0), each=n) #Usim=1 in first copy, Usim=0 in second copy
    combined$pU <- c(pU, 1 - pU) #when Usim=1, pU is prob of U=1; when Usim=0, pU is prob of U=0
    
    Final <- glm(Y ~ X + C + Usim, family = binomial(link = "logit"), weights=pU, data=combined)
    est[i] <- coef(Final)[2]
  }
  
  out <- list(exp(median(est)), exp(quantile(est, c(.025, .975))), hist(exp(est)))
  return(out)
}

# Apply functions ------------------------------------------------------------------------------------

# Known, correct bias parameters   
fun(cU = U0, cUX = U1, cUC = U2, cUY = U3)
# ORyx = 2.02 (1.98, 2.06)

# Incorret bias parameters
fun(cU = 2 * U0, cUX = 2 * U1, cUC = 2 * U2, cUY = 2 * U3)
# ORyx = 0.80 (0.79, 0.82)








