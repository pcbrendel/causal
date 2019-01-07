# A Regression weighting approach to misclassification bias
# by Paul Brendel

# Create data --------------------------------------------------------------------------------------

set.seed(1234)
n <- 100000

C <- rbinom(n, 1, .5)
X  <- rbinom(n, 1, plogis(-.5 + .5 * C))
Y  <- rbinom(n, 1, plogis(-.5 + log(2) * X + .5 * C))
Xstar <- ifelse(X == 1 & Y == 1, rbinom(n, 1, .75), 
                (ifelse(X == 1 & Y == 0, rbinom(n, 1, .65),
                        (ifelse(X == 0 & Y == 1, rbinom(n, 1, .25), rbinom(n, 1, .35))))))

#note that P(Y=1|X=1,C=c,U=u)/P(Y=1|X=0,C=c,U=u) should equal expit(log(2))
#thus odds(Y=1|X=1,C=c,U=u)/odds(Y=1|X=0,C=c,U=u) = ORyx = exp(log(2)) = 2

df <- data.frame(X, Xstar, Y, C)
rm(C, X, Y, Xstar)

# Compare biased model to bias-free model ------------------------------------------------------------

no_bias <- glm(Y ~ X + C, family = binomial(link = "logit"), data = df) #model with no bias
exp(coef(no_bias)[2])
exp(coef(no_bias)[2] + summary(no_bias)$coef[2, 2] * qnorm(.025))
exp(coef(no_bias)[2] + summary(no_bias)$coef[2, 2] * qnorm(.975))
#ORyx = 2.04 (1.99, 2.09)

mc_bias <- glm(Y ~ Xstar + C, family = binomial(link = "logit"), data = df) #model with misclassification bias
exp(coef(mc_bias)[2])
exp(coef(mc_bias)[2] + summary(mc_bias)$coef[2, 2] * qnorm(.025))
exp(coef(mc_bias)[2] + summary(mc_bias)$coef[2, 2] * qnorm(.975))
#ORyx = 1.27 (1.24, 1.31)

# Model P(X=1|Xstar,C,Y) and obtain regression parameters -----------------------------------------------

x_model <- glm(X ~ Xstar + C + Y, family = binomial(link = "logit"), data = df)
X0 <- coef(x_model)[1]
X1 <- coef(x_model)[2]
X2 <- coef(x_model)[3]
X3 <- coef(x_model)[4]

# Create bootstrap function --------------------------------------------------------------------------

fun <- function (x1_0, x1_xstar, x1_c, x1_y) {
  set.seed(1234)
  est <- vector()
  nreps <- 10 #can vary number of bootstrap samples
  
  for(i in 1:nreps){
    bdf <- df[sample(1:nrow(df), n, replace = TRUE), ] #random samping with replacement
    
    pX <- plogis(x1_0 + x1_xstar * bdf$Xstar + x1_c * bdf$C + x1_y * bdf$Y) #model the probability of X
    
    combined <- bdf[rep(seq_len(nrow(bdf)), 2), ] #duplicate data
    combined$Xsim <- rep(c(1, 0), each = n) #Xsim=1 in first copy, Xsim=0 in second copy
    combined$pX <- c(pX, 1 - pX) #when Xsim=1, pX is prob of X=1; when Xsim=0, pX is prob of X=0
    
    final <- glm(Y ~ Xsim + C, family = binomial(link = "logit"), weights = pX, data = combined)
    est[i] <- coef(final)[2]
  }
  
  out <- list(exp(median(est)), exp(quantile(est, c(.025, .975))), hist(exp(est)))
  return(out)
}

# Apply functions ------------------------------------------------------------------------------------

# Known, correct bias parameters 

fun(x1_0 = X0, x1_xstar = X1, x1_c = X2, x1_y = X3)
# ORyx = 2.04 (2.03, 2.05)

# Incorret bias parameters

# Double each bias parameter
fun(x1_0 = 2*X0, x1_xstar = 2*X1, x1_c = 2*X2, x1_y = 2*X3)
# ORyx = 2.69 (2.67, 2.71)





