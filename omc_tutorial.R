# Adjusting for Outcome Misclassification

library(tidyverse)

set.seed(1234)
n <- 100000

# DERIVE DATA
c <- rbinom(n, 1, 0.5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c))
ystar <- rbinom(n, 1, plogis(-1 + log(1.25) * x + log(5) * y))

df <- data.frame(X = x, Y = y, C = c, Ystar = ystar)

rm(c1, x, y, ystar)

# INSPECT MODELS
nobias_model <- glm(Y ~ X + C,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 1.97 (1.87, 2.07)

bias_model <- glm(Ystar ~ X + C,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.44 (1.39, 1.50)

# OBTAIN BIAS PARAMETERS
y_model <- glm(Y ~ X + Ystar + C1,
               data = df,
               family = binomial(link = "logit"))
summary(y_model)

# ADJUST
# imputation approach
set.seed(1234)
est <- vector()
nreps <- 10

for (i in 1:nreps) {
  bdf <- df[sample(seq_len(n), n, replace = TRUE), ]

  bdf$Ypred <- rbinom(
    n, 1, plogis(coef(y_model)[1] + coef(y_model)[2] * bdf$X +
                   coef(y_model)[3] * bdf$Ystar + coef(y_model)[4] * bdf$C1)
  )

  final <- glm(Ypred ~ X + C1, family = binomial(link = "logit"), data = bdf)
  est[i] <- exp(coef(final)[2])
}

round(median(est), 2)
round(quantile(est, c(.025, .975)), 2)
# 1.98 (1.88, 2.07)

hist(exp(est))

# weighting approach
set.seed(1234)
est <- vector()
nreps <- 10

for (i in 1:nreps) {
  bdf <- df[sample(seq_len(n), n, replace = TRUE), ]

  y_probability <- plogis(
    coef(y_model)[1] + coef(y_model)[2] * bdf$X +
      coef(y_model)[3] * bdf$Ystar + coef(y_model)[4] * bdf$C1
  )

  combined <- bind_rows(bdf, bdf)
  combined$Ybar <- rep(c(1, 0), each = n)
  combined$y_weight <- c(y_probability, 1 - y_probability)

  final <- glm(Ybar ~ X + C1, family = binomial(link = "logit"),
               data = combined, weights = y_weight)
  est[i] <- exp(coef(final)[2])
}

round(median(est), 2)
round(quantile(est, c(.025, .975)), 2)
# 1.96 (1.94, 1.98)

hist(exp(est))
