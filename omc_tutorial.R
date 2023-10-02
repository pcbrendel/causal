# Adjusting for Outcome Misclassification

library(tidyverse)

# DERIVE DATA

set.seed(1234)
n <- 100000

c <- rbinom(n, 1, 0.5)
x <- rbinom(n, 1, plogis(-2 + log(1.5) * c))
y <- rbinom(n, 1, plogis(-2.5 + log(2) * x + log(1.5) * c))
ystar <- rbinom(n, 1, plogis(-1 + log(1.25) * x + log(5) * y))

df <- data.frame(X = x, Y = y, C = c, Ystar = ystar)
rm(c, x, y, ystar)

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
# 1.44 (1.39, 1.49)

# OBTAIN BIAS PARAMETERS

y_model <- glm(Y ~ X + Ystar + C,
               data = df,
               family = binomial(link = "logit"))
summary(y_model)

# ADJUST

# imputation approach
adjust_omc_imp_loop <- function(
  coef_0, coef_x, coef_ystar, coef_c, nreps, plot = FALSE
) {

  est <- vector()
  for (i in 1:nreps) {

    bdf <- df[sample(seq_len(n), n, replace = TRUE), ]


    bdf$Ypred <- rbinom(
      n, 1, plogis(coef_0 + coef_x * bdf$X +
                     coef_ystar * bdf$Ystar + coef_c * bdf$C)
    )

    final_model <- glm(Ypred ~ X + C,
                       family = binomial(link = "logit"),
                       data = bdf)
    est[i] <- exp(coef(final_model)[2])
  }

  out <- list(
    estimate = round(median(est), 2),
    ci = round(quantile(est, c(.025, .975)), 2)
  )

  if (plot) {
    out$hist <- hist(exp(est))
  }

  return(out)

}

# using known, correct bias parameters
set.seed(1234)
correct_results <- adjust_omc_imp_loop(
  coef_0 =     coef(y_model)[1],
  coef_x =     coef(y_model)[2],
  coef_ystar = coef(y_model)[3],
  coef_c =     coef(y_model)[4],
  nreps = 10
)

correct_results$estimate
correct_results$ci
# 1.98 (1.88, 2.07)

# using incorret bias parameters
set.seed(1234)
incorrect_results <- adjust_omc_imp_loop(
  coef_0 =     coef(y_model)[1] * 2,
  coef_x =     coef(y_model)[2] * 2,
  coef_ystar = coef(y_model)[3] * 2,
  coef_c =     coef(y_model)[4] * 2,
  nreps = 10
)

incorrect_results$estimate
incorrect_results$ci
# 3.62 (3.35, 3.87)

# weighting approach
adjust_omc_wgt_loop <- function(
  coef_0, coef_x, coef_ystar, coef_c, nreps, plot = FALSE
) {

  est <- vector()
  for (i in 1:nreps) {

    bdf <- df[sample(seq_len(n), n, replace = TRUE), ]

    y_probability <- plogis(
      coef_0 + coef_x * bdf$X + coef_ystar * bdf$Ystar + coef_c * bdf$C
    )

    combined <- dplyr::bind_rows(bdf, bdf)
    combined$Ybar <- rep(c(1, 0), each = n)
    combined$y_weight <- c(y_probability, 1 - y_probability)

    final_model <- glm(Ybar ~ X + C, family = binomial(link = "logit"),
                       data = combined, weights = combined$y_weight)
    est[i] <- exp(coef(final_model)[2])
  }

  out <- list(
    estimate = round(median(est), 2),
    ci = round(quantile(est, c(.025, .975)), 2)
  )

  if (plot) {
    out$hist <- hist(exp(est))
  }

  return(out)

}

# using known, correct bias parameters
set.seed(1234)
correct_results <- adjust_omc_wgt_loop(
  coef_0 =     coef(y_model)[1],
  coef_x =     coef(y_model)[2],
  coef_ystar = coef(y_model)[3],
  coef_c =     coef(y_model)[4],
  nreps = 10
)

correct_results$estimate
correct_results$ci
# 1.96 (1.94, 1.98)

# using incorret bias parameters
set.seed(1234)
incorrect_results <- adjust_omc_wgt_loop(
  coef_0 =     coef(y_model)[1] * 2,
  coef_x =     coef(y_model)[2] * 2,
  coef_ystar = coef(y_model)[3] * 2,
  coef_c =     coef(y_model)[4] * 2,
  nreps = 10
)

incorrect_results$estimate
incorrect_results$ci
# 3.66 (3.60, 3.73)
