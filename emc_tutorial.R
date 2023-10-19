# Adjusting for Exposure Misclassification

library(tidyverse)

# DERIVE DATA

set.seed(1234)
n <- 100000

c <- rbinom(n, 1, .5)
x  <- rbinom(n, 1, plogis(-.5 + .5 * c))
y  <- rbinom(n, 1, plogis(-.5 + log(2) * x + .5 * c))
xstar <- rbinom(n, 1, plogis(-1 + log(5) * x + log(1.25) * y))

df <- data.frame(X = x, Y = y, C = c, Xstar = xstar)
rm(c, x, y, xstar)

# INSPECT MODELS

nobias_model <- glm(Y ~ X + C, family = binomial(link = "logit"), data = df)
exp(coef(nobias_model)[2])
exp(coef(nobias_model)[2] + summary(nobias_model)$coef[2, 2] * qnorm(.025))
exp(coef(nobias_model)[2] + summary(nobias_model)$coef[2, 2] * qnorm(.975))
# 2.04 (1.99, 2.09)

bias_model <- glm(Y ~ Xstar + C, family = binomial(link = "logit"), data = df)
exp(coef(bias_model)[2])
exp(coef(bias_model)[2] + summary(bias_model)$coef[2, 2] * qnorm(.025))
exp(coef(bias_model)[2] + summary(bias_model)$coef[2, 2] * qnorm(.975))
# 1.61 (1.57, 1.65)

# OBTAIN BIAS PARAMETERS

x_model <- glm(X ~ Xstar + Y + C, family = binomial(link = "logit"), data = df)
summary(x_model)

# ADJUST

# imputation approach
adjust_emc_imp_loop <- function(
  coef_0, coef_xstar, coef_y, coef_c, nreps, plot = FALSE
) {

  est <- vector()
  for (i in 1:nreps) {

    bdf <- df[sample(seq_len(n), n, replace = TRUE), ]

    bdf$Xpred <- rbinom(
      n, 1, plogis(coef_0 + coef_xstar * bdf$Xstar +
                     coef_y * bdf$Y + coef_c * bdf$C)
    )

    final_model <- glm(Y ~ Xpred + C,
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
correct_results <- adjust_emc_imp_loop(
  coef_0 =     coef(x_model)[1],
  coef_xstar = coef(x_model)[2],
  coef_y =     coef(x_model)[3],
  coef_c =     coef(x_model)[4],
  nreps = 10
)

correct_results$estimate
correct_results$ci
# 2.04 (2.01, 2.08)

# using incorret bias parameters
set.seed(1234)
incorrect_results <- adjust_emc_imp_loop(
  coef_0 =     coef(x_model)[1] * 2,
  coef_xstar = coef(x_model)[2] * 2,
  coef_y =     coef(x_model)[3] * 2,
  coef_c =     coef(x_model)[4] * 2,
  nreps = 10
)

incorrect_results$estimate
incorrect_results$ci
# 2.85 (2.80, 2.89)

# weighting approach
adjust_emc_wgt_loop <- function(
  coef_0, coef_xstar, coef_y, coef_c, nreps, plot = FALSE
) {

  est <- vector()
  for (i in 1:nreps) {

    bdf <- df[sample(seq_len(n), n, replace = TRUE), ]

    x_probability <- plogis(
      coef_0 + coef_xstar * bdf$Xstar + coef_y * bdf$Y + coef_c * bdf$C
    )

    combined <- dplyr::bind_rows(bdf, bdf)
    combined$Xbar <- rep(c(1, 0), each = n)
    combined$x_weight <- c(x_probability, 1 - x_probability)

    final_model <- glm(Y ~ Xbar + C, family = binomial(link = "logit"),
                       data = combined, weights = combined$x_weight)
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
correct_results <- adjust_emc_wgt_loop(
  coef_0 =     coef(x_model)[1],
  coef_xstar = coef(x_model)[2],
  coef_y =     coef(x_model)[3],
  coef_c =     coef(x_model)[4],
  nreps = 10
)

correct_results$estimate
correct_results$ci
# 2.05 (2.04, 2.06)

# using incorret bias parameters
set.seed(1234)
incorrect_results <- adjust_emc_wgt_loop(
  coef_0 =     coef(x_model)[1] * 2,
  coef_xstar = coef(x_model)[2] * 2,
  coef_y =     coef(x_model)[3] * 2,
  coef_c =     coef(x_model)[4] * 2,
  nreps = 10
)

incorrect_results$estimate
incorrect_results$ci
# 2.85 (2.84, 2.88)