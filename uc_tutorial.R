# Adjusting for Uncontrolled Confounding

library(tidyverse)

# DERIVE DATA

set.seed(1234)
n <- 100000

c <- rbinom(n, 1, 0.5)
u <- rbinom(n, 1, 0.5)
x <- rbinom(n, 1, plogis(-0.5 + 0.5 * c + 1.5 * u))
y <- rbinom(n, 1, plogis(-0.5 + log(2) * x + 0.5 * c + 1.5 * u))

df <- data.frame(X = x, Y = y, C = c, U = u)
rm(x, y, c, u)

# INSPECT MODELS

nobias_model <- glm(Y ~ X + C + U,
                    family = binomial(link = "logit"),
                    data = df)
exp(coef(nobias_model)[2])
c(exp(coef(nobias_model)[2] + summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(coef(nobias_model)[2] + summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 2.02 (1.96, 2.09)

biased_model <- glm(Y ~ X + C,
                    family = binomial(link = "logit"),
                    data = df)
exp(coef(biased_model)[2])
c(exp(coef(biased_model)[2] + summary(biased_model)$coef[2, 2] * qnorm(.025)),
  exp(coef(biased_model)[2] + summary(biased_model)$coef[2, 2] * qnorm(.975)))
# 3.11 (3.03, 3.20)

# OBTAIN BIAS PARAMETERS

u_model <- glm(U ~ X + Y + C,
               family = binomial(link = "logit"),
               data = df)
summary(u_model)

# ADJUST

# imputation approach
adjust_uc_imp_loop <- function(
  coef_0, coef_x, coef_c, coef_y, nreps, plot = FALSE
) {

  est <- vector()
  for (i in 1:nreps){
    # bootstrap sample
    bdf <- df[sample(seq_len(n), n, replace = TRUE), ]

    # impute u
    bdf$Upred <- rbinom(n, 1, plogis(coef_0 + coef_x * bdf$X +
                                       coef_c * bdf$C + coef_y * bdf$Y))

    final_model <- glm(Y ~ X + C + Upred,
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
correct_results <- adjust_uc_imp_loop(
  coef_0 = coef(u_model)[1],
  coef_x = coef(u_model)[2],
  coef_c = coef(u_model)[3],
  coef_y = coef(u_model)[4],
  nreps = 10
)

correct_results$estimate
correct_results$ci
# 2.03 (1.95, 2.10)

# using incorret bias parameters
set.seed(1234)
incorrect_results <- adjust_uc_imp_loop(
  coef_0 = coef(u_model)[1] * 2,
  coef_x = coef(u_model)[2] * 2,
  coef_c = coef(u_model)[3] * 2,
  coef_y = coef(u_model)[4] * 2,
  nreps = 10
)

incorrect_results$estimate
incorrect_results$ci
# 0.81 (0.78, 0.84)

# weighting approach
adjust_uc_wgt_loop <- function(
  coef_0, coef_x, coef_c, coef_y, nreps, plot = FALSE
) {

  est <- vector()
  for (i in 1:nreps){
    # bootstrap sample
    bdf <- df[sample(seq_len(n), n, replace = TRUE), ]

    # the probability of U for each observation
    prob_u <- plogis(coef_0 + coef_x * bdf$X + coef_c * bdf$C + coef_y * bdf$Y)

    # duplicate data
    combined <- dplyr::bind_rows(bdf, bdf)
    # assign values for U
    combined$Ubar <- rep(c(1, 0), each = n)
    # create weight
    # when Ubar=1, u_weight=P(U=1); when Ubar=0, u_weight=P(U=0)
    combined$u_weight <- c(prob_u, 1 - prob_u)

    final_model <- glm(Y ~ X + C + Ubar, family = binomial(link = "logit"),
                       data = combined, weights = combined$u_weight)
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
correct_results <- adjust_uc_wgt_loop(
  coef_0 = coef(u_model)[1],
  coef_x = coef(u_model)[2],
  coef_c = coef(u_model)[3],
  coef_y = coef(u_model)[4],
  nreps = 10
)

correct_results$estimate
correct_results$ci
# 2.03 (1.98, 2.07)

# using incorret bias parameters
set.seed(1234)
incorrect_results <- adjust_uc_wgt_loop(
  coef_0 = coef(u_model)[1] * 2,
  coef_x = coef(u_model)[2] * 2,
  coef_c = coef(u_model)[3] * 2,
  coef_y = coef(u_model)[4] * 2,
  nreps = 10
)

incorrect_results$estimate
incorrect_results$ci
# 0.81 (0.79, 0.82)