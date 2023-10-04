# Adjusting for Exposure Misclassification & Outcome Misclassification

library(dplyr)
library(nnet)

# DERIVE DATA

set.seed(1234)
n <- 100000

c <- rbinom(n, 1, .5)
x  <- rbinom(n, 1, plogis(-.5 + .5 * c))
y  <- rbinom(n, 1, plogis(-.5 + log(2) * x + .5 * c))
xstar <- rbinom(n, 1, plogis(-1 + log(5) * x + log(1.25) * y))
ystar <- rbinom(n, 1, plogis(-1 + log(1.25) * x + log(5) * y))

df <- data.frame(X = x, Y = y, C = c, Xstar = xstar, Ystar = ystar)
rm(c, x, y, xstar, ystar)

# INSPECT MODELS

nobias_model <- glm(Y ~ X + C,
                    family = binomial(link = "logit"),
                    data = df)

exp(summary(nobias_model)$coef[2, 1])
c(exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(nobias_model)$coef[2, 1] +
        summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 2.04 (1.99, 2.09)

bias_model <- glm(Ystar ~ Xstar + C,
                  family = binomial(link = "logit"),
                  data = df)

exp(summary(bias_model)$coef[2, 1])
c(exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(summary(bias_model)$coef[2, 1] +
        summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.26 (1.23, 1.29)

# OBTAIN BIAS PARAMETERS

x_model <- glm(X ~ Xstar + Ystar + C, #???
               data = df,
               family = binomial(link = "logit"))
y_model <- glm(Y ~ X + Ystar + C, #???
               data = df,
               family = binomial(link = "logit"))

xy_model <- multinom(
  paste(X, Y) ~ Xstar + Ystar + C,
  data = df
)
summary(xy_model)$coefficients

# multinomial approach
adjust_multinom_emc_omc_loop <- function(
  x1y0_model_coefs,
  x0y1_model_coefs,
  x1y1_model_coefs,
  nreps,
  plot = FALSE
) {

  est <- vector()
  for (i in 1:nreps) {

    bdf <- df[sample(seq_len(n), n, replace = TRUE), ]

    p_x1y0 <- exp(
      x1y0_model_coefs[1] +
        x1y0_model_coefs[2] * bdf$Xstar +
        x1y0_model_coefs[3] * bdf$Ystar +
        x1y0_model_coefs[4] * bdf$C
    )
    p_x0y1 <- exp(
      x0y1_model_coefs[1] +
        x0y1_model_coefs[2] * bdf$Xstar +
        x0y1_model_coefs[3] * bdf$Ystar +
        x0y1_model_coefs[4] * bdf$C
    )
    p_x1y1 <- exp(
      x1y1_model_coefs[1] +
        x1y1_model_coefs[2] * bdf$Xstar +
        x1y1_model_coefs[3] * bdf$Ystar +
        x1y1_model_coefs[4] * bdf$C
    )

    denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

    x0y0_pred <- 1 / denom
    x1y0_pred <- p_x1y0 / denom
    x0y1_pred <- p_x0y1 / denom
    x1y1_pred <- p_x1y1 / denom

    df_xy_pred <- data.frame(
      X0Y0 = x0y0_pred,
      X1Y0 = x1y0_pred,
      X0Y1 = x0y1_pred,
      X1Y1 = x1y1_pred
    )
    df_xy_pred4 <- bind_rows(
      df_xy_pred,
      df_xy_pred,
      df_xy_pred,
      df_xy_pred
    )

    combined <- bind_rows(bdf, bdf, bdf, bdf) %>%
      bind_cols(df_xy_pred4) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pXY = case_when(Xbar == 0 & Ybar == 0 ~ X0Y0,
                             Xbar == 1 & Ybar == 0 ~ X1Y0,
                             Xbar == 0 & Ybar == 1 ~ X0Y1,
                             Xbar == 1 & Ybar == 1 ~ X1Y1))

    final_model <- glm(Ybar ~ Xbar + C,
                       family = binomial(link = "logit"),
                       data = combined,
                       weights = combined$pXY)
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

set.seed(1234)
correct_results <- adjust_multinom_emc_omc_loop(
  x1y0_model_coefs = c(
    summary(xy_model)$coefficients[2, 1],
    summary(xy_model)$coefficients[2, 2],
    summary(xy_model)$coefficients[2, 3],
    summary(xy_model)$coefficients[2, 4]
  ),
  x0y1_model_coefs = c(
    summary(xy_model)$coefficients[1, 1],
    summary(xy_model)$coefficients[1, 2],
    summary(xy_model)$coefficients[1, 3],
    summary(xy_model)$coefficients[1, 4]
  ),
  x1y1_model_coefs = c(
    summary(xy_model)$coefficients[3, 1],
    summary(xy_model)$coefficients[3, 2],
    summary(xy_model)$coefficients[3, 3],
    summary(xy_model)$coefficients[3, 4]
  ),
  nreps = 100
)

correct_results$estimate
correct_results$ci
# 2.04 (2.03, 2.05)

# two model approach
adjust_emc_omc_loop <- function(
  x_model_coefs,
  y_model_coefs,
  nreps,
  plot = FALSE
) {

  est <- vector()
  for (i in 1:nreps) {

    bdf <- df[sample(seq_len(n), n, replace = TRUE), ]

    bdf$Xpred <- rbinom(
      n,
      1,
      plogis(x_model_coefs[1] + x_model_coefs[2] * bdf$Xstar +
               x_model_coefs[3] * bdf$Ystar + x_model_coefs[4] * bdf$C) #???
    )
    bdf$Ypred <- rbinom(
      n,
      1,
      plogis(y_model_coefs[1] + y_model_coefs[2] * bdf$Xpred +
               y_model_coefs[3] * bdf$Ystar + y_model_coefs[4] * bdf$C) #???
    )

    final_model <- glm(Ypred ~ Xpred + C,
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

set.seed(1234)
correct_results <- adjust_emc_omc_loop(
  x_model_coefs = c(
    x_model$coef[1],
    x_model$coef[2],
    x_model$coef[3],
    x_model$coef[4]
  ),
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4]
  ),
  nreps = 100
)

correct_results$estimate
correct_results$ci
# 2.05 (1.99, 2.09)