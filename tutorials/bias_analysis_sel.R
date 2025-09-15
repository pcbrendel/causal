# Adjusting for Selection Bias

# DERIVE DATA

set.seed(1234)
n <- 100000

c <- rbinom(n, 1, .5)
x  <- rbinom(n, 1, plogis(-.5 + .5 * c))
y  <- rbinom(n, 1, plogis(-.5 + log(2) * x + .5 * c))
s <- rbinom(n, 1, plogis(-.5 + 1.5 * x + 1.5 * y))

df <- data.frame(X = x, Y = y, C = c, S = s)
rm(c, x, y, s)

# INSPECT MODELS

nobias_model <- glm(Y ~ X + C,
                    family = binomial(link = "logit"),
                    data = df)
exp(coef(nobias_model)[2])
c(exp(coef(nobias_model)[2] + summary(nobias_model)$coef[2, 2] * qnorm(.025)),
  exp(coef(nobias_model)[2] + summary(nobias_model)$coef[2, 2] * qnorm(.975)))
# 2.04 (1.99, 2.09)

bias_model <- glm(Y ~ X + C,
                  family = binomial(link = "logit"),
                  data = df[sample(seq_len(n), n, replace = TRUE, df$S), ])
exp(coef(bias_model)[2])
c(exp(coef(bias_model)[2] + summary(bias_model)$coef[2, 2] * qnorm(.025)),
  exp(coef(bias_model)[2] + summary(bias_model)$coef[2, 2] * qnorm(.975)))
# 1.34 (1.31, 1.38)

# OBTAIN BIAS PARAMETERS

s_model <- glm(S ~ X + Y, data = df, family = binomial(link = "logit"))
summary(s_model)

# ADJUST

adjust_sel_loop <- function(
  coef_0, coef_x, coef_y, nreps, plot = FALSE
) {

  est <- vector()
  for (i in 1:nreps){

    bdf <- df[sample(seq_len(n), n, replace = TRUE, df$S), ]

    prob_s <- plogis(coef_0 + coef_x * bdf$X + coef_y * bdf$Y)

    final_model <- glm(Y ~ X + C,
                       family = binomial(link = "logit"),
                       weights = (1 / prob_s),
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
correct_results <- adjust_sel_loop(
  coef_0 = coef(s_model)[1],
  coef_x = coef(s_model)[2],
  coef_y = coef(s_model)[3],
  nreps = 10
)

correct_results$estimate
correct_results$ci
# 2.05 (2.01, 2.07)

# using incorret bias parameters
set.seed(1234)
incorrect_results <- adjust_sel_loop(
  coef_0 = coef(s_model)[1] * 2,
  coef_x = coef(s_model)[2] * 2,
  coef_y = coef(s_model)[3] * 2,
  nreps = 10
)

incorrect_results$estimate
incorrect_results$ci
# 3.92 (3.84, 3.96)