library(multibias)
library(tidyverse)
library(nhanesA)
library(haven)

save_path <- '/Users/pbrendelprojects/analyses/'
data_path <- '/Users/pbrendel/data/nhanes/'


# load data ----
# https://wwwn.cdc.gov/nchs/nhanes/search/DataPage.aspx?Component=Dietary&Cycle=2021-2023

nhanesTables("Q", year = 2013)

# demographic data
# SEQN, RIAGENDR (gender), RIDAGEYR (age at interview), RIDRETH1 (race/hispanic origin),
# DMDEDUC2 (education), INDHHIN2 (household income), SDDSRVYR (Data release cycle)
df_demo_1314 <- nhanes('DEMO_H')
df_demo_1516 <- nhanes('DEMO_I')
df_demo_1718 <- nhanes('DEMO_J')

# nutrition data
# SEQN, DR1TALCO (alcohol, gm), DR1DRSTZ (recall status)
# WTDRD1 - Dietary day one sample weight
# WTDR2D - Dietary two-day sample weight
df_nutrition_1314_day1 <- nhanes('DR1TOT_H') |>
  rename(alcohol = DR1TALCO, recall_status = DR1DRSTZ)
df_nutrition_1314_day2 <- nhanes('DR2TOT_H') |>
  rename(alcohol = DR2TALCO, recall_status = DR2DRSTZ)
df_nutrition_1516_day1 <- nhanes('DR1TOT_I') |>
  rename(alcohol = DR1TALCO, recall_status = DR1DRSTZ)
df_nutrition_1516_day2 <- nhanes('DR2TOT_I') |>
  rename(alcohol = DR2TALCO, recall_status = DR2DRSTZ)
df_nutrition_1718_day1 <- nhanes('DR1TOT_J') |>
  rename(alcohol = DR1TALCO, recall_status = DR1DRSTZ)
df_nutrition_1718_day2 <- nhanes('DR2TOT_J') |>
  rename(alcohol = DR2TALCO, recall_status = DR2DRSTZ)

# alcohol data
# ALQ120Q - How often drink alcohol over past 12 mos
# ALQ121 - Past 12 mo how often have alcohol drink
df_alcohol_1314 <- nhanes('ALQ_H') |>
  select(SEQN, ALQ120Q) |>
  rename(alcohol_12mo = ALQ120Q)
df_alcohol_1516 <- nhanes('ALQ_I') |>
  select(SEQN, ALQ120Q) |>
  rename(alcohol_12mo = ALQ120Q)
df_alcohol_1718 <- nhanes('ALQ_J') |>
  select(SEQN, ALQ121) |>
  rename(alcohol_12mo = ALQ121) |>
  mutate(
    alcohol_12mo = case_when(
      alcohol_12mo == 'Never in the last year' ~ 0,
      alcohol_12mo == 'Every day' ~ 365,
      alcohol_12mo == 'Nearly every day' ~ 0,
      alcohol_12mo == '3 to 4 times a week' ~ 4*52,
      alcohol_12mo == '2 times a week' ~ 2*52,
      alcohol_12mo == 'Once a week' ~ 52,
      alcohol_12mo == '2 to 3 times a month' ~ 3*12,
      alcohol_12mo == 'Once a month' ~ 12,
      alcohol_12mo == '7 to 11 times in the last year' ~ 11,
      alcohol_12mo == '3 to 6 times in the last year' ~ 6,
      alcohol_12mo == '1 to 2 times in the last year' ~ 2,
      TRUE ~ NA
    )
  )
# the 17-18 data shows MANY more every-day drinkers
summary(df_alcohol_1718$alcohol_12mo)

# smoking data
# SMQ020 - Smoked at least 100 cigarettes in life
# SMQ040 - Do you now smoke cigarettes?
df_smoking_1314 <- nhanes('SMQ_H')
df_smoking_1516 <- nhanes('SMQ_I')
df_smoking_1718 <- nhanes('SMQ_J')

# mortality data
# seqn, eligstat (eligibility), mortstat (mortality status), ucod_leading
# ucod_leading of 2 = Malignant neoplasms (C00-C97)
# https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/
# https://ehsanx.github.io/EpiMethods/accessing7.html

df_mortality_1314 <- read_fwf(
  file = paste0(data_path, 'NHANES_2013_2014_MORT_2019_PUBLIC.dat'),
  col_types = "iiiiiiii",
  fwf_cols(SEQN = c(1, 6),
           eligstat = c(15, 15),
           mortstat = c(16, 16),
           ucod_leading = c(17, 19),
           diabetes = c(20, 20),
           hyperten = c(21, 21),
           permth_int = c(43, 45),
           permth_exm = c(46, 48)),
  na = c("", "."))

df_mortality_1516 <- read_fwf(
  file = paste0(data_path, 'NHANES_2015_2016_MORT_2019_PUBLIC.dat'),
  col_types = "iiiiiiii",
  fwf_cols(SEQN = c(1,6),
           eligstat = c(15,15),
           mortstat = c(16,16),
           ucod_leading = c(17,19),
           diabetes = c(20,20),
           hyperten = c(21,21),
           permth_int = c(43,45),
           permth_exm = c(46,48)),
  na = c("", "."))

df_mortality_1718 <- read_fwf(
  file = paste0(data_path, 'NHANES_2017_2018_MORT_2019_PUBLIC.dat'),
  col_types = "iiiiiiii",
  fwf_cols(SEQN = c(1,6),
           eligstat = c(15,15),
           mortstat = c(16,16),
           ucod_leading = c(17,19),
           diabetes = c(20,20),
           hyperten = c(21,21),
           permth_int = c(43,45),
           permth_exm = c(46,48)),
  na = c("", "."))

# join, prep data ----

df_demo <- df_demo_1314 |>
  bind_rows(df_demo_1516) |>
  bind_rows(df_demo_1718) |>
  select(SEQN, SDDSRVYR, RIAGENDR, RIDAGEYR, RIDRETH1, DMDEDUC2, INDHHIN2) |>
  rename(release = SDDSRVYR,
         gender = RIAGENDR,
         age = RIDAGEYR,
         race = RIDRETH1,
         education = DMDEDUC2,
         income = INDHHIN2
         )

# daily alcohol (averaged over two days)
df_nutrition <- select(df_nutrition_1314_day1, SEQN, alcohol, recall_status, WTDRD1, WTDR2D) |>
  bind_rows(select(df_nutrition_1314_day2, SEQN, alcohol, recall_status, WTDRD1, WTDR2D)) |>
  bind_rows(select(df_nutrition_1516_day1, SEQN, alcohol, recall_status, WTDRD1, WTDR2D)) |>
  bind_rows(select(df_nutrition_1516_day2, SEQN, alcohol, recall_status, WTDRD1, WTDR2D)) |>
  bind_rows(select(df_nutrition_1718_day1, SEQN, alcohol, recall_status, WTDRD1, WTDR2D)) |>
  bind_rows(select(df_nutrition_1718_day2, SEQN, alcohol, recall_status, WTDRD1, WTDR2D)) |>
  filter(recall_status == 'Reliable and met the minimum criteria') |>
  rename(weight_day1 = WTDRD1,
         weight_day2 = WTDR2D
         ) |>
  group_by(SEQN, weight_day1, weight_day2) |>
  summarize(alcohol = mean(alcohol))

df_alcohol <- df_alcohol_1314 |>
  bind_rows(df_alcohol_1516) |>
  bind_rows(df_alcohol_1718)

df_smoking <- df_smoking_1314 |>
  bind_rows(df_smoking_1516) |>
  bind_rows(df_smoking_1718) |>
  select(SEQN, SMQ020) |>
  rename(smoked_100_cigs = SMQ020)

df_mortality <- df_mortality_1314 |>
  bind_rows(df_mortality_1516) |>
  bind_rows(df_mortality_1718) |>
  select(SEQN, eligstat, mortstat, ucod_leading)

df <- df_demo |>
  inner_join(df_nutrition, by = 'SEQN') |>
  inner_join(df_alcohol, by = 'SEQN') |>
  inner_join(df_smoking, by = 'SEQN') |>
  inner_join(df_mortality, by = 'SEQN') |>
  mutate(gender_female = if_else(gender == "Female", 1, 0)) |>
  mutate(no_alcohol_12mo = if_else(alcohol_12mo == 0, 1, 0)) |>
  mutate(smoked_100cigs = case_when(
    smoked_100_cigs == "Yes" ~ 1,
    smoked_100_cigs == "No" ~ 0,
    TRUE ~ NA)
  ) |>
  rename(seqn = SEQN,
         alcohol_day_total = alcohol) |>
  select(seqn, release, gender_female, age, race, education, income,
         weight_day1, weight_day2, alcohol_day_total, alcohol_12mo,
         smoked_100cigs, eligstat, mortstat, ucod_leading
         )

write_csv(df, paste0(data_path, 'nhanes.csv'))

df_filtered <- df |>
  mutate(alcohol_extreme = case_when(
    alcohol > 14*1.5 & gender_female == 1 ~ 1,
    alcohol > 28*1.5 & gender_female == 0 ~ 1,
    TRUE ~ 0)
    ) |>
  filter(age >= 18) |>
  filter(eligstat == 1) |>
  filter(alcohol_12mo > 0)

# inspect data ----
df_filtered |>
  group_by(alcohol_extreme, mortstat) |>
  summarize(count = n()) |>
  ungroup() |>
  mutate(proportion = count / sum(count))

table(df_filtered$mortstat)
table(df_filtered$ucod_leading)

summary(df$alcohol)
summary(df_filtered$alcohol)

table(df_filtered$alcohol_extreme, df_filtered$mortstat)
summary(df$age)
table(df_filtered$gender_female, useNA = "ifany")

table(df$race, useNA = "ifany")
table(df$education, useNA = "ifany")
table(df$income, useNA = "ifany")

table(df_filtered$smoked_100_cigs, useNA = "ifany")

# no age/sex pattern among extreme alcohol
ggplot(df_filtered,
       aes(x = as.factor(alcohol_extreme), y = age, color = as.factor(gender_female))
       ) +
  geom_boxplot()

# no age/sex pattern among extreme alcohol
ggplot(df_filtered,
       aes(x = as.factor(mortstat), y = age, color = as.factor(gender_female))
) +
  geom_boxplot()


# multibias ----
# selection bias:
# creating less opportunity for cancer-specific
# exposure misclassification: under-reporting alcohol use
# uncontrolled confounding: smoking (many missing values)

base_mod <- glm(mortstat ~ alcohol_extreme + gender_female + age,
                data = df_filtered,
                family = binomial(link = 'logit')
                )

summary(base_mod)
exp(coef(base_mod))
exp(confint(base_mod))

# 1
set.seed(1234)

df_observed <- data_observed(
  data = df_filtered,
  exposure = 'alcohol_extreme',
  outcome = 'mortstat',
  confounders = c('age', 'gender_female')
)

df_temp <- df_filtered |>
  filter(!is.na(smoked_100_cigs))

df_val1 <- data_validation(
  data = df_temp,
  true_exposure = 'alcohol_extreme',
  true_outcome = 'mortstat',
  confounders = c('age', 'gender_female', 'smoked_100_cigs')
)

adjust_uc(df_observed, df_val1)

# 2
set.seed(1234)

df_temp2 <- df_filtered |>
  filter(!is.na(smoked_100_cigs)) |>
  mutate(alcohol_adj = if_else(
    income == '$100,000 and Over' | education == 'College graduate or above',
    alcohol * 1.5,
    alcohol
    )
  ) |>
  mutate(alcohol_extreme_adj = case_when(
    alcohol_adj > 14*2 & gender_female == 1 ~ 1,
    alcohol_adj > 14*4 & gender_female == 0 ~ 1,
    TRUE ~ 0
    )
  )

table(df_temp2$alcohol_extreme, df_temp2$alcohol_extreme_adj)

df_val2 <- data_validation(
  data = df_temp2,
  true_exposure = 'alcohol_extreme_adj',
  true_outcome = 'mortstat',
  confounders = c('age', 'gender_female', 'smoked_100_cigs'),
  misclassified_exposure = 'alcohol_extreme'
)

adjust_uc_em(df_observed, df_val2)

# 3

df_temp3 <- df_filtered |>
  filter(!is.na(smoked_100_cigs)) |>
  mutate(alcohol_adj = if_else(
    income == '$100,000 and Over' | education == 'College graduate or above',
    alcohol * 1.5,
    alcohol
    )
  ) |>
  mutate(alcohol_extreme_adj = case_when(
    alcohol_adj > 14*2 & gender_female == 1 ~ 1,
    alcohol_adj > 14*4 & gender_female == 0 ~ 1,
    TRUE ~ 0
    )
  ) |>
  mutate(
    weight = if_else(WTDR2D == 0, WTDRD1, WTDR2D)
  )

# weight_sum <- sum(df_temp3$weight)
# df_temp3$pS <- df_temp3$weight / weight_sum

# get the selected subjects: sample with replacement until reach original N
selected_sample <- sample(
  x = df_temp3$SEQN,
  size = nrow(df_temp3),
  replace = TRUE,
  prob = df_temp3$weight
)

df_selected_sample <- data.frame()
for (id in selected_sample) {
  df_selected_sample <- rbind(df_selected_sample, df_temp3[df_temp3$SEQN == id, ])
}

# get the un-selected subjects
not_selected_sample <- df_temp3$SEQN[!(df_temp3$SEQN %in% selected_sample)]
df_not_selected_sample <- df_temp3[df_temp3$SEQN %in% not_selected_sample, ]

# check
length(not_selected_sample) + length(unique(selected_sample)) == nrow(df_temp3)

# make df
df_selected_sample$selection <- 1
df_not_selected_sample$selection <- 0
df_temp3b <- rbind(df_selected_sample, df_not_selected_sample)


df_val3 <- data_validation(
  data = df_temp3b,
  true_exposure = 'alcohol_extreme_adj',
  true_outcome = 'mortstat',
  confounders = c('age', 'gender_female', 'smoked_100_cigs'),
  misclassified_exposure = 'alcohol_extreme',
  selection = 'selection'
)

adjust_uc_em_sel(df_observed, df_val3)
