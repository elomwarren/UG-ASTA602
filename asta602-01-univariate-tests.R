# Departement of Statistics and Actuarial Science, University of Ghana
# asta602-univariate-tests.R
# APPLIED MULTIVARIATE ANALYSIS
# @author: Elom Kodjoh-Kpakpassou
# Email: elwkodjoh-kpakpassou@st.ug.edu.gh


# Clear global environment
rm(list = ls())

# For reproducibility
set.seed(123)

# install.packages("misty") # Z-test and T-test
# install.packages("car") # levene's test
# install.packages("nortest") # normality test


# ######################
# I - UNIVARIATE TESTS #
# ######################

##########
# IMPORTS
library(nortest)
library(misty)
library(car)

##########################
# --- 0. Normality tests #
##########################

###########################
# --- 0.1 Shapiro-Walk test
###########################

# DATA: mtcars {datasets}

# Visually

hist(mtcars$mpg)

qqnorm(mtcars$mpg)
qqline(mtcars$mpg,
    distribution = qnorm,
    probs = c(0.25, 0.75), qtype = 7
)


# Small to moderate samples - 3 < n < 5000

shapiro.test(mtcars$mpg)

#############################################
# --- 0.2 Anderson Darling test for normality
#############################################
nortest::ad.test(mtcars$mpg)

#############################################
# --- 0.3 Cramer-von Mises test for normality
#############################################
nortest::cvm.test(mtcars$mpg)

# ...

#######################################
# --- 1. Test of equality of variance #
#######################################

# H0: Equality of variance

####################################
# --- 1.1 Bartlett Test - Parametric
####################################
# Assumption: Normality

# DATA: InsectSprays {datasets}

plot(count ~ spray, data = InsectSprays)
bartlett.test(InsectSprays$count, InsectSprays$spray)
# bartlett.test(count ~ spray, data = InsectSprays) # same result

######################################
# --- 1.2 Levene Test - Non-parametric
######################################

# Levene’s test is less sensitive to departures from normal distribution than the Bartlett’s test.

car::leveneTest(count ~ spray, data = InsectSprays)

##########################################################################
# --- 2. One Sample / Two Samples / Paired Samples Mean Comparison Tests #
##########################################################################

################
# --- 2.1 Z-test
################
# ASSUMPTIONS
# Normality (small sample)
# All data points must be independent.
# Population variance known or
# n >= 30 and unknown population variance (sample variance as estimate of the population variance)

#############################
# --- 2.1.1 ONE SAMPLE Z-TEST
#############################

# NOTE: alternative = "two.sided" is default

# Two-sided one-sample z-test, population mean = 20, population SD = 6
misty::test.z(mtcars$mpg, sigma = 6, mu = 20)


# p=0.932 > 0.05: We fail to reject the null hypothesis.
# We do not have enough evidence to say the sample mean differs from the population mean.

#############################
# --- 2.1.2 TWO-SAMPLE Z-TEST
#############################

# Two-sided two-sample z-test, population SD = 6, equal SD assumption
misty::test.z(mpg ~ vs, data = mtcars, sigma = 6)

# p=0.000 < 0.05: We reject the null hypothesis.
# The data suggest that there is a significant difference between the two sample means.

#################################
# --- 2.1.3 PAIRED SAMPLES Z-TEST
#################################

# Two-sided paired-sample z-test, population SD of difference score = 1.2
# IQ dataset

IQ1 <- c(127, 98, 105, 83, 133, 90, 107, 98, 91, 100, 88, 96, 110, 87, 88, 88, 105, 95, 79, 106)
IQ2 <- c(137, 108, 115, 93, 143, 100, 117, 108, 101, 110, 98, 106, 120, 97, 98, 100, 115, 111, 89, 116)


misty::test.z(IQ1, IQ2, sigma = 1.4, paired = TRUE)
# ...

################
# --- 2.2 T-test
################

# ASSUMPTIONS
# All data points must be independent.
# Unknown Population Variance
# Small Sample Size - n < 30

#############################
# --- 2.2.1 ONE SAMPLE T-TEST
#############################
# Two-sided one-sample t-test, population mean = 20
t.test(mtcars$mpg, mu = 20) # 'stats' package
misty::test.t(mtcars$mpg, mu = 20)
# ...

#############################
# --- 2.2.2 TWO-SAMPLE T-TEST
#############################
# Two-sided two-sample t-test
misty::test.t(mpg ~ vs, data = mtcars)
# ...

#################################
# --- 2.2.3 PAIRED-SAMPLE T-TEST
#################################

# Assumption
# Observed difference comes from a normal population

# Two-sided paired-sample t-test, population SD of difference score = 1.2
# IQ dataset
IQ1 <- c(127, 98, 105, 83, 133, 90, 107, 98, 91, 100, 88, 96, 110, 87, 88, 88, 105, 95, 79, 106)
IQ2 <- c(137, 108, 115, 93, 143, 100, 117, 108, 101, 110, 98, 106, 120, 97, 98, 100, 115, 111, 89, 116)

# We assume normality for this data

misty::test.t(IQ1, IQ2, paired = TRUE)
# ...

############################################################################################
# --- 3. One-Way ANalysis Of VAriance (ANOVA) - Multiple groups (samples) mean comparison #
############################################################################################

# https://www.r-bloggers.com/2020/10/anova-in-r/
# https://statsandr.com/blog/anova-in-r/#anova-in-r

# ASSUMPTIONS
# (Continuous outcome, One Categorical predictor with at least 2 levels)
# (Independence) - by design
# Normality of residuals
# Equality of variances (homoscedasticity)


# DATA: Penguins dataset {datasets}
# ?penguins # Description of the dataset
# outcome: flipper_len
# categorical predictor: species
#  research question “Is the length of the flippers different between the 3 species of penguins?”

summary(penguins)
pgdat <- penguins[, c("species", "flipper_len")]

###############################
# --- Residuals Normality check
res_aov <- aov(flipper_len ~ species,
    data = pgdat
)

# graphical layout - 1 row, 2 columns
par(mfrow = c(1, 2)) # combine plots

# histogram
hist(res_aov$residuals)

# QQ-plot
car::qqPlot(res_aov$residuals,
    id = FALSE # id = FALSE to remove point identification
)

# reset graphical layout
par(mfrow = c(1, 1))

#################################
# --- Equality of variances check

# Visually (informal)
boxplot(
    flipper_len ~ species,
    data = pgdat
)

# Levene's test
car::leveneTest(
    flipper_len ~ species,
    data = pgdat
)
# p-value = 0.7188: we fail to reject the null hypothesis

# ###### ANOVA test

# 1. oneway.test() function:
oneway.test(
    flipper_len ~ species,
    data = pgdat,
    var.equal = TRUE # assuming equal variances
)

# 2. With the summary() and aov() functions:
res_aov <- aov(
    flipper_length_mm ~ species,
    data = dat
) # run earlier

summary(res_aov)

##### INTERPRETATION
# Given that the p-value is smaller than 0.05, we reject the null hypothesis, so we reject the hypothesis that all means are equal. Therefore, we can conclude that at least one species is different than the others in terms of flippers length (p-value < 2.2e-16).
