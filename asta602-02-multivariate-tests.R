# Department of Statistics and Actuarial Science, University of Ghana
# asta602-multivariate-tests.R
# APPLIED MULTIVARIATE ANALYSIS
# @author: Elom Kodjoh-Kpakpassou
# Email: elwkodjoh-kpakpassou@st.ug.edu.gh


#______________________________________________________
# List of packages to install -- to be run once
# After packages are installed

pkgs <- c(
    "MVTests", # Multivariate tests
    "mvhtests", # Multivariate tests
    "Hotelling", # Hotelling's T^2 Test and Variants
    "heplots", # Visualization for multivariate tests
    "energy", # Multivariate normal test
    "psych", # mardia test of Multivariate normality
    "dplyr", # data manipulation
    "codingMatrices", # for contrast matrices (Repeated measures)
    "readr", # import datasets
    "ggplot2", # graphics
    "gridExtra", # graphics tools
    "ICSNP" # Multivariate tests
)

install.packages(pkgs)

# Code in this block are to be run once. NO NEED TO REINSTALL THESE PACKAGES AGAIN

#______________________________________________________

# Clear global environment
rm(list = ls())

# For code reproducibility
set.seed(123)


# #########################
# II - MULTIVARIATE TESTS #
# #########################

##########
# IMPORTS
#########

# IMPORT REQUIRED PACKAGES

library(Hotelling)
library(dplyr) # data manipulation
library(MVTests)
library(mvhtests)
library(Hotelling)
library(heplots)
library(energy)
library(psych)
library(codingMatrices)
library(readr)
library(ggplot2)
library(gridExtra)
library(MASS) # in-built package to import
library(ICSNP)


#############################################
# --- 0. MULTIVARIATE NORMALITY (MVN) tests #
#############################################

# Section 7.7, P190, Applied Multivariate Statistics with R Second Edition

######
# DATA: Candy data from Applied Multivariate Statistics with R Second Edition Table 7.4 P193
######

candy_data <- data.frame(
    Name = c(
        "100 Grand", "3 Musketeers", "5th Avenue", "Almond Joy", "Andes Mints",
        "Baby Ruth", "Butterfinger", "Cadbury Dairy Milk", "Charleston Chew",
        "Dove Smooth Milk Choc.", "Goobers", "Heath Toffee", "Hershey's bar",
        "Hershey's Skor", "Junior Mints", "Kit Kat", "M&M's, peanut",
        "M&M's, plain", "Milk Duds", "Milky Way", "Mounds", "Mr Goodbar",
        "Nestle Crunch", "Oh Henry!", "Payday", "Raisinets", "Reese's Fast Break",
        "Reese's Nutrageous", "Reese's Peanut Butter cups", "Reese's Pieces",
        "Reese's Sticks", "Rolo", "Snickers", "Symphony", "Twix",
        "Whatchamacalit", "Whoppers", "Zero Candy Bar"
    ),
    Calories = c(
        190, 240, 260, 220, 200, 275, 270, 260, 230, 220, 200, 210, 210,
        200, 220, 207, 250, 230, 230, 240, 240, 250, 220, 280, 240, 190,
        260, 260, 210, 200, 220, 220, 230, 223, 250, 237, 190, 200
    ),
    Fat = c(
        8, 7, 12, 13, 13, 13, 11, 15, 6, 13, 13, 13, 13, 12, 4, 10, 13, 9, 8,
        9, 13, 16, 11, 17, 13, 8, 12, 16, 13, 9, 13, 10, 11, 13, 12, 11, 7, 7
    ),
    Satfat = c(
        5, 5, 5, 8, 11, 7, 6, 9, 5, 8, 5, 7, 8, 7, 3, 7, 5, 6, 5, 7, 10,
        7, 7, 7, 3, 5, 5, 5, 5, 8, 5, 7, 4, 8, 7, 8, 7, 5
    ),
    Carbs = c(
        30, 42, 38, 26, 22, 39, 43, 28, 43, 24, 20, 24, 26, 25, 45, 26, 30,
        34, 38, 37, 29, 25, 30, 36, 27, 32, 35, 28, 24, 25, 23, 33, 32, 24,
        33, 30, 31, 34
    ),
    Sugar = c(
        22, 36, 29, 20, 20, 32, 29, 28, 30, 22, 17, 23, 24, 24, 42, 20, 25,
        31, 27, 31, 21, 22, 24, 32, 21, 28, 30, 22, 21, 21, 17, 29, 27, 23,
        24, 23, 24, 29
    ),
    Sodium = c(
        90, 90, 120, 50, 20, 138, 135, 0, 30, 25, 15, 135, 35, 130, 35, 22,
        25, 35, 135, 75, 55, 50, 60, 65, 120, 15, 190, 100, 150, 55, 130, 80,
        115, 42, 100, 144, 100, 105
    )
)

# Visualize the data
View(candy_data)

# Remove the 'Name' column from the candy dataframe
# RUN ONLY ONE OF THE FOLLOWING

candy <- candy_data %>% dplyr::select(-Name) # %>% is a dplyr function
# OR similarly
candy <- candy_data[, -1] # need index of 'Name' column (Here it is 1)

#############################################
# --- 0.1 Using Chi-square / Gamma plots
#############################################

# Using Squared Generalized Distances - Mahalanobis distance or Chi-Square Distance (D_j ^2)

# ?mahalanobis # help on function
# Original Mahalanobis distances
mah_og <- mahalanobis(candy, colMeans(candy), var(candy))


# Sorting while keeping indices
sorted_idx <- order(mah_og) # returns the indices that would sort the squared distances in "mah_og"
mah <- mah_og[sorted_idx] # reorder the squared distances in "mah_og" by the index in "sorted_idx"
sorted_names <- candy_data$Name[sorted_idx] # Sort the names of the candies by the index in "sorted_idx"

rm(candy_data, mah_og, sorted_idx) # clear objects to keep global environment clean

# Why not use the function sort() ?
# sort() only gives sorted values — it does not give the indices of how the original elements were rearranged.

# There are p = 6 data columns (calories, fat, saturated fat, carbohydrates, sugar,
# and sodium), so the ordered Mahalanobis distances from the overall six-dimensional
# mean are compared to a chi-squared with 6 df. This is only an approximation because
# we are using estimates of the means and variances as their true values

################################################################
# --- QQ plot of Mahalanobis distances of candy nutrition values
# (vs Chi-squared quantiles with p=6 df)
################################################################

# Number of variables (degrees of freedom for chi-squared)
p <- ncol(candy)

# Theoretical quantiles (chi-squared)
chisq_quantiles <- qchisq(ppoints(length(mah)), df = p)

# Identify outliers (e.g., top 5% of χ²-distribution)
threshold <- qchisq(0.95, df = p) # 95th percentile cutoff
outliers <- mah > threshold

rm(threshold)

# Create QQ plot
qqplot(chisq_quantiles, mah,
       main = "QQ Plot of Mahalanobis Distances vs Chi-Squared",
       xlab = "Theoretical Quantiles (Chi-Squared)",
       ylab = "Mahalanobis Distances",
       col = "red", pch = 19
)
abline(0, 1, col = "black", lwd = 3) # Straight line of slope 1

# qqline(mah,
#     distribution = function(x) qchisq(x, df = p),
#     col = "green", lwd = 2
# )


# Add labels for outliers identified
text(
    x = chisq_quantiles[outliers],
    y = mah[outliers],
    labels = sorted_names[outliers], # Use sorted_names to identify outliers
    pos = 2, # Position labels to the left of points. see ?text
    cex = 0.9, # Label size
    col = "black" # Label color
)
# "Oh Henry!" and "Junior Mints" are outliers

# COMMENT:
# All values appear to fall on the QQ line connecting the upper and lower quantiles of the chi-squared distribution and those of the data.
# From this plot, we see there is no evidence of a lack of fit, and the data values appear to follow a multivariate normal distribution.

##########################
# --- 0.2 Correlation test
##########################

r_Q <- cor(mah, chisq_quantiles)
r_Q # 0.9936129
# The high correlation coefficient suggests that the data follows a multivariate normal distribution.

# _________________________________

# Remove some objects
rm(
    chisq_quantiles, mah,
    outliers, p,
    sorted_names
)
# _________________________________

#####################
# --- 0.3 Energy test
#####################
# Applied Multivariate Statistics with R Second Edition, P196

# The energy test is programmed in the library of the same name

# Linux system note: install GSL (GNU Scientific Library) system package
# install.packages("energy")
# library(energy)
energy::mvnorm.etest(candy, R = 999)

# p-value = 0.8238
# The energy test does not indicate a lack of fit to the multivariate normal distribution

# _________________________________

#####################
# --- 0.4 Mardia test
#####################
# Applied Multivariate Statistics with R Second Edition, P195

# library(psych)
psych::mardia(candy)
# No evidence of extreme multivariate skewness or kurtosis in these data

# ______________________________________________________________________________

#######################
# --- 1. BOX'S M Test #
#######################

# SAMPLES FROM MVN

# MVTests::BoxM
# mvhtests::Mtest.cov
# heplots::boxM

data(iris) # {datasets} or {MVTests}
# Variables: Sepal.Length, Sepal.Width, Petal.Length, Petal.Width
# Groups: 3 levels of Species ( Setosa, Versicolor. Virginica)

# --- {MVTests}
MVTests::BoxM(data = iris[, 1:4], group = iris[, "Species"])

# K = 140.943, df = 20, p-value < 0.05
# INTERPRETATION: We reject the null hypothesis.
# We conclude that we have at least two unequal covariance matrices.

# --- {mvhtests}
mvhtests::Mtest.cov(x = as.matrix(iris[, 1:4]), ina = iris[, "Species"], a = 0.05)
# ...

# --- {heplots}
heplots::boxM(iris[, 1:4], iris[, "Species"])
# ...

# ______________________________________________________________________________

####################################
# --- 2. ONE SAMPLE Hotelling's T2 #
####################################

# DATA ~ MULTIVARIATE NORMAL

# DATA:
# TABLE 5.1 -- SWEET DATA from "Applied Multivariate Statistical Analysis", 5th Ed, P215
sweet_data <- data.frame(
    Sweat.Rate = c(
        3.7, 5.7, 3.8, 3.2, 3.1, 4.6, 2.4, 7.2, 6.7, 5.4,
        3.9, 4.5, 3.5, 4.5, 1.5, 8.5, 4.5, 6.5, 4.1, 5.5
    ),
    Sodium = c(
        48.5, 65.1, 47.2, 53.2, 55.5, 36.1, 24.8, 33.1, 47.4, 54.1,
        36.9, 58.8, 27.8, 40.2, 13.5, 56.4, 71.6, 52.8, 44.1, 40.9
    ),
    Potassium = c(
        9.3, 8.0, 10.9, 12.0, 9.7, 7.9, 14.0, 7.6, 8.5, 11.3,
        12.7, 12.3, 9.8, 8.4, 10.1, 7.1, 8.2, 10.9, 11.2, 9.4
    )
)

# Normality check
energy::mvnorm.test(sweet_data, R = 199)
# ...


# TEST:
# We have assumed that the sweat data are multivariate normal
# H0: mu' = [4, 50, 10] at alpha = .10
mu0 <- c(4, 50, 10)
MVTests::OneSampleHT2(sweet_data, mu0 = mu0, alpha = 0.10)

# T2 = 9.738 ~ 9.74, dfs = [3, 17], p-value = 0.06492834 < 0.1
# INTERPRETATION: We reject H0 at the 10% level of significance

# No info on T2
mvhtests::hotel1T2(as.matrix(sweet_data), M = mu0, a = 0.10)
# ...

ICSNP::HotellingsT2(sweet_data, mu = mu0)
# p-value consistent with MVTests
# T2 = 2.9045 != 9.738 -- T2 specification

# _________________________________

rm(mu0) # remove object mu0

# ______________________________________________________________________________

####################################
# --- 3. TWO-SAMPLE Hotelling's T2 #
####################################

# 2 INDEPENDENT SAMPLES
# MVN
# EQUAL COVARIANCE MATRICES

# DATA:
# TABLE 6.9 -- CARAPACE MEASUREMENTS (in mm) for PAINTED TURTLES
# from "Applied Multivariate Statistical Analysis", 5th Ed, P339
# Hint from exercise: Logarithmic transformations of the observations

turtles_female <- data.frame(
    Length = c(
        98, 103, 103, 105, 109, 123, 123, 133, 133, 133, 134, 136,
        138, 138, 141, 147, 149, 153, 155, 155, 158, 159, 162, 177
    ),
    Width = c(
        81, 84, 86, 86, 88, 92, 95, 99, 102, 102, 100, 102,
        98, 99, 105, 108, 107, 107, 115, 117, 115, 118, 124, 132
    ),
    Height = c(
        38, 38, 42, 42, 44, 50, 46, 51, 51, 51, 48, 49,
        51, 51, 53, 57, 55, 56, 63, 60, 62, 63, 61, 67
    )
)

turtles_male <- data.frame(
    Length = c(
        93, 94, 96, 101, 102, 103, 104, 106, 107, 112, 113, 114,
        116, 117, 117, 119, 120, 120, 121, 125, 127, 128, 131, 135
    ),
    Width = c(
        74, 78, 80, 84, 85, 81, 83, 83, 82, 89, 88, 86,
        90, 90, 91, 93, 89, 93, 95, 93, 96, 95, 95, 106
    ),
    Height = c(
        37, 35, 35, 39, 38, 37, 39, 39, 38, 40, 40, 40,
        43, 41, 41, 41, 40, 44, 42, 45, 45, 45, 46, 47
    )
)

n_turtles <- length(turtles_female$Length)
sex <- rep(1:2, each = n_turtles)
turtles <- rbind(turtles_female, turtles_male) # putting the female and male datasets together
# first comes the rows in 'turtles_female' then the rows in 'turtles_male'



# MVN check
energy::mvnorm.test(turtles_female, R = 199)
energy::mvnorm.test(turtles_male, R = 199)

# Box'M test
MVTests::BoxM(data = turtles, group = sex)
# ...

# TEST:

# --- {MVTests}
MVTests::TwoSamplesHT2(
    data = turtles,
    group = sex,
    alpha = 0.05,
    Homogenity = TRUE # Equal Covariance Matrix
)
#

# --- {mvhtests}
# Equal Covariance Matrix
mvhtests::hotel2T2(
    as.matrix(turtles_female),
    as.matrix(turtles_male),
    a = 0.05
)
# Without Assumption of Eq. Cov. Matrix
mvhtests::james(
    as.matrix(turtles_female),
    as.matrix(turtles_male),
    a = 0.05
)
# ...

# --- {Hotelling}
two_sample_test <- Hotelling::hotelling.test(
    as.matrix(turtles_female),
    as.matrix(turtles_male),
    var.equal = TRUE, # Equal Covariance Matrix
)
two_sample_test

# ...

# --- {ICSNP}
# equal covariance matrix
ICSNP::HotellingsT2(
    turtles_female,
    turtles_male,
    mu = NULL,
    test = "f"
)
# unequal covariance matrix
ICSNP::HotellingsT2(
    turtles_female,
    turtles_male,
    mu = NULL,
    test = "chi"
)
# ...

rm(sex, n_turtles)

# ______________________________________________________________________________

#######################################
# --- 4. PAIRED-SAMPLE Hotelling's T2 #
#######################################

# PAIRED DIFFERENCES ~ MVN

# library(ICSNP)

# DATA1: ICSNP {ICSNP}
# ?LASERI
data(LASERI)

# The LASERI data in the ICSNP library is a data frame with 32 measurements made on each of 223 healthy Finnish subjects.
# The subjects were monitored while in a supine position and then again with their heads elevated on a motorized table.
# We will concentrate on four measurements and their average differences: average heart rate (HRT1T4); average cardiac output (COT1T4); average systemic vascular resistance index (SVRIT1T4); and average pulse wave velocity (PWVT1T4).
# Each of these variables is expressed as a difference of the pre- and post-tilt values. Testing equality of means of pre- and post- values is the same as testing whether each of the variables in the data has a zero population mean

# --- Using ICSNP::HotellingsT2 -- USES THE DIFFERENCES of paired differences to calculate T2
# LASERI data already formatted as differences
ICSNP::HotellingsT2(LASERI[, 25:28])

# p < 2.2e-16
# Considerable evidence that the mean measurements are very different in the
# supine and tilted positions. We conclude tilting the subject results in very different
# measures.

# DATA2: TABLE 6.1 -- EFFLUENT DATA from "Applied Multivariate Statistical Analysis", 5th Ed, P215
# Also Example P11 - Slides LA, Multivariate Methods, Lecture Notes

LABCommercial <- data.frame(
    BOD = c(
        6, 6, 18, 8, 11, 34, 28, 71, 43, 33, 20
    ),
    SS = c(
        27, 23, 64, 44, 30, 75, 26, 124, 54, 30, 14
    )
)

LABState <- data.frame(
    BOD = c(
        25, 28, 26, 35, 15, 44, 42, 54, 34, 29, 39
    ),
    SS = c(
        15, 13, 22, 29, 31, 64, 30, 64, 56, 20, 21
    )
)

Effluent_diff <- CommercialLAB - StateLAB

# MVN check
energy::mvnorm.etest(Effluent_diff, R = 199)

# TEST:

# --- Using original data
# Raw result
MVTests::Mpaired(T1 = CommercialLAB, T2 = StateLAB)
# T2 = 14.90961 != 13.6 in book -- specification of T2 in package
# df=[2, 9], F=6.709324
# Critical value: matter of scaling
# p-value = 0.01645694 < 0.05
# INTERPRETATION
# ...

# More concise result
summary(MVTests::Mpaired(T1 = CommercialLAB, T2 = StateLAB))

# --- Using Differences of paired observations
ICSNP::HotellingsT2(Effluent_diff)
# T2 = 6.7093 != 13.6 -- ...
# ...

# No info on T2
mvhtests::hotel1T2(as.matrix(Effluent_diff), M = c(0, 0), a = 0.05)

# ______________________________________________________________________________


###################################
# --- 5. REPEATED MEASURES DESIGN #
###################################

# ASSUMPTION: MVN

# DATA:
# TABLE 6.2 -- SLEEPING-DOG DATA from "Applied Multivariate Statistical Analysis", 5th Ed, P281

# Treatment 1 = high C02 pressure without H
# Treatment 2 = low C02 pressure without H
# Treatment 3 = high C02 pressure with H
# Treatment 4 = low C02 pressure with H

sleep_dog <- data.frame(
    Treatment1 = c(
        426, 253, 359, 432, 405, 324, 310, 326, 375, 286,
        349, 429, 348, 412, 347, 434, 364, 420, 397
    ),
    Treatment2 = c(
        609, 236, 433, 431, 426, 438, 312, 326, 447, 286,
        382, 410, 377, 473, 326, 458, 367, 395, 556
    ),
    Treatment3 = c(
        556, 392, 349, 522, 513, 507, 410, 350, 547, 403,
        473, 488, 447, 472, 455, 637, 432, 508, 645
    ),
    Treatment4 = c(
        600, 395, 357, 600, 513, 539, 456, 504, 548, 422,
        497, 547, 514, 446, 468, 524, 469, 531, 625
    )
)

# TEST:
mvhtests::rm.hotel(as.matrix(sleep_dog), a = 0.05)
# ...

# User-defined function for test
# library(codingMatrices)
RepeatedMeasures.test <- function(data, alpha) {
    p <- ncol(data) # number of variables
    # https://cran.r-project.org/web/packages/codingMatrices/vignettes/codingMatrices.pdf#page=6.36
    C <- codingMatrices::mean_contrasts(MASS::contr.sdif(p)) # Contrast matrix (Successive)
    n <- length(data[, 1]) # size of the data
    xbar <- colMeans(data) # mean vector
    S <- var(data)
    T2 <- n * t(C %*% xbar) %*% solve((C %*% S %*% t(C))) %*% (C %*% xbar) # T2 statistic
    dfs <- c(
        numdf = p - 1,
        denomdf = n - p + 1
    )
    CriticalValue <- ((n - 1) * dfs[1] / dfs[2]) * qf(1 - alpha, df1 = dfs[1], df2 = dfs[2]) # Critical Value
    p.value <- pf((dfs[2] * T2) / ((n - 1) * dfs[1]), df1 = dfs[1], df2 = dfs[2], lower.tail = FALSE)
    return(
        list(
            num.vars = p, size = n, degrees.freedom = dfs, true.mean = xbar,
            dispersion = S, T.squ.stat = T2, CriticalValue = CriticalValue
            p.value = p.value
        )
    )
}

RepeatedMeasures.test(sleep_dog, alpha = 0.05)

# ______________________________________________________________________________

#########################
# --- 6. ONE-WAY MANOVA #
#########################

# https://www.r-bloggers.com/2021/11/manovamultivariate-analysis-of-variance-using-r/

# ASSUMPTIONS
# MVN
# EQUAL COVARIANCE MATRICES

# MANOVA Hypotheses
# Null hypothesis: group mean vectors are same for all groups
# Alternative hypothesis: group mean vectors are not same for all groups

# DATA2: Exercise P40, LA, Multivariate Methods, Lecture Notes
# treatment_manova <- data.frame(
#     treatment = rbind(rep(1, 3), rep(2, 4))
# )

# DATA1:
# dataset of various plant varieties (plant_var) and their associated phenotypic
# measurements for plant heights (height) and canopy volume (canopy_vol).
# We want to see if plant heights and canopy volume are
# associated with different plant varieties



# library(readr)
df <- readr::read_csv("https://reneshbedre.github.io/assets/posts/ancova/manova_data.csv")

# library(gridExtra)

p1 <- ggplot2::ggplot(df, aes(x = plant_var, y = height, fill = plant_var)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    theme(legend.position = "top")

p2 <- ggplot2::ggplot(df, aes(x = plant_var, y = canopy_vol, fill = plant_var)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    theme(legend.position = "top")

gridExtra::grid.arrange(p1, p2, ncol = 2)

# Create matrix of dependent vars
dep_vars <- cbind(df$height, df$canopy_vol)

# _________________________________
# Box's M test
MVTests::BoxM(dep_vars, group = df$plant_var)


# _________________________________
# Perform one-way MANOVA
fit <- stats::manova(dep_vars ~ plant_var, data = df)
summary(fit)

# p < 0.001
# plant varieties have a statistically significant association with
# both combined plant height and canopy volume.


