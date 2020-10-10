source("MarkovLinearModelExamples.R")

##############################
#
# This file contains an implementation of the examples given in
# "Estimations of means and variances in a Markov linear model" 
# by A. Gutierrez and S. Mueller
#
##############################


##############################
#### Example 11
##############################

set.seed(42)
number_col <- 2   # number of columns
number_lines <- numeric(number_col)
number_lines <- c(2, 2)           # number of machines in each column


### Definition of the transition kernels
Q <- list(number_col)

Q[[1]] <- c(1/2, 1/2)
Q[[2]] <- matrix(c(3/4, 1/4, 1/4, 3/4), nrow=2, byrow=TRUE)

### Definition of the quality matrix
V <- matrix( c(2,1,1,1), nrow=2, byrow = TRUE)      # matrix of variances
ES <- matrix( c(0,1,-2,2), nrow=2, byrow = TRUE)    # matrix of means

### Sampling
observation <- sample_path_quality(number_col, number_lines, Q, ES,V)
n <- 1000
for (i in 2:n){
  observation<-bind_rows(observation, sample_path_quality(number_col, number_lines, Q, ES,V))
}


Q_estimator(observation, number_col, number_lines) -> Q_est
T_U(observation, number_col, number_lines) -> TU
round(apply(TU, 2, function(x) x[1]-x[-1]),2)

T_QU(observation, number_col, number_lines, Q) -> TQU
round(apply(TQU, 2, function(x) x[1]-x[-1]),2)

T_QU(observation, number_col, number_lines, Q_est) -> TQU_est
round(apply(TQU_est, 2, function(x) x[1]-x[-1]),2)

S_U(observation, number_col, number_lines)->SU
round(apply(SU, 2, function(x) x[1]-x[-1]),2)

S_QU(observation, number_col, number_lines, Q)->SQU
round(apply(SQU, 2, function(x) x[1]-x[-1]),2)

S_QU(observation, number_col, number_lines, Q_est)->SQUest
round(apply(SQUest, 2, function(x) x[1]-x[-1]),2)


##############################
#### Example 13  - Tooth growth
##############################

data(ToothGrowth)
str(ToothGrowth)

# covariates <- supp and dose as factor
# response <- length

ToothGrowth -> data
data <- data[c(2, 3, 1)]
data$supp <-
  factor(data$supp,
         levels = c("VC", "OJ"),
         labels = c("VC", "OJ"))
data$dose <-
  factor(
    data$dose,
    levels = c(0.5, 1.0, 2.0),
    labels = c("Dose 0.5", "Dose 1.0", "Dose 2.0")
  )

str(data)

qplot(
  supp,
  len,
  data = data,
  facets =  ~ dose,
  main = "Tooth growth of guinea pigs",
  xlab = "Supplement type",
  ylab = "Tooth length"
) +
  geom_violin(aes(fill = supp)) +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  theme(legend.position = "none") 

### Classical linear regression with factors

fit <- aov(len ~ supp + dose, data = data)
summary(fit)
fit$coefficients

### Markov linear regression

number_col <- 2   # number of columns
number_lines <- numeric(number_col)  # number of lines in each column
number_lines <- c(2,3)   # first colum supp, second dose


n <- nrow(data)
observation <- data
observation$supp <- as.numeric(observation$supp)
observation$dose <- as.numeric(observation$dose)
names(observation)[number_col+1] <- "b"
T_U(observation, number_col, number_lines)-> TU

apply(TU, 2, function(x) x[-1]-x[1])

fit$coefficients
# we obtain the same results as with the linear regression

# we can check if the experiment is balanced

Q_estimator(observation, number_col, number_lines)



T_QU(observation, number_col, number_lines, Q_estimator(observation, number_col, number_lines)) -> TQU
TQU -TU



# We obtain the same results, since balanced experiment and no bias.

######## Add variance
S_U(observation, number_col, number_lines)-> SU
SU
apply(SU, 2, function(x) x[1]-x[-1])

S_QU(observation, number_col, number_lines, Q_estimator(observation, number_col, number_lines))-> SQU
SQU
apply(SU, 2, function(x) x[1]-x[-1])




data %>% group_by(supp) %>%
  summarize(variance =var(len)) -> res
(res$variance[1]-res$variance[2]) *29/30

data %>% group_by(dose) %>%
  summarize(variance =var(len)) -> res
(res$variance[1]-res$variance[2]) * 19 /20
(res$variance[1]-res$variance[3]) * 19 /20


##############################
#### Example 14  - Biomass response
##############################


gdn <-
  read.csv("http://www.math.montana.edu/courses/s217/documents/gundalebachnordin_2.csv")
gdn$Species <- factor(gdn$Species)
gdn$Treatment <- factor(gdn$Treatment)

qplot(
  Species,
  Massperha,
  data = gdn,
  facets =  ~ Treatment,
  main = "Influence of Nitrogen on growth ",
  xlab = "Species",
  ylab = "Mass per ha"
) +
  geom_violin(aes(fill = Species)) +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  theme(legend.position = "none")


### Look at classical linear regression with factors

fit <- aov(Massperha ~ Treatment + Species, data = gdn)
summary(fit)
fit$coefficients

fitlog <- aov(log(Massperha) ~ Treatment + Species, data = gdn)
summary(fitlog)

### Markov linear regression approach

number_col <- 2   # number of columns, her set to 4
number_lines <- numeric(number_col)
number_lines <-
  c(2, 3)   # first colum contains supp and second the factor  dose


data <- data.frame(gdn[c(1, 2, 4)])
n <- nrow(data)
observation <- data
observation$Species <- as.numeric(observation$Species)
observation$Treatment <- as.numeric(observation$Treatment)
names(observation)[number_col + 1] <- "b"
T_U(observation, number_col, number_lines) -> TU

apply(TU, 2, function(x)
  x[-1] - x[1])

fit$coefficients

# we can check if the experiment is balanced

Q_estimator(observation, number_col, number_lines)

T_QU(
  observation,
  number_col,
  number_lines,
  Q_estimator(observation, number_col, number_lines)
) -> TQU
TQU - TU



# We obtain the sames result; experiment is balanced and no bias.

######## Estimation of the variance

S_U(observation, number_col, number_lines) -> SU
SU
apply(SU, 2, function(x)
  x[1] - x[-1])

S_QU(
  observation,
  number_col,
  number_lines,
  Q_estimator(observation, number_col, number_lines)
) -> SQU
SQU
apply(SU, 2, function(x)
  x[1] - x[-1])

##############################
#### Example 15  - CASchools
##############################
library("AER")
library("MASS")
data(CASchools)

CASchools$STR <- CASchools$students/CASchools$teachers 
CASchools$score <- (CASchools$read + CASchools$math)/2

cor(CASchools$STR, CASchools$score)
cor(CASchools$STR, CASchools$english)

### Plot the CASchool data
ggplot(data = CASchools, aes(x = STR, y = score)) +
  geom_point(alpha = .8, aes(size = english)) +
  labs(title = "The California test score data")


### Create categorical predictor variables
CASchools$STRCat <-
  cut_number(CASchools$STR,
             5,
             labels = 1:5)
CASchools$englishCat <-
  cut_number(
    CASchools$english,
    5,
    labels = c("English 1", "English 2", "English 3", "English 4", "English 5")
  )

ggplot(data = CASchools, aes(STRCat, score)) +
  geom_violin(aes(fill = STRCat)) +
  geom_jitter(shape = 16,
              position = position_jitter(0.1),
              size = 0.5) +
  labs(title = "Influence on the score") +
  facet_wrap(~ englishCat, ncol = 5)

#### multivariate linear regression
mult.modcat <- lm(score ~ STRCat + englishCat, data = CASchools)
summary(mult.modcat)

### Markov linear regression approach
CASchools$englishCat <-
  cut_number(CASchools$english,
             5,
             labels = 1:5)
data <- (CASchools[c(18,17,16)])
number_col <- 2

number_lines <- c(5, 5)
data[1, ] <- as.numeric(data[1, ])
data[2, ] <- as.numeric(data[2, ])
names(data)[number_col + 1] <- "b"
summary(data)
T_U(data, number_col, number_lines) -> TU
apply(TU, 2, function(x)
  x[-1] - x[1])

Q_estimator(data, number_col, number_lines) -> Q_est
Q_est

T_QU(data, number_col, number_lines, Q_est) -> TQU
apply(TQU, 2, function(x)
  x[-1] - x[1])

######## Estimation of the variance

S_U(data, number_col, number_lines) -> SU
SU
apply(SU, 2, function(x)
  x[1] - x[-1])

S_QU(data,
     number_col,
     number_lines,
     Q_est) -> SQU
SQU
apply(SQU, 2, function(x)
  x[1] - x[-1])

## Variances for each predictor
data %>% 
    group_by(englishCat) %>%
    summarise(variance = sd(b)^2)

data %>% 
  group_by(STRCat) %>%
  summarise(variance = sd(b)^2)
