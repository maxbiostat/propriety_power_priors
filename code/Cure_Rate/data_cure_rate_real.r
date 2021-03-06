## Historical data, E1673 trial
historical.data <- read.csv("../data/e1673.jasa.dat", sep = " ")
stand_age_0 <- (historical.data$age -mean(historical.data$age))/sd(historical.data$age)
N_0 <- length(stand_age_0)
X_0 <- cbind(rep(1, N_0), # x0 - dummy for intercept
           stand_age_0, # x1 - (standardised) age 
           historical.data$gender, # x2 - gender (sex)
           historical.data$PS # x3 performance status
)
Y_0_cens <- historical.data$stime
Delta_0 <- historical.data$censor


## New data, E1684 trial
full.data <- read.csv("../data/e1684_and_e1690_data.csv", sep = "\t")
new.data <- subset(full.data,
                          study  == "1684" & survtime > 0)

stand_age <- ( new.data$age - mean(new.data$age))/sd(new.data$age)
N <- length(stand_age)
X <- cbind(rep(1, N), # x0 - dummy for intercept
             stand_age, # x1 - (standardised) age 
             new.data$sex, # x2 - gender (sex)
             new.data$perform # x3 performance status
)
Y_cens <- new.data$survtime
Delta <- new.data$scens

## Hyperparameters

mu0 <- 0
s0 <- 100

d0 <- 1
tau0 <- 0.01

## \pi_A parameters
eta <- 100
nu <- 100
