library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(deSolve)
library(optimr)

population <- 146748590 * 0.3 #0.3
lock_start <- 89 #23 march 90
lock_end <- 167 #9 june 170



get_data_russia <- function() {
  data_confirmed <- read.csv("https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
  data_recovered <- read.csv("https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv")
  data_deaths <- read.csv("https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
  data_confirmed <- pivot_longer(data_confirmed, cols = seq(5, ncol(data_confirmed), 1), names_to = "Date")
  data_recovered <- pivot_longer(data_recovered, cols = seq(5, ncol(data_recovered), 1), names_to = "Date")
  data_deaths <- pivot_longer(data_deaths, cols = seq(5, ncol(data_deaths), 1), names_to = "Date")
  data_confirmed <- data_confirmed %>% filter(Country.Region == "Russia")
  data_recovered <- data_recovered %>% filter(Country.Region == "Russia")
  data_deaths <- data_deaths %>% filter(Country.Region == "Russia")
  data_confirmed <- data_confirmed[10:nrow(data_confirmed), ]
  data_recovered <- data_recovered[10:nrow(data_recovered), ]
  data_deaths <- data_deaths[10: nrow(data_deaths), ]
  data_russia <- tibble(I = data_confirmed$value - data_recovered$value - data_deaths$value, 
                        R = data_recovered$value, 
                        S = population - data_confirmed$value, 
                        D = data_deaths$value,
                        cumsum = data_confirmed$value,
                        total = data_confirmed$value,
                        day = seq(1, nrow(data_confirmed)))
  gather(data_russia, type, value, I:total, factor_key=TRUE)
}

plot_data <- function(data, t) {
  ggplot(data = data %>% filter(type == t), aes(x = day, y = value, col = type)) +
  geom_point(size = 0.1)
}

plot_model <- function(data, model, t) {
  ggplot(data = data %>% filter(type == t), aes(x = day, y = value, col = type)) +
    geom_line(size = 0.2) +
    geom_line(data = model %>% filter(type == t), aes(x = day, y = value, col = type),
              linetype = 'dotted', size = 1) + xlab("day") +
          ylab("value") + geom_vline(xintercept = 380, linetype = 'dotted', size = 0.5, alpha = 0.5) +
    ggtheme)
    
}

theme_set(theme_light())

plot_cs <- function(data, model) {
  ggplot(data = data, aes(x = day, y = value)) + geom_line(size = 0.2) + 
    geom_line(data = model, aes(x = day, y = value), linetype = 'dotted', 
              size = 1)
}

data_russia <- get_data_russia()
plot_data(data_russia)

lfn2 <- function(p) {
  out <- SEIRD(S0, E0, I0, R0, D0, alpha1 = p[1]  ,
               alpha2 = p[2]  ,
               alpha3 = p[3]  ,
               beta1 = p[4],
               gamma1 = p[5],
               eta1 = p[6]  ,
               n_days = 390,             #388
               zeta1 = p[7]   , 
               zeta2 = p[8]    , 
               lock_start = p[9], 
               lock_end = p[10],
               delta1 = p[11] )
  
  rss <- 
         (sum(((data_russia %>% filter(type == 'I'))$value -
           (out %>% filter(type == 'I'))$value) ^ 2)) / (population * population)
  print(rss)
  return(rss)
}


SEIRD <- function(S0, E0, I0, R0, D0, alpha1, alpha2, alpha3, beta1, gamma1,
                  eta1, n_days, zeta1, zeta2, lock_start, lock_end, delta1) {
  out <- tibble(day = 1, S1 = S0 + E0, S = S0, E = E0, I = I0, R = R0, D = D0, 
                cumsum = I0)
  
  total <- I0
  alpha = alpha1
  for (i in 2:n_days) {
    S0i <- S0
    E0i <- E0
    I0i <- I0
    R0i <- R0
    D0i <- D0

    if (i < lock_start) alpha = alpha1
    if (i >= lock_start && i < lock_end) alpha = max(alpha2, alpha - zeta1)
    if (i >= lock_end) {
      alpha = min(alpha3, alpha + zeta2)
      
    } 
    if (i >= 430) delta1 = 0

    S0 <- max(0, S0i - alpha * S0i * I0i * 1 / population + delta1 * R0i)
    E0 <- max(0, E0i + alpha * I0i * S0i * 1 / population - gamma1 * E0i)
    I0 <- max(0, I0i + gamma1 * E0i - beta1 * I0i - eta1 * I0i)
    R0 <- max(0, R0i + beta1 * I0i - delta1 * R0i)
    D0 <- max(0, D0i + eta1 * I0i)
    total <- I0 + D0 + R0
    out <- out %>% add_row(day = i, S1 = S0, S = S0 + E0, E = E0, I = I0, R = R0,
                           D = D0, cumsum = total)
  }
  return(gather(out, type, value, S1:cumsum, factor_key=TRUE))
}

S0 <- population
E0 <- 0
I0 <- 2
R0 <- 0
D0 <- 0

out <- SEIRD(S0, E0, I0, R0, D0, alpha1 = 0.7673478   ,
             alpha2 = 0.4152731   ,
             alpha3 = 0.9476324   ,
             beta1 = 0.370029, 
             gamma1 = 0.27,
             eta1 = 0.0007751027   ,
             n_days = 500,        
             zeta1 = 0.03674178    , 
             zeta2 = 0.001360923     , 
             lock_start = 89, 
             lock_end = 189,
             delta1 = 0.004786315 )

p= c(0.7673478,
     0.4152731,
     0.9476324,
     0.370029,
     0.27,
     0.0007751027,
     0.03674178,
     0.001360923,
     89,
     189,
     0.004786315)

upper = c(2, 2, 2, 1, 1, 1, 1, 1, 100, 200, 0.01)
lower = c(0, 0, 0, 0, 0, 0, 0, 0, 50, 150, 0)

S_ = out %>% filter(type == 'S') + out %>% filter(type == 'E')
plot_model(data = data_russia, model = out, c('I',  'E', 'D'))
plot_model(data = data_russia, model = out, c('R'))

ggplot(data_russia, aes(x = day, y = cs_russia)) + geom_point()



meth0 <- c("Nelder-Mead")
res1 <- opm(p,
            lfn2, "grfwd", method=meth0)
res1
