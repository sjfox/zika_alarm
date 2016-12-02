## Generate R0s
library(tidyverse)
library(stringr)


temps <- read_tsv("../data/weather/tx_county_temps.txt")
temps <- temps %>% filter(is.na(Notes)) %>%
  rename(avg_max_temp = `Avg Daily Max Air Temperature (C)`,
         avg_min_temp=`Avg Daily Min Air Temperature (C)`) %>%
  mutate(avg_temp = (avg_max_temp + avg_min_temp) / 2,
         county = tolower(str_replace_all(County, pattern = " County, TX", ""))) %>%
  arrange(`Month Code`) %>%
  mutate(Month = factor(Month, levels = unique(Month))) %>%
  dplyr::select(County, month=Month, avg_temp) %>%
  spread(key = month, value = avg_temp) %>% 
  separate(col = County, into = c("Area.Name", "County"), sep = -12) %>% 
  dplyr::select(-County)



tx_county <- read_csv(file = "data/texas_county_info.csv")

tx_county <- tx_county %>% select(county, gdp, mosquito_abundance, ex_inc_period) %>%
  mutate(county=tolower(str_replace_all(county, pattern = " County, Texas", "")),
         mosq_expected = exp(-1.79 - (.07 / .5) * (log(gdp)-10)) * mosquito_abundance) %>%
  left_join(temps,by=c("county"))

eip <- function(temp, tau=4.9, b0=2.9, bt=-0.08){
  mu <- exp(b0 + bt * temp)
  exp(mu + tau^(-1)/2)
  # exp(-tau * (log(t) - mu)^2 / 2 ) * (1 / t) * sqrt(tau / 2 / pi)
}

calc_r0 <- function(eip, exp_mosq, c_over_r=9*.77, mosq_hum_trans=.634, alpha=.63, mosq_mort=1/14) {
  c_over_r *(exp_mosq * mosq_hum_trans * alpha^2 * exp(-mosq_mort * eip)) / mosq_mort
}


rnot_dat <- tx_county %>% mutate_at(vars(Jan:Dec), funs(round(calc_r0(eip(.), mosq_expected), digits=5)))

devtools::use_data(rnot_dat)



