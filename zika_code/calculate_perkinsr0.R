## Want a function that 
# 1 Takes in a list of months
# 2 Calculates all rnotts for those months
# 3 Saves to csv 

#############################################################
# Step 1 load Perkins Nature Microbiology data and Texas data 
#############################################################
source('r0functions.R')
load("../data/perkins_sims/functions_R0_AR_random_draws.RData")
county_parms <- read.csv("../csvs/county_r0_parameters.csv")
tx_temps <- read.csv("../csvs/tx_county_temps.csv")


#####################################################################
# Functions for calculating rnot according to a specified temperature 
#####################################################################

select_month <- function(temps, desired.months) {
  # Down select the temperatures 
  temps %>% dplyr::select(county, which(colnames(temps) %in% 
                    eval(parse(text = "desired.months"))))
}

align_data <- function(temps.selected, county_parms) {
  full.data <- left_join(county_parms, temps.selected, by = "county")
}


############################################################ 
desired.months = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
desired.month.decrease = "Aug"

csv.path = "../csvs/rnot_temp/"

calculate_rnots_bymonth(temps = tx_temps, desired.months = desired.months, county_parms = county_parms)

# Original perksin methods-returns median/quantile of rnott estimates
# have to specify temperature as a month 
rnot_ests <- rnot_calc_counties(county_parms$mosquito.abundance, 
                                county_parms$gdp, 
                                temperature, 
                                a=a, b=b, c.r=c.r, 
                                mort.fun.list=mort.fun, 
                                eip.fun.list=eip.fun, 
                                scam.est.list=scam.est.list)


############################
econ <- seq(9, 12, by=0.1)
which(econ > 9)
which(econ > 12)
library(scam)
county_parms <- read.csv("../csvs/county_r0_parameters.csv")
range(log(county_parms$gdp)): r# range 9-12, econ 32-61


# Function that determines which scam functions have a negative relationship 
index = seq(1:1000)

get_negative_slopes = function(gdp.min, gdp.max, scam.est.list, index) {
  answer <- lapply(X = index, FUN = function(x) {
    prediction = round(predict(scam.est.list[[x]], newdata = data.frame(econ = econ)), digits = 1)
    
    #is it flat? 
    if (prediction[1] - prediction[31] == 0) {
      answer[[x]] = "flat"
    }
    else if ( which.min(c(prediction[1], prediction[31])) == 2) {
      answer[[x]] = "decrease"
    }
    else {
      answer[[x]] = "increase"
    }
  })
}

flat.functions <- scam.est.list[which(answer == "flat")]
plot(econ, predict(flat.functions[[1]], newdata=data.frame(econ=econ)), type="l", ylim=c(-2, 5))
for(ii in 2:239){
  lines(econ, predict(flat.functions[[ii]], newdata=data.frame(econ=econ)), type="l")
}

decrease.functions <- scam.est.list[which(answer == "decrease")]
plot(econ, predict(decrease.functions[[1]], newdata=data.frame(econ=econ)), type="l", ylim=c(-2, 5))

decrease.index = which(answer == "decrease")
for(ii in 2:604){
  lines(econ, predict(decrease.functions[[ii]], newdata=data.frame(econ=econ)), type="l")
}
