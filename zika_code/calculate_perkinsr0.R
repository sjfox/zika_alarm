## Want a function that 
# 1 Takes in a list of months
# 2 Calculates all rnotts for those months
# 3 Saves to csv 

#############################################################
# Step 1 load Perkins Nature Microbiology data and Texas data 
#############################################################
source('r0functions.R')
load("../data/perkins_sims/functions_R0_AR_random_draws.RData")
load("~/Downloads/parms_fxns_r0.RData")
county.parms <- read.csv("../csvs/county_r0_parameters.csv")
tx.temps <- read.csv("../csvs/tx_county_temps.csv")


#####################################################################
# Functions for calculating rnot according to a specified temperature 
#####################################################################

select_month <- function(temps, desired.months) {
  # Down select the temperatures 
  temps %>% dplyr::select(county, which(colnames(temps) %in% 
                    eval(parse(text = "desired.months"))))
}

align_data <- function(temps.selected, county.parms) {
  full.data <- left_join(county.parms, temps.selected, by = "county")
}


############################################################## 
## This first function will calculate full R0 distribution range for each of the specified months
## Will write a csv of the distribution, with a 254x1000 matrix for each month 

desired.months = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")

csv.path = "../csvs/rnot_temp/"

calculate_rnots_bymonth(temps = tx.temps, desired.months = desired.months, county.parms = county.parms)

##############################################################
# Returns 254x1000 dataframe 
# Only works for one month-have to specify the month with temperature  
# Can use the output of this directly in the next step (R0Script_Clean)
rnot.ests <- rnot_calc_cty_dist(county.parms$mosquito.abundance, 
                                county.parms$gdp, 
                                temperature = tx.temps$Aug, 
                                a=a, b=b, c.r=c.r, 
                                mort.fun.list=mort.fun, 
                                eip.fun.list=eip.fun, 
                                scam.est.list=scam.est.list)

############################
