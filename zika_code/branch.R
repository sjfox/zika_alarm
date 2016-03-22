
# prop_p -- probability an I infects a new individual in a time period: how to determine this? 
# recov_p -- probability an I recovers in a time period: human recovery 
# disc_p -- probability we discover an I: symptom ratio 
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative):

#Parameters 
prop_p <- 1.7/7
recov_p <- 1.0/7
disc_p <- .01
d_thres <- 5
e_thresh <- 150
prob_symp <- 1
incub_p <- 1
dis_prob_symp <- 0.01
dis_prob_asymp <- 0.00 
intro_rate <- .000


trial_1 <- run_branch(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate)
prob_ext(prop_p = prop_p,recov_p = recov_p, incub_p = incub_p, prob_symp = prob_symp, d_thres,e_thresh = e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=100)
  




