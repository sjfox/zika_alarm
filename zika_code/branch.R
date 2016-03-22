
# prop_p -- probability an I infects a new individual in a time period: how to determine this? 
# recov_p -- probability an I recovers in a time period: human recovery 
# disc_p -- probability we discover an I: symptom ratio 
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative):

#Parameters 
prop_p <- 2/7
recov_p <- 1.0/7
disc_p <- .01
d_thresh <- 5
e_thresh <- 150

# Execute Run
prob_ext(prop_p, recov_p, disc_p, d_thresh, e_thresh)
  
