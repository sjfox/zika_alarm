import scipy

# Three kinds of counters:
# UI - Undiscovered Infecteds
# DI - Discovered Infecteds
# D - Cumulative number of discoveries
#
# Simulate until we have either:
# - we run out of infecteds
# - e_thresh number of infecteds, and d_thresh number of discoveries
#
# prop_p -- probability an I infects a new individual in a time period
# recov_p -- probability an I recovers in a time period
# disc_p -- probability we discover an I
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative)
def run_branch(prop_p, recov_p, disc_p, d_thresh, e_thresh):
    (UI,DI,D) = 1,0,0
    I = UI + DI
    time_record = [(I,D)]
    while ( (I < e_thresh) and (I > 0) ) or ( (I > 0) and (D < d_thresh) ):
        newDI_draws = scipy.rand(DI)
        newDI_count = scipy.sum( newDI_draws < prop_p) 
        decDI_draws = scipy.rand(DI)
        decDI_count = scipy.sum( decDI_draws < recov_p) 

        newUI_draws = scipy.rand(UI)
        newUI_count = scipy.sum( newUI_draws < prop_p) 
        decUI_draws = scipy.rand(UI)
        decUI_count = scipy.sum( decUI_draws < recov_p) 
        disUI_draws = scipy.rand(UI)
        disUI_count = scipy.sum( disUI_draws < disc_p) 

        removeUI = scipy.sum( (decUI_draws < recov_p) | ( disUI_draws < disc_p) )
        removeUDI = scipy.sum( (decUI_draws < recov_p) & ( disUI_draws < disc_p) )

        UI = UI + newUI_count + newDI_count - removeUI
        DI = DI - decDI_count + disUI_count - removeUDI
        D = D + disUI_count
        I = UI + DI
        time_record.append( (I,D) )
    return time_record

# Estimate the probability we will reach e_thresh number of I's (as opposed to
# the infection dying out) given we have discovered at least d_thresh cases
#
# prop_p -- probability an I infects a new individual in a time period
# recov_p -- probability an I recovers in a time period
# disc_p -- probability we discover an I
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative)
def prob_ext(prop_p, recov_p, disc_p, d_thresh, e_thresh, nsamples=10000):
    escapes = 0
    i = 0
    while i < nsamples:
        record = run_branch(prop_p, recov_p, disc_p, d_thresh, e_thresh)
        I,D = record[-1]
        if D < d_thresh:
            continue
        if I > e_thresh:
            escapes += 1
        i += 1 
        print 'Estimate', float(escapes) / i

# R0 is prop_p / recov_p (propagation probability / recovery probability)
# disc_p is probability we discover a case
# disc_threshold is how many cases we have "discovered"
# and e_threshold is the escape number
# estimates the probability escape will happen, if we have discovered at least d_thresh cases
prob_ext(2.0/7, 1.0/7, .01, 5, 150)

