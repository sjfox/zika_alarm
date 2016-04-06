###################################################
## Determining incubation period boxes and infectious period boxes
## simulate distributions
# incubation_time <- function(nboxes, prob){
#   box = 1
#   t = 0
#   while(box <= nboxes){
#     if(runif(1) < prob){
#       box = box+1
#     }
#     t = t + 1
#   }
#   t
# }
# rincubation <- function(num_reps, ...) {
#   raply(.n = num_reps, .expr = incubation_time(...) )
# }

fitDist <- function(prob, size, mean) {
  draws <- rnbinom(n=100000, size = size, prob = prob) + size
  mean - mean(draws)
}

## Look at weibull distribution from lessler 
qweibull(c(0.025, 0.5, 0.975), shape=1.97, scale=10.87)
hist(rweibull(100000, shape=1.97, scale=10.87), breaks=100)

## First fit the weibll viral clearance distribution
## Lessler et al biorxiv key times paper
## mean = 9.88, cI : 
n_boxes <- 3
infectious_best_fit <- uniroot(f = fitDist, interval = c(0.001,.999), size= n_boxes, mean=9.88)
infectious_times <- rnbinom(100000, size = n_boxes, prob = infectious_best_fit$root) + n_boxes
hist(infectious_times, breaks=50)
quantile(infectious_times, c(.025, .5, .975))

## Now determine the incubation period boxes
## Use the serial interval from brownstein paper
## 10-23 days 
fitDistBoth <- function(prob, size, min, max, infectious_times) {
  draws <- rnbinom(n=100000, size = size, prob = prob) + size
  qs <- quantile((draws+infectious_times/2), probs = c(.025, .975))
  sum(qs - c(min,max))^2
}
inc_boxes = 6
incubation_best_fit <- optimize(f = fitDistBoth, interval = c(0.001,.999), size= inc_boxes, min=10, max=23, infectious_times=infectious_times)
incubation_times <- rnbinom(n = 100000, size = inc_boxes, prob = incubation_best_fit$minimum) + inc_boxes
quantile((incubation_times+infectious_times/2), probs = c(.025, .5, .975))
hist(incubation_times+infectious_times/2, breaks=100, xlim=c(0,45))

infectious_best_fit$root
incubation_best_fit$minimum

mean(infectious_times)

