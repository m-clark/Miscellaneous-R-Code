# The following provides some conceptual code for the Chinese restaurant and 
# Indian buffet process for categorical and continuous latent variables 
# respectively.  For more detail, see the Bayesian nonparametric section of my
# structural equation modeling document. 
# https://m-clark.github.io/docs/sem/bayesian-nonparametric-models.html


# To start, a couple functions demonstrating the Chinese restaurant processes.  The first
# is succinct and more conceptual, but notably slower.


crp <- function(alpha, n) {
  table_assignments = 1
  
  for (i in 2:n){
    table_counts = table(table_assignments)       # counts of table assignments
    nt = length(table_counts)                     # number of tables  
    table_prob = table_counts/(i-1+alpha)         # probabilities of previous table assignments
    
    # sample assignment based on probability of current tables and potential next table
    current_table_assignment = sample(1:(nt+1), 1, prob=c(table_prob, 1-sum(table_prob)))
    
    # concatenate new to previous table assignments
    table_assignments = c(table_assignments, current_table_assignment)  
  }
  table_assignments
}


# similar to the restaurant function here https://github.com/mcdickenson/shinyapps
crpF <- function(alpha, n) {
  table_assignments = c(1, rep(NA, n-1))
  table_counts = 1
  
  for (i in 2:n){
    init =  c(table_counts, alpha)
    table_prob = init/sum(init)
    current_table_assignment = sample(seq_along(init), 1, prob=table_prob)
    table_assignments[i] = current_table_assignment
    if(current_table_assignment==length(init)){
      table_counts[current_table_assignment] = 1
    } else {
      table_counts[current_table_assignment] = table_counts[current_table_assignment] + 1
    }
  }
  table_assignments
}




# zoom these to see properly
par(mfrow=c(5, 5))
replicate(25, {
  out = crp(alpha=1, n=500)
  plot(table(out), bty='n', lwd=5)
})

replicate(25, {
  out = crpF(alpha=1, n=500)
  plot(table(out), bty='n', lwd=5)
})
layout(1)



# library(microbenchmark)
# test  = microbenchmark(crp(alpha = 1, n=1000), 
#                        crpF(alpha = 1, n=1000), times=100)
# test
# ggplot2::autoplot(test)

library(tidyverse)
n = 100
crp_1 = crp(alpha=1, n=n)
crp_4 = crp(alpha=4, n=n)
crp_1_mat = matrix(0, nrow=n, ncol=n_distinct(crp_1))

for (i in 1:n_distinct(crp_1)){
  crp_1_mat[, i] = ifelse(crp_1==i, 1, 0)
}


crp_4_mat = matrix(0, nrow=n, ncol=n_distinct(crp_4))

for (i in 1:n_distinct(crp_4)){
  crp_4_mat[, i] = ifelse(crp_4==i, 1, 0)
}



d3heatmap::d3heatmap(crp_1_mat, Rowv=F, Colv=F, yaxis_font_size='0pt', colors='Blues', width=400)
d3heatmap::d3heatmap(crp_4_mat, Rowv=F, Colv=F, yaxis_font_size='0pt', colors='Blues', width=400)



# The following demonstrates the Indian buffet process for continuous latent
# variable settings.

ibp = function(alpha, N){
  # preallocate assignments with upper bound of N*alpha number of latent factors
  assignments = matrix(NA, nrow=N, ncol=N*alpha) 
  
  # start with some dishes/assigments
  dishes = rpois(1, alpha)      
  zeroes = ncol(assignments) - dishes   # fill in the rest of potential dishes
  assignments[1,] = c(rep(1, dishes), rep(0, zeroes))
  
  for(i in 2:N){
    prev = i-1
    # esoteric line that gets the last dish sampled without a search for it
    last_previously_sampled_dish = sum(colSums(assignments[1:prev,,drop=F]) > 0)    
    
    # initialize 
    dishes_previously_sampled = matrix(0, nrow=1, ncol=last_previously_sampled_dish)
    
    # calculate probability of sampling from previous dishes
    dish_prob = colSums(assignments[1:prev, 1:last_previously_sampled_dish, drop=F]) / i
    dishes_previously_sampled[1,] = rbinom(n=last_previously_sampled_dish, size=1, prob=dish_prob)
    
    # sample new dish and assign based on results
    new_dishes = rpois(1, alpha/i)
    zeroes = ncol(assignments) - (last_previously_sampled_dish + new_dishes)
    assignments[i,] = c(dishes_previously_sampled, rep(1,new_dishes), rep(0, zeroes))
  }
  
  # return only the dimensions sampled
  last_sampled_dish = sum(colSums(assignments[1:prev,]) > 0) 
  return(assignments[, 1:last_sampled_dish])
}



set.seed(123)
ibp_1 = ibp(1, 100)
ibp_4 = ibp(4, 100)

d3heatmap::d3heatmap(ibp_1, Rowv=F, Colv=F, yaxis_font_size='0pt', colors='Blues', width=400)
d3heatmap::d3heatmap(ibp_4, Rowv=F, Colv=F, yaxis_font_size='0pt', colors='Blues', width=400)
