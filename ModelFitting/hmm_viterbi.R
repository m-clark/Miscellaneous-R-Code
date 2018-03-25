
# Description -------------------------------------------------------------

# This function duplicates hmm_viterbi.py, which comes from the Viterbi
# algorithm wikipedia page. This first function is just to provide R code that
# is similar in case anyone is interested in comparison, but the original used
# lists of tuples and thus was very inefficient and provided output that wasn't
# succinct.  The second function takes a vectorized approach and returns a
# matrix in a much more straightforward fashion.  Both will provide the same
# result as the Python code.



# Original in R -----------------------------------------------------------

viterbi <- function(obs, states, start_p, trans_p, emit_p) {
  V = vector('list', length(obs))
  
  for (st in seq_along(states)) {
    V[[1]][[states[st]]] = list("prob"= start_p[st] * emit_p[[st]][obs[1]], "prev"= NULL)
  }
  
  for (t in 2:length(obs)) {
    
    for (st in seq_along(states)) {
      max_tr_prob = numeric()
      
      for (prev_st in states) {
        max_tr_prob[prev_st] = V[[t-1]][[prev_st]][["prob"]] * trans_p[[prev_st]][[st]]
      }
      max_tr_prob = max(max_tr_prob)
      
      for (prev_st in states) {
        flag =  V[[t-1]][[prev_st]][["prob"]] * trans_p[[prev_st]][[st]] == max_tr_prob
        if (flag) {
          max_prob = max_tr_prob * emit_p[[st]][obs[t]]
          V[[t]][[states[st]]] = list('prob' = max_prob, 'prev' = prev_st)
        }
      }
    }
  }
  
  # I don't bother duplicating the text output code
  df_out = rbind(Healthy = sapply(V, function(x) x$Healthy$prob),
                 Fever = sapply(V, function(x) x$Fever$prob))
  colnames(df_out) = obs
  print(df_out)
  m = paste0('The steps of states are: ', 
             paste(rownames(df_out)[apply(df_out, 2, which.max)], collapse = ' '), 
             paste('\nHighest probability: ', max(df_out[,ncol(df_out)])))
  message(m)
  V
}


#: Data Setup --------------------------------------------------------------

obs = c('normal', 'cold', 'dizzy')
states = c('Healthy', 'Fever')
start_p = c('Healthy'= 0.6, 'Fever'= 0.4)
trans_p = list(
  'Healthy' = c('Healthy'= 0.7, 'Fever'= 0.3),
  'Fever' = c('Healthy'= 0.4, 'Fever'= 0.6)
)
emit_p = list(
  'Healthy' = c('normal'= 0.5, 'cold'= 0.4, 'dizzy'= 0.1),
  'Fever' = c('normal'= 0.1, 'cold'= 0.3, 'dizzy'= 0.6)
)


#: Demo --------------------------------------------------------------------

test = viterbi(obs, states, start_p, trans_p, emit_p)
# test

set.seed(123)
obs = sample(obs, 6, replace = T)
test = viterbi(obs, states, start_p, trans_p, emit_p)
# test


# Vectorized --------------------------------------------------------------

viterbi_sane <- function(obs, states, start_p, trans_mat, emit_mat) {
  prob_mat = matrix(NA, nrow = length(states), ncol = length(obs))
  colnames(prob_mat) = obs
  rownames(prob_mat) = states

  prob_mat[,1] = start_p * emit_mat[,1]

  for (t in 2:length(obs)) {
    prob_tran = prob_mat[,t-1] * trans_mat
    max_tr_prob =  apply(prob_tran, 2, max)
    prob_mat[,t] =  max_tr_prob * emit_mat[, obs[t]]
  }
  
  print(prob_mat)
  m = paste0('The steps of states are: ', 
             paste(states[apply(prob_mat, 2, which.max)], collapse = ' '), 
             paste('\nHighest probability: ', max(prob_mat[,ncol(prob_mat)])))
  message(m)
}


#: Data Setup --------------------------------------------------------------

# pick data
obs = c('normal', 'cold', 'dizzy')
set.seed(123)
obs = sample(obs, 6, replace = T)

# need matrices now
emit_mat = do.call(rbind, emit_p)
trans_mat = do.call(rbind, trans_p)


#: Demo --------------------------------------------------------------------

viterbi_sane(obs, states, start_p, trans_mat, emit_mat)



# A final demo ------------------------------------------------------------

# from the hidden markov model wikipedia page

states = c('Rainy', 'Sunny')

observations = c('walk', 'shop', 'clean')

start_probability = c('Rainy'= 0.6, 'Sunny'= 0.4)

transition_probability = rbind(
  'Rainy' = c('Rainy'= 0.7, 'Sunny'= 0.3),
  'Sunny' = c('Rainy'= 0.4, 'Sunny'= 0.6)
)

emission_probability = rbind(
  'Rainy' = c('walk'= 0.1, 'shop'= 0.4, 'clean'= 0.5),
  'Sunny' = c('walk'= 0.6, 'shop'= 0.3, 'clean'= 0.1)
)

viterbi_sane(observations, states, start_probability, transition_probability, emission_probability)
