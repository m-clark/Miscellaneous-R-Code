## An implementation of Gibbs sampling for topic models for the example in 
## section 4 of Steyvers and Griffiths (2007): 
## http://cocosci.berkeley.edu/tom/papers/SteyversGriffiths.pdf; a very clear 
## intro in my opinion. The core of the functions's code comprises mostly 
## cosmetic changes to that found here: 
## https://bcomposes.wordpress.com/2012/05/29/gibbs-sampler-for-toy-topic-model-example/
## Added are the creation of a function with several arguments, plotting etc. 
## Packages required to run everything here are corrplot, abind, coda, reshape2,
## and ggplot2

##################
### Data Setup ###
##################
vocab = factor(c("river","stream","bank","money","loan"))
K = 2  # n of topics
v = length(vocab)  # number of unique words
d = 16 # number of documents 

# topic 1 gives equal probability to money loan and bank (zero for river and stream)
# topic 2 gives equal probability to river stream and bank (zero for money and loan)

# document term matrix; each doc consists of a mix of 16 tokens of the vocab;
# the first few regard financial banks, the last several water banks, and the
# docs in between possess a mixed vocab
dtm = matrix(c(0,0,4,6,6,
               0,0,5,7,4,
               0,0,7,5,4,
               0,0,7,6,3,
               0,0,7,2,7,
               0,0,9,3,4,
               1,0,4,6,5,
               1,2,6,4,3,
               1,3,6,4,2,
               2,3,6,1,4,
               2,3,7,3,1,
               3,6,6,1,0,
               6,3,6,0,1,
               2,8,6,0,0,
               4,7,5,0,0,
               5,7,4,0,0), ncol=v, byrow=T)

row.names(dtm) = paste0('doc', 1:d); colnames(dtm) = vocab

# matrix of words in each document
wordmat = t(apply(dtm, 1, function(row) rep(vocab, row)))

# initialize random topic assignments to words
T0 = apply(wordmat, c(1,2), function(token) sample(1:2, 1))

# word by topic matrix of counts containing the number of times word w is assigned to topic j
C_wt = t(sapply(vocab, function(word) cbind(sum(T0[wordmat == word] == 1), 
                                            sum(T0[wordmat == word] == 2))
                ))
row.names(C_wt) = vocab

# topic by document matrix of counts containing the number of times topic j is assigned to a word in document d
C_dt = t(apply(T0, 1, table))


################
### Function ###
################
# note that this function is not self contained in that it uses some of the
# objects created above

tmod = function(alpha, beta, nsim=2000, warmup=nsim/2, thin=10, verbose=T, plot=F, dotsortext='dots'){
  ### Arguments: 
  # alpha and beta- hyperparameters
  # nsim- the total number of simulations
  # warmup- the number of initial sims to discard
  # thin- keep every thin simulation after warmup 
  # verbose- to print every 100th iteration and total time
  # plot- to visualize every x sim; the plot will show current estimates of 
  # theta, phi, and topic assignments for each term in the docs; can be
  # visualized as dots or the actual terms; for the latter you may need to
  # fiddle with your viewing area size so that terms don't appear to run
  # together; plotting will increase runtime
  
  ### requires corrplot and abind packages.
  
  # initialize 
  Z = T0                                                   # topic assignments
  saveSim = seq(warmup+1, nsim, thin)                      # iterations to save
  thetaList = list()                                       # saved theta estimates
  phiList = list()                                         # saved phi estimates
  p = proc.time()
  # for every simulation
  for (s in 1:nsim) {
    if(verbose && s %% 100 == 0) {                         # paste every 100th iteration if desired
      secs = round((proc.time()-p)[3],2)
      min = round(secs/60, 2)
      message(paste0('Iteration number: ', s, '\n', 'Total time: ', 
                     ifelse(secs <= 60, paste(secs, 'seconds'), 
                            paste(min, 'minutes'))))
    }
    if(plot > 0 && s >1 &&  s %% plot == 0){                  # plot every value of argument
      require(corrplot)
      
      layout(matrix(c(1,1,2,2,3,3,3,3,3,3), ncol=10))
      corrplot(theta, is.corr=F, method='color', tl.cex=.75, 
               tl.col='gray50', cl.pos='n', addgrid=NA)
      corrplot(phi, is.corr=F, method='color', tl.cex=1, 
               tl.col='gray50', cl.pos='n', addgrid=NA)
      if(dotsortext == 'dots'){
        Zplot = Z
        Zplot[Zplot==2] = -1
        corrplot(Zplot, is.corr=F, method='circle', tl.cex=.75, 
                 tl.col='gray50', cl.pos='n', addgrid=NA) 
      } else {
        cols = apply(Z, c(1,2), function(topicvalue) ifelse(topicvalue==1, '#053061', '#67001F'))
        plot(1:nrow(wordmat), 1:ncol(wordmat), type="n", axes=F, xlab='', ylab='')
        text(col(wordmat), rev(row(wordmat)), wordmat,  col = cols, cex=.75)
      }
    }
      
    
    # for every document and every word in the document
    for (i in 1:d) {
      for (j in 1:length(wordmat[i,])) {
        word.id = which(vocab == wordmat[i,j])
        topic.old = Z[i,j]
        
        # Decrement counts before computing equation (3) in paper noted above
        C_dt[i, topic.old] = C_dt[i, topic.old] - 1
        C_wt[word.id, topic.old] = C_wt[word.id, topic.old] - 1
        
        # Calculate equation (3) for each topic
        vals = prop.table(C_wt + beta, 2)[word.id,] * prop.table(C_dt[i,] + alpha)
        
        # Sample the new topic from the results for (3);         
        # note, sample function does not require you to have the probs sum to 1
        # explicitly, i.e.  prob=c(1,1,1) is the same as prob=c(1/3,1/3,1/3)
        Z.new = sample(1:K, 1, prob=vals)
        
        # Set the new topic and update counts
        Z[i,j] = Z.new
        C_dt[i, Z.new] = C_dt[i, Z.new] + 1
        C_wt[word.id, Z.new] = C_wt[word.id, Z.new] + 1
      }
    }
    
    theta = prop.table(C_dt + alpha,1)                     # doc topic distribution
    phi = prop.table(C_wt + beta,2)                        # word topic distribution
    
    # save simulations
    if (s %in% saveSim){
      thetaList[[paste(s)]] = theta
      phiList[[paste(s)]] = phi
    }
  }
  
  layout(1) # reset plot window
  
  # value
  results = list(theta=theta, phi=phi, 
                 thetaSims = abind::abind(thetaList, along=3),    # abind creates arrays from the lists
                 phiSims = abind::abind(phiList, along=3)
                 )
}



###########
### Run ###
###########
# values of alpha and beta as suggested in paper
alpha = K/50
beta = .01

topicModelDemo = tmod(alpha = alpha, beta = beta, nsim = 5500, warmup=500, thin=10, 
                      plot=5)



###############
### Results ###
###############
str(topicModelDemo, 1)
sapply(topicModelDemo[1:2], round, 2)

## compared to reported paper estimates for phi and topicmodels package (note
## delta is beta above; the vignette actually references Steyvers & Griffiths
## paper)
library(topicmodels)
ldaout = LDA(dtm, k=2, method='Gibbs', control=list(alpha=alpha, delta=.01, iter=5500, 
                                                    burnin=500, thin=10, initialize='random'))
phiEstimatesCompared = data.frame(round(topicModelDemo$phi, 2), 
                                  paper=cbind(c(0,0,.39,.32,.29), c(.25,.4,.35,0,0)),
                                  ldapack=round(t(posterior(ldaout)$terms),2))
phiEstimatesCompared


## interval estimates
library(coda)
phiEstimatesTopic1 = as.mcmc(t(topicModelDemo$phiSims[,1,]))
phiEstimatesTopic2 = as.mcmc(t(topicModelDemo$phiSims[,2,]))
lapply(list(Topic1 = round(phiEstimatesTopic1,2), Topic2 = round(phiEstimatesTopic2, 2)), HPDinterval)
summary(phiEstimatesTopic1)

## symmetrized or mean Kullback Liebler divergence for topic (dis)similarity
KLdivs = .5*sum(apply(topicModelDemo$phi, 1, function(row) row[1]*log2(row[1]/row[2]) + row[2]*log2(row[2]/row[1])))
KLdivs

# can compare to following, see various other packages for KL-divergence
# LaplacesDemon::KLD(topicModelDemo$phi[,1], topicModelDemo$phi[,2])$mean.sum.KLD
# mean(c(entropy::KL.empirical(topicModelDemo$phi[,1], topicModelDemo$phi[,2]),
#      entropy::KL.empirical(topicModelDemo$phi[,2], topicModelDemo$phi[,1])))

## symmetrized or mean Kullback Liebler divergence for document (dis)similarity
# example
docs2compare = c(1,16)
.5*sum(apply(topicModelDemo$theta[docs2compare,], 2, 
             function(topicprob) topicprob[1]*log(topicprob[1]/topicprob[2]) + 
               topicprob[2]*log(topicprob[2]/topicprob[1])))

# create a function to work on theta to produce Kullback-Liebler or Jensen-Shannon divergence
divKLJS = function(input, method='KL'){
  if(method=='KL'){
    .5*sum(apply(topicModelDemo$theta[input,], 2, 
                 function(topicprob) topicprob[1]*log2(topicprob[1]/topicprob[2]) + 
                   topicprob[2]*log2(topicprob[2]/topicprob[1])
                 )
           )
  } else{
    .5*sum(apply(topicModelDemo$theta[input,], 2, 
                 function(topicprob) topicprob[1]*log2(topicprob[1]/mean(topicprob)) + 
                   topicprob[2]*log2(topicprob[2]/mean(topicprob))
                 )
           )
  }
}

# Now do for all
docpairs = combn(1:d, 2)
KLdivs = apply(docpairs, 2, divKLJS)
mat0 = matrix(0, d, d)
mat0[lower.tri(mat0)] <- KLdivs
KLdivs = mat0; dimnames(KLdivs) = list(row.names(dtm))
corrplot(KLdivs, diag=T, type='lower', is.corr=F, main='Kullback-Liebler Divergence')


# Jensen Shannon divergence;
JSdivs = apply(docpairs, 2, divKLJS, method='JS')
mat0[lower.tri(mat0)] <- JSdivs
JSdivs = mat0; dimnames(JSdivs) = list(row.names(dtm))
corrplot(JSdivs, diag=T, type='lower', is.corr=F, main='Jensen-Shannon Divergence')

# comparison: not recommended to install cummeRbund just for this as it requires ~ dozen or more dependencies
# corrplot(as.matrix(cummeRbund::JSdist(t(topicModelDemo$theta))), diag=T, type='lower', is.corr=F)

par(mar=c(5, 4, 4, 2) + 0.1) # reset plot margins



########################
### Diagnostics etc. ###
########################
## examples of diagnostics with the word topic probability estimates; these are from coda
plot(phiEstimatesTopic1, ask=F)
acfplot(phiEstimatesTopic1)
cumuplot(phiEstimatesTopic1)


## traceplots using ggplot2
library(reshape2)
gphi = melt(topicModelDemo$phiSims)
colnames(gphi) = c('word', 'topic', 'iter', 'value')
str(gphi)

library(ggplot2)
ggplot(aes(x=iter, y=value, color=factor(topic), group=factor(topic)), data=gphi) +
  geom_line(alpha=.5, show_guide=F) +
  facet_grid(~word) +
  theme_minimal()

  
