#####################################################################################
### The following function will go to the xkcd website, scrape the html, look for ###
### a particular part based on a div tag, snag it, and write it out to a file for ###
### subsequent text analysis.                                                     ###
###                                                                               ###
### The function requires only one argument, number, indicating the number of     ###
### comic strips.                                                                 ###  
#####################################################################################

### Inspiration from: 
### https://github.com/CabbagesAndKings/xkcd-Topics/blob/master/scripts/getTranscripts.sh

xkcdscrape = function(number){
  # create url string
  url = paste0("http://xkcd.com/", number,'/') 
  
  # creates a single element list argument of unparsed html
  x = scrape(url, parse=F)
  x = as.character(x[[1]]) # not necessary, but simplifies things a bit
  
  # find the point in the text matching the regex
  transcript = regexpr('<div id="transcript" style="display: none">(.*?)</div>', x)
  
  # writes what is sent to console (via cat) to a xkcd#.txt file
  sink(paste0('xkcdRScrape/xkcd', number, '.txt')) # open connection to a file of the name formed by paste0; create the folder first.
  cat(unlist(regmatches(x, transcript))) # send to the console the point specified by transcript in the html text
  sink()  # close connection
}

library(parallel)
cl = makeCluster(3)      # nubmer of cores
clusterEvalQ(cl, library(scrapeR)) # load scrapeR package on the cores
clusterExport(cl, 'xkcdscrape') # export the function to the cores

n = 1283  # number of xkcd comics by late Oct 2013

# a standard non-parallel way to do it without an explicit loop and more straightforward/easier code
# sapply(1:n, xkcdscrape)  # each number 1:n is fed to the function

# the parallized version (timed via system.time function)
system.time({
  parSapply(cl, 1:n, xkcdscrape)
  }) 



### Example loop for comparison; slow
# Check website for number published
# n = 1283
# library(scrapeR)
# system.time({
# for (i in 1:n){
#   url = paste0("http://xkcd.com/", i,'/')
#   x = scrape(url, parse=F)
#   transcript = regexpr('<div id="transcript" style="display: none">(.*?)</div>', x[[1]])
#   sink(paste0('xkcdRScrape/xkcd', i))
#   cat(regmatches(x[[1]], transcript))
#   sink()
# }
# })