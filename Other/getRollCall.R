### This function allows one to get multiple rollcall data objects (and in 
### parallel) via the readKH function in the pscl package.  See the readKH 
### function helpfile, voteview.com, and also the readme in rollcall folder in 
### the Datasets repo, where the results of this are available for all 
### legislatures. A little messy but works ok.  With 7 threads it took ~90 
### seconds to get all data possible, ~2 minutes for csv.  So if
### you're not requesting quite a few, you wouldn't need to set cores very high.

getRollCall = function(congress, HoRSenBoth='Both', matrix=F, cores=2, Dir=NULL, csv=F){
  # Arguments: 
  # congress- integer vector
  # HoRSenBoth- character, branch of Congress; House, Senate, or Both
  # matrix- logical, if matrix instead of rollcall format is desired. See pscl::convertCodes. See csv arg.
  # cores- integer, number of cores/threads for parallel processing (Windows)
  # Dir- character, where to write out RData files of rollcall objects, or csv if matrix and csv=T
  # csv- logical, write matrix of binary votes to file; only if matrix=T and Dir provided
  
  # requires stringr and pscl packages
  
  # initial checks
  if (csv && !matrix) message('matrix must be TRUE for csv files to be written. *csv files will not be written.')
  if (csv && is.null(Dir)) message('Dir must be specified for csv. *csv files will not be written.')
  
  ### create urls
  # initial check for first 9 since the early files are numbered 01 02 etc.
  congressforurl = as.character(congress)
  congressforurl[congress %in% 1:9] = paste0(0, congressforurl[congress %in% 1:9])
  
  urls = paste0('ftp://voteview.com/', c('hou', 'sen'), rep(congressforurl, e=2), 'kh.ord')

  # deal with current congress as it comes from different website; made to work in subsequent years
  currentYear = as.numeric(substr(Sys.Date(), 1, 4))
  curfutureCongress = data.frame(Year=2013:currentYear, 
                                 Congress=rep(113:200, e=2, length=length(2013:currentYear)))
  currCongress = curfutureCongress[curfutureCongress$Year==currentYear, 'Congress']
  
  current = ifelse(currCongress %in% congress, T, F)
  
  # relace url for current congress with different web address
  if (current) {
    urls[grep(currCongress, urls)] = c(
             paste0('http://adric.sscnet.ucla.edu/rollcall/static/','H', currCongress,'.ord'),
             paste0('http://adric.sscnet.ucla.edu/rollcall/static/','S', currCongress,'.ord'))
  }
  
  # reduce as requested
  if (HoRSenBoth=='House') {
    urls = urls[grep('hou|H', urls)]
  } else if (HoRSenBoth=='Senate') {
    urls = urls[grep('sen|S', urls)]
  }
  
  
  ### autocreate descriptions for objects and files just for giggles
  firsts = c(1, as.numeric(paste0(2:10, 1)))
  seconds = c(2, as.numeric(paste0(2:10, 2)))
  thirds = c(3, as.numeric(paste0(2:10, 3)))
  
  suppressPackageStartupMessages(require(stringr))
  descriptions = rep(NA, length(urls))
  descripN = as.numeric(str_extract(urls, '[0-9]+|[0-9]+'))
  descripB = str_extract(urls, 'hou|H|sen|S')
  
  if (any(descripN %in% c(firsts, seconds, thirds))){
    idx1 = which(descripN %in% firsts)
    idx2 = which(descripN %in% seconds)
    idx3 = which(descripN %in% thirds)
    
    ends = c(rep('st', length(idx1)), rep('nd', length(idx2)), rep('rd', length(idx3)))
    idx = c(idx1, idx2, idx3)
    descriptions = rep(NA, length(urls))
    descriptions[idx] = paste0(descripN[idx], ends, ' U.S. ', 
                               ifelse(descripB[idx] %in% c('hou','H'), 'House of Representatives', 'Senate')) 
  }
  
  idx = is.na(descriptions)
  descriptions[idx] = paste0(descripN[idx], 'th', ' U.S. ', 
                             ifelse(descripB[idx] %in% c('hou','H'), 'House of Representatives', 'Senate'))
  
  
  ### set up parallel
  library(parallel)
  cl = makeCluster(cores)
  clusterEvalQ(cl, library(pscl))

  suppressPackageStartupMessages(require(pscl))
  rcList = clusterMap(cl, fun=readKH, urls, desc=descriptions)
  
  rcenv = new.env() # so export doesn't look in global
  
  ### write out RData files
  if (!is.null(Dir) & !csv){
    filelist = paste0(Dir, descriptions, '.RData')
    descriptionsBrief = paste0(descripB, descripN)
    
    clusterExport(cl, c('rcList','filelist', 'descriptionsBrief'), envir = rcenv)

    parSapply(cl, 1:length(rcList), function(i) {
      assign(descriptionsBrief[i], rcList[[i]]); save(list=descriptionsBrief[i], file=filelist[i])
    })
  }
  
  
  ### convert to binary 
  if (matrix){
    clusterExport(cl, c('rcList'), envir=rcenv)
    rcList = parLapply(cl, rcList, convertCodes)
    if (csv && !is.null(Dir)){
      filelist = paste0(Dir, descriptions, '.csv')
      clusterMap(cl, write.csv, rcList, filelist)
    } 
  }
  
  on.exit(stopCluster(cl))
  rcList
}
