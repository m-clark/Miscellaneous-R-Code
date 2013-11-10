#####################################################################################
### The following function will go to the xkcd website, scrape the html, look for ###
### a particular part based on a div tag, snag it, and write it out to a file for ###
### subsequent text analysis.                                                     ###
#####################################################################################

### Inspiration from: 
### https://github.com/CabbagesAndKings/xkcd-Topics/blob/master/scripts/getTranscripts.sh
### See the R version in my repo.  Being new to Python I thought this would be a nice challenge.
### I make no claims as to its efficiency.

###############
### Preface ###
###############
import os
os.getcwd()

# change working directory to wherever you want the text files to be stored
os.chdir('C:/Users/mclark19/Desktop/CSR/Clients/Me/xkcdscrape/')

import time # for time comparisons if desired
import urllib2
import BeautifulSoup  

Soup = BeautifulSoup.BeautifulSoup 

###############
### Example ###
###############
### read in web contents
url = urllib2.urlopen("http://www.xkcd.com/1/")
xkcd1 = url.read()

### create a beautifulSoup object
out = Soup(xkcd1)

### Inspect in a standard html form
print(out.prettify())

### example of primary process: extract the div with id='transcript'
out.find('div', {'id' : 'transcript'}) 

######################################
### Main Process via Standard Loop ###
######################################

### Preliminaries
# number of xkcd comics
n = 1283

# test n
# n = 100

### set initial time
t0 = time.time() 

### Run the loop over all comics
for i in xrange(1, n+1):
    filename = 'xkcdpyScrape/xkcd' +  str(i) + '.txt'
    if i == 404:    # this deals with an xkcd joke for comic #404
        out = " "
    else:
        # read in the web contents
        url = urllib2.urlopen('http://www.xkcd.com/' + str(i))
        xkcd = url.read()
        
        # create a soup object
        xkcdsoup = Soup(xkcd)
        
        # scrape the div with id = 'transcript'
        out = xkcdsoup.find('div', {'id' : 'transcript'})
    
    # Write out the file
    f = open(filename, 'w+')
    f.write(str(out))
    f.close()

time.time() - t0  # roughly 20 seconds for 100, about 4.5 minutes for all

####################
### Parallelized ###
####################

##################################################################################
### Create scrapexkcd function that will do what's in the loop above given a   ###
### start and end point. See the loop above for details. Arguments include the ###
### start and end number specific to the comics one desires to download.       ###
##################################################################################

def scrapexkcd(start, end):
    import urllib2
    import BeautifulSoup 
    for i in xrange(start, end):
        filename = 'xkcdpyScrape/parallel/xkcd' +  str(i) + '.txt'
        if i == 404:
            out = " "
        else:
            url = urllib2.urlopen('http://www.xkcd.com/' + str(i))
            xkcd = url.read()
            xkcdsoup = BeautifulSoup.BeautifulSoup(xkcd)
            out = xkcdsoup.find('div', {'id' : 'transcript'})
        f = open(filename, 'w+')
        f.write(str(out))
        f.close()

start = 1
end = 1283 + 1 

# test end
end = 100

import pp # for parallelization

# set number of workers
ncpus = 3

# Create jobserver
job_server = pp.Server(ncpus=ncpus)

# Divide the task into subtasks
parts = 3
step = (end - start) / parts + 1

jobs = []

t0 = time.time()

for index in xrange(parts):
    # create start and endpoints for the function
    starti = start + index*step
    endi = min(start + (index+1)*step, end)
    # Submit a job which will scrape the site
    job_server.submit(scrapexkcd, (starti, endi))

#wait for jobs in all groups to finish 
job_server.wait()

time.time() - t0  # time since t0

# print stats
job_server.print_stats()