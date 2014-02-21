#######################
### Getting Started ###
#######################

# create a simple object
UScapt = 'D.C.'
UScapt

# 3 ways to create the same object
myvar1 = c(1, 2, 3)
myvar2 = 1:3
myvar3 = seq(from = 1, to=3, by=1)

myvar1
myvar2
myvar3


#################################
### Working with the language ###
#################################

### Functions ###
# Simple function to calculate the mean
mymean = function(x){
  if (!is.numeric(x)){                 # If x is not numeric...
    stop('STOP! Does not compute.')    # Stop with message
    }
  
  return(sum(x)/length(x))             # Otherwise calculate the mean
}

var1 = 1:5
mymean(var1)
mymean('thisisnotanumber')             # error!


### Loops etc. ###
temp = matrix(rnorm(100), ncol=10)     # a 10 x 10 matrix with values from N(0,1)
means = vector()                       # Initialize the means vector

for (i in 1:ncol(temp)){
  means[i] = mean(temp[, i])           # for the ith column, calculate its mean, assign to the ith element of means vector
}

means <- colMeans(temp)                # much easier/better


### Conditionals ###
x <- c(1, 3, 5)
y <- c(6, 4, 2)
z <- ifelse(x < y, 'Cake', 'Death')    # If x is less than y, z = 'Cake', otherwise 'Death'
z

#######################################
### Importing and working with data ###
#######################################

### Creation of Data ###

x = sample(c('Heads', 'Tails', 'Edge', 'Blows Up'), 5, 
           replace=T, prob=c(.45, .45, .05, .05))
x2 = rbinom(5, 1, .5)
x3 = rnorm(50, mean=50, sd=10)


cormat=matrix(c(1, .5, .5, 1), nrow=2) #type cormat at the command prompt if you want to look at it after creation

library(MASS) # the following function is found in the MASS library

xydata=mvrnorm(40, mu=c(0, 0), Sigma=cormat, empirical=T) #empirical = F will produce data that will randomly deviate from the assumed correlation matrix

head(xydata) #take a look at it
cor(xydata)

# Create data set with factor variables of time, subject and gender
mydata <- expand.grid(time = factor(1:4), subject = factor(1:10))      # creates a data.frame with crossed factors of time and subject
mydata$gender <- factor(rep(0:1, each=4), labels=c('Male', 'Female'))  # add a gender factor
mydata$y <- as.numeric(mydata$time) +                                  # add a response with subject specific effect
  rep(rnorm(10, 50, 10), e=4) +
  rnorm(mydata$time, 0, 1)
  
library(lattice)
dotplot(y ~ time, data=mydata, xlab='Time', ylab='DV') #example plot

# Other plot variations
dotplot(y ~ time, data=mydata, groups=gender, xlab='Time', ylab='DV',
        key = simpleKey(levels(mydata$gender)))
dotplot(y ~ time|gender, data=mydata, groups=subject, xlab='Time', ylab='DV',
        pch=19, cex=1)


### Importing Data ###

## not run, examples only
# mydata <- read.table('c:/trainsacrossthesea.txt', header = TRUE, sep=',', row.names='id')
# 
# library(foreign)
# statadat <- read.dta('c:/rangelife.dta', convert.factors = TRUE)
# spssdat <- read.spss('c:/ourwaytofall.sav', to.data.frame = TRUE)

#example read from web source
mydata <- read.table('http://csr.nd.edu/assets/22641/testwebdata.txt', header=T, sep='')


### Working with the Data ###

# Indexing
state2 <- data.frame(state.x77)
str(state2)  #object structure


head(state2, 10) #first 10 rows
tail(state2, 10) #last 10 rows
state2[,3]       #third column
state2[14,]      #14th row
state2[3,6]      #3rd row, 6th column
state2[,'Frost'] #variable by name

# Data Sets & Subsets
mysubset = subset(state2, state.region== 'South') #note that state.region is a separate R object

mysubset = state2[state.region== 'South',]        # alternate method

mysubset = state2[,c(1:2, 7:8)]                                # grab columns 1 and 2 and 7 and 8
mysubset = state2[,c('Population', 'Income', 'Frost', 'Area')] # grab columns by name
mysubset = state2[grep('^I.*a$', rownames(state2)),]           # get any States starting with 'I' and ends with 'a' using regular expressions
head(mysubset)


library(ggplot2)
qplot(x=HS.Grad, y=Murder, data=state2[state.region=='South',])

# Merging and Reshaping

mydat <- data.frame(id=factor(1:12), group=factor(rep(1:2, e=3)) )
x = rnorm(12)
y = sample(70:100, 12)


# add columns
mydat$grade = y  # add y via extract operator
df <- data.frame(id=mydat$id, y)
mydat2 <- merge(mydat, df, by='id', sort=F) # using merge
mydat3 <- cbind(mydat, x)                   # using cbind

# add rows
df <- data.frame(id=factor(13:24), group=factor(rep(1:2, e=3)), grade=sample(y))
mydat2 <- rbind(mydat, df)

# Reshape (I actually prefer the functions in the reshape2 package for this sort of thing)
mydata <- data.frame(id=factor(rep(1:6, e=2)), time=factor(rep(1:2, 6)), y)
mydataWide <- reshape(mydata, v.names='y', direction='wide')
mydataLong <- reshape(mydataWide, direction='long')


# Miscellaneous
# Create a factor
gender <- c(rep('male', 20), rep('female', 20))
gender <- factor(gender) #levels are now noted, class changes from character to factor

gender <- factor(rep(c(0, 1), each=20), labels=c('Male', 'Female'))


#############################
### Initial Data Analysis ###
#############################

### IDA Numeric ###
mean(state2$Life.Exp)
sd(state2$Life.Exp)
summary(state2$Population)
table(state.region) #region is a separate vector in the base R environment

library(psych)
describe(state2)

# Get statistics by region using different approaches
tapply(state2$Frost, state.region, describe)   
by(state2$Frost, state.region, describe)
describeBy(state2$Frost, state.region)

# Correlation matrix
cor(state2)

### IDA Graphic ###
# histogram
hist(state2$Illiteracy)

# bar plot
barplot(table(state.region), col=c('lightblue', 'mistyrose', 'papayawhip', 'lavender'))

# scatter plot
plot(state2$Illiteracy ~ state2$Population)

# stripchart
stripchart(state2$Illiteracy~state.region, data=state2, col=rainbow(4), method='jitter')


######################
### Basic Modeling ###
######################

# income predicted by illiteracy rate
mod1 = lm(Income ~ Illiteracy, data=state2)
summary(mod1)  # model summary

# add population and high school graduation rate; ~. means keep all predictors from mod1
mod2 = update(mod1, ~. + Population + HS.Grad)
summary(mod2)

# remove population from previous model
mod3 = update(mod2, ~. - Population)
summary(mod3)

# compare models 1 and 2
anova(mod1, mod2)

# compare models 2 and 3
anova(mod3, mod2)

# extract coefficients, residuals, get confidence intervals, plot model diagnostics
mod1$coef
coef(mod1)
mod1$res
confint(mod1)
par(mfrow=c(2,2))
plot(mod1, ask=F)
par(mfrow=c(1,1))

####################################
### Visualization of Information ###
####################################
x <- rnorm(100)
hist(x)
boxplot(x)

qplot(depth, data=diamonds, geom='bar',fill=cut, xlim=c(55, 70) ) #used to create margin histogram

library(MASS) #for the data

hist(Cars93$MPG.highway, col='lightblue1', main='Distance per Gallon 1993', xlab='Highway MPG')

par(mfrow=c(1,2))
hist(Cars93$MPG.highway, col='lightblue1', prob= T)

# add density lines
lines(density(Cars93$MPG.highway), col='lightblue4')

# add boxplot as second graphic
boxplot(MPG.highway ~ Origin, col='burlywood3', data=Cars93)
par(mfrow=c(1,1))

# get more from your scatter plots
library(car)
scatterplot(prestige ~ income|type, data=Prestige)
scatterplot(vocabulary ~ education, jitter=list(x=1,y=1), data=Vocab, col=c('green','red','gray80'), pch=19)
