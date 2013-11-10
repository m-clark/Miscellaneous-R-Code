############################################################################################
### The following code relates to a short course given by the Center for Social Research ###
### at Notre Dame.  The course page and associated handout can be found at:              ###
### http://csr.nd.edu/statistical-services/non-credit-short-courses/introduction-to-r/   ###
############################################################################################
testtest
#######################
### Getting Started ###
#######################

UScapt = "D.C."
UScapt

myvar1 = c(1,2,3)
myvar2 = 1:3
myvar3 = seq(from = 1, to=3, by=1)

myvar1
myvar2
myvar3

#library(Rcmdr)

#################################
### Working with the language ###
#################################

### Functions ###
mymean = function(x){
  if (!is.numeric(x)){
    stop("STOP! Does not compute.")
    }
  return(sum(x)/length(x))
}
var1 = 1:5
mymean(var1)
mymean("saywhatnow")


### Loops etc. ###
temp = matrix(rnorm(100), ncol=10)
means = vector()
for (i in 1:ncol(temp)){
  means[i] = mean(temp[,i])
}

means <- colMeans(temp)


### Conditionals ###
x <- c(1,3,5)
y <- c(6,4,2)
z <- ifelse(x < y, "Cake","Death")
z

#######################################
### Importing and working with data ###
#######################################

### Creation of Data ###

x = sample(c("Heads", "Tails", "Edge","Blows Up"), 5, replace=T, prob=c(.45,.45,.05,.05))
x2 = rbinom(5,1,.5)
x3 = rnorm(50, mean=50,sd=10)


cormat=matrix(c(1,.5,.5,1), nrow=2) #type cormat at the command prompt if you want to look at it after creation

library(MASS) # the following function is found in the MASS library

xydata=mvrnorm(40, mu=c(0,0),Sigma=cormat, empirical=T) #empirical = F will produce data that will randomly deviate from the assumed correlation matrix

head(xydata) #take a look at it
cor(xydata)

# Create data set with factor variables of time, subject and gender
mydata <- expand.grid(time = factor(1:4), subject = factor(1:10))
mydata$gender <- factor(rep(0:1, each=4), labels=c("Male", "Female"))
mydata$y <- as.numeric(mydata$time) + 
  rep(rnorm(10,50,10),e=4) +
  rnorm(mydata$time,0,1)
  
library(lattice)
dotplot(y~time, data=mydata,xlab="Time", ylab="DV") #example plot
# Other plot variations
dotplot(y~time,data=mydata, groups=gender, xlab="Time", ylab="DV",
        key = simpleKey(levels(mydata$gender)))
dotplot(y~time|gender,data=mydata, groups=subject, xlab="Time", ylab="DV",
        pch=19, cex=1)


### Importing Data ###

##not run
# mydata <- read.table("c:/trainsacrossthesea.txt", header = TRUE, sep=",",row.names="id")
# 
# library(foreign)
# statadat <- read.dta("c:/rangelife.dta", convert.factors = TRUE)
# spssdat <- read.spss("c:/ourwaytofall.sav", to.data.frame = TRUE)

#example read from web source
mydata <- read.table("http://csr.nd.edu/assets/22641/testwebdata.txt", header=T,sep="")


### Working with the Data ###

# Indexing
state2 <- data.frame(state.x77)
str(state2)  #object structure


head(state2, 10) #first 10 rows
tail(state2, 10) #last 10 rows
state2[,3] #third column
state2[14,] #14th row
state2[3,6] #3rd row, 6th column
state2[,"Frost"] #variable by name

# Data Sets & Subsets
mysubset = subset(state2, state.region== "South") #note that state.region is a separate object

mysubset = state2[state.region== "South",]

mysubset = state2[,c(1:2,7:8)]
mysubset = state2[,c("Population","Income","Frost","Area")]

library(ggplot2)
qplot(x=HS.Grad, y=Murder, data=state2[state.region=="South",])

# Merging and Reshaping
mydat <- data.frame(id=factor(1:12), group=factor(rep(1:2,e=3)) )
x = rnorm(12)
y = sample(70:100,12)
x2 = rnorm(12)

# add columns
mydat$grade = y  #add y via extract operator
df <- data.frame(id=mydat$id,y)
mydat2 <- merge(mydat, df, by="id", sort=F) #using merge
mydat3 <- cbind(mydat,x) #using cbind

# add rows
df <- data.frame(id=factor(13:24), group=factor(rep(1:2,e=3)), grade=sample(y))
mydat2 <- rbind(mydat, df)

# Reshape
mydat <- data.frame(id=factor(rep(1:6,e=2)), time=factor(rep(1:2,6)), y)
mydat2 <- reshape(mydat, v.names="y",direction="wide")
mydat3 <- reshape(mydat2, direction="long")


# Miscellaneous
# Create a factor
gender <- c(rep("male",20), rep("female", 20))
gender <- factor(gender) #levels are now noted, class changes from character to factor

gender <- factor(rep(c(0,1), each=20), labels=c("Male","Female"))


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

tapply(state2$Frost, state.region, describe)
by(state2$Frost, state.region,describe)
describeBy(state2$Frost, state.region)

cor(state2)

### IDA Graphic ###
hist(state2$Illiteracy)
barplot(table(state.region), col=c("lightblue", "mistyrose", "papayawhip", "lavender"))
plot(state2$Illiteracy ~ state2$Population)
stripchart(state2$Illiteracy~state.region, data=state2, col=rainbow(4), method="jitter")


######################
### Basic Modeling ###
######################

mod1 = lm(Income~Illiteracy, data=state2)
summary(mod1)

mod2 = update(mod1, ~. + Population + HS.Grad); summary(mod2)
mod3 = update(mod2, ~. - Population); summary(mod3)
anova(mod1,mod2)
anova(mod3,mod2)


mod1$coef
coef(mod1)
mod1$res
confint(mod1)


####################################
### Visualization of Information ###
####################################
x <- rnorm(100)
hist(x)
boxplot(x)

qplot(depth, data=diamonds,geom="bar",fill=cut, xlim=c(55, 70) ) #used to create margin histogram

library(MASS) #for the data
hist(Cars93$MPG.highway, col="lightblue1", main="Distance per Gallon 1993",xlab="Highway MPG")

par(mfrow=c(1,2))
hist(Cars93$MPG.highway, col="lightblue1",prob= T)
lines(density(Cars93$MPG.highway),col="lightblue4")
boxplot(MPG.highway~Origin, col="burlywood3",data=Cars93)
par(mfrow=c(1,1))

library(car)
scatterplot(prestige~income|type, data=Prestige)
scatterplot(vocabulary~education, jitter=list(x=1,y=1), data=Vocab, col=c("green","red","gray80"), pch=19)
