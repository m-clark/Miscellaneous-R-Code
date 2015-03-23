# A simple function to fill NAs by the last value; not something one should
# regularly do but it does come up from time to time for variables such as age, year
# etc. Just an exercise in recursion.

y = rnorm(100)
y2 = y
y2[sample(100, 30)] = NA

fillbyLast = function(x){
  if (length(x) == 1){
    return(x)
  } else {
    if (is.na(x[2])){
      x[2] = x[1]
    }
    c(x[1], fillbyLast(x[2:length(x)]))
  }
}

y3 = fillbyLast(y2)
cbind(y, y2, y3)

