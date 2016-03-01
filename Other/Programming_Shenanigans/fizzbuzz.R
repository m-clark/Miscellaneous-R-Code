#################################################################
### See for example: http://c2.com/cgi/wiki?FizzBuzzTest      ###
### The stated goal is to take a sequence and change anything ###  
### that is a multiple of 3 to 'Fizz', multiples of 5 to      ###
### 'Fuzz', and any multiples of both to 'FizzBuzz'.          ###
### R makes it easy to generalize to any sequence/multiples   ###
### and does so very efficiently.                             ###
#################################################################


## clean, clear, gets the job done; also generalizes beyond common example
fizzbuzz = function(min, max, num1, num2){
  x_ = min:max
  x = x_
  x[x_%%num1 == 0] = "Fizz"
  x[x_%%num2 == 0] = "Buzz"
  x[x_%%(num1*num2) == 0] = "FizzBuzz"
  x
}

fizzbuzz(1, 100, 3, 5)
fizzbuzz(-50, 50, 7, 4)

## One-liner, just because you can
fizzbuzzOneLine = function(min, max, num1, num2){
  sapply(min:max, function(val) ifelse(val%%(num1*num2) == 0, "FizzBuzz", 
                                       ifelse(val%%num1==0, 'Fizz', 
                                              ifelse(val%%num2==0, 'Buzz', val))))
}

fizzbuzzOneLine(1, 100, 3, 5)

## recursive for even more flexibility
fizzbuzzRecursive = function(x, nums=c(15,3,5), nams=c('FizzBuzz','Fizz', 'Buzz')){
  if (length(nums) == 0) {paste(x)}
  else {
    x = sapply(x, function(val) ifelse(!is.na(suppressWarnings(as.numeric(val))) && as.numeric(val) %% nums[1] == 0, nams[1], val))
    fizzbuzzRecursive(x, nums[-1], nams[-1])
  }
}

fizzbuzzRecursive(1:100, nums=c(15,3,5), nams=c('FizzBuzz','Fizz', 'Buzz'))


cbind(fizzbuzz(1, 100, 3, 5), fizzbuzzOneLine(1, 100, 3, 5), fizzbuzzRecursive(1:100, nums=c(15,3,5), nams=c('FizzBuzz','Fizz', 'Buzz')))
cbind(fizzbuzz(-50, 50, 7, 4), fizzbuzzOneLine(-50, 50, 7, 4), fizzbuzzRecursive(-50:50, nums=c(28,7,4), nams=c('FizzBuzz','Fizz', 'Buzz')))

fizzbuzzRecursive(-50:50, nums=c(28, 21, 12, 7, 4, 3), nams=c('FizzBuzz','FizzKapow', 'BuzzKapow', 'Fizz', 'Buzz', 'Kapow'))

