#################################################################
### See for example: http://c2.com/cgi/wiki?FizzBuzzTest      ###
### The stated goal is to take a sequence and change anything ###  
### that is a multiple of 3 to 'Fizz', multiples of 5 to      ###
### 'Fuzz', and any multiples of both to 'FizzBuzz'.          ###
### R makes it easy to generalize to any sequence/multiples   ###
### and does so very efficiently.                             ###
#################################################################



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

