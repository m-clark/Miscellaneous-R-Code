#################################################################
### See for example: http://c2.com/cgi/wiki?FizzBuzzTest      ###
### The stated goal is to take a sequence and change anything ###
### that is a multiple of 3 to 'Fizz', multiples of 5 to      ###
### 'Fuzz', and any multiples of both to 'FizzBuzz'.          ###
#################################################################


# R example. R makes it easy to generalize to any sequence/multiples and does so very efficiently.

#fizzbuzz = function(min, max, num1, num2){
#  x_ = min:max
#  x = x_
#
#  x[x_%%num1 == 0] = "Fizz"
#  x[x_%%num2 == 0] = "Buzz"
#  x[x_%%(num1*num2) == 0] = "FizzBuzz"
#  x
# }

# fizzbuzz(1, 100, 3, 5)
# fizzbuzz(-50, 50, 7, 4)



# The following is one way to FizzBuzz with a standard loop approach.

function fizzbuzz1(nmin, nmax, num1, num2)
    x0 = [nmin:nmax]
    x = Array(Any,  length(x0))

    for i in 1:length(x0)
        if x0[i]%(num1*num2)==0
            x[i] = "FizzBuzz"
        elseif x0[i]%num2==0
            x[i] = "Buzz"
        elseif x0[i]%num1==0
            x[i] = "Fizz"
        else
            x[i] = x0[i]
        end
    end

    return x
end

print(fizzbuzz1(1, 100, 3, 5)')
print(fizzbuzz1(-50, 50, 7, 4)')



# One can get fizzbuzz without a loop, but initializing x as x0 as in the R
# script will make it one type regardless of the explicit 'Any' declaration
# ('any' type declaration will be overwritten by Int or Float64, even using
# convert). Could add a ! double 'or' (as in the next example) or maybe
# construct a special type, but that's unsatisfactory. No real code efficiency
# is gained over the loop.

function fizzbuzz2(nmin, nmax, num1, num2)
    x0 = [nmin:nmax]
    x = Array(Any, length(x0))

    check1 = x0%num1.==0
    check2 = x0.%num2.==0
    check3 = x0.%(num1*num2).==0
    checkAll = check1 | check2 | check3

    x[check1] = "Fizz"
    x[check2] = "Buzz"
    x[check3] = "FizzBuzz"
    x[!checkAll] = x0[!checkAll]
    return x
end

print(fizzbuzz2(1, 100, 3, 5)')
print(fizzbuzz2(-50, 50, 7, 4)')



# A compromise approach that retains the 'any' type. Aside from the loop it's
# pretty much the same as the R script.

function fizzbuzz3(nmin, nmax, num1, num2)
    x0 = [nmin:nmax]
    x  = Array(Any, length(x0))

    for i in 1:length(x0)
        x[i] = x0[i]
    end

    x[x0%num1.==0] = "Fizz"
    x[x0%num2.==0] = "Buzz"
    x[x0%(num1*num2).==0] = "FizzBuzz"
    return x
end

print(fizzbuzz3(1, 100, 3, 5)')
print(fizzbuzz3(-50, 50, 7, 4)')
