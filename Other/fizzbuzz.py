# Python Fizz Buzz. Replicate of the R code.

def fizzbuzz(min, max, num1, num2):
	x = list(range(min, max+1))
	for i in range(len(x)):
		if  x[i] % (num1*num2)==0 : x[i] = 'FizzBuzz'
		elif x[i] % (num1)==0 : x[i] = 'Fizz'
		elif x[i] % (num2)==0 : x[i] = 'Buzz'
	return(x)
	
fizzbuzz(1, 100, 3, 5)
fizzbuzz(-50, 50, 7, 4)

