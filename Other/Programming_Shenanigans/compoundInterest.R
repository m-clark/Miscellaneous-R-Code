# a (yearly) compound interest calculator to demonstrate recursion in a simple fashion

compoundInterest <- function(principal, years, rate, yearlybump) {
  if(years==0) return(principal)
  
  if(years>0) compoundInterest(principal=principal*(1+rate) + yearlybump, years=years-1, rate=rate, yearlybump=yearlybump)
}

compoundInterest(principal=10000, years=10, rate=.05, yearlybump=1000)

