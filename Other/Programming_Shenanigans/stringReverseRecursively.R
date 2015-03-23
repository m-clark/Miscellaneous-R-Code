# reverse a string recursively

revString = function(string){
  require(stringr)
  if (nchar(string)==1) string
  else paste0(str_sub(string, -1), revString(str_sub(string, end=-2)))
}


string = paste0(letters, collapse = '')

revString(string)

# base R approach
revString = function(string){
  if (nchar(string)==1) string
  else paste0(substr(string, nchar(string), nchar(string)), 
              revString(substr(string, 1, nchar(string)-1)))
}


