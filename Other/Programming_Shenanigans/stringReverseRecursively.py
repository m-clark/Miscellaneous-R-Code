def reverseString(aStr):
  if len(aStr) == 1:
    return aStr
  return  aStr[-1] + reverseString(aStr[:-1])

import string

string.ascii_letters
reverseString(string.ascii_letters)
