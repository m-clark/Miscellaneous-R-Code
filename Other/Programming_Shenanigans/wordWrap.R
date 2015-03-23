wordWrap = function(text, lineLength){
  require(stringr)
  idx = str_locate(str_sub(text, start = lineLength),' ')[,'start']
  if(is.na(idx)) str_trim(text)
  else paste0(str_trim(str_sub(text, 1, lineLength+idx-1)), '\n', wordWrap(str_sub(text, start=lineLength+idx), lineLength)) 
}


li = 'Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.'

cat(wordWrap(li, 80))
