# Function to check whether package is installed
is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
} 

# Check if sleuth and ggplot2 exist
code <- 0
if (!is.installed('ggplot2')){
  code <- code + 1
}

if (!is.installed('sleuth')){
  code <- code + 2
}

quit(status = code)