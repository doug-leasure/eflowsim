# cleanup
rm(list=ls()); gc(); cat("\014"); try(dev.off(), silent=T)

# working directory
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path)))

# package documentation
devtools::document()

# restart R
rstudioapi::restartSession()

# install package
install.packages(getwd(), repo=NULL, type='source')

# load package
library(eflowsim)
