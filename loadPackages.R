# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# loadPackages.R
# 
# Checks if required packages are installed. If not, installs them.
# Then loads all required packages.
# 
# 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

cranPackages <- c("coda",
                  "tidyverse",
                  "reshape2",
                  "ggforce",
                  "ggplot2",
                  "GGally",
                  "TMB",
                  "raster",
                  "grid",
                  "RColorBrewer",
                  "HapEstXXR",
                  "parallel",
                  "stringr",
                  "wesanderson",
                  "scales",
                  "beepr",
                  "tmbstan",
                  "here",
                  "bookdown" )

for( pkg in cranPackages )
  while(!require(pkg, character.only = TRUE) )
    install.packages( pkg, repos = "https://mirror.its.sfu.ca/mirror/CRAN/" )


githubPackages <- c(ggsidekick = "seananderson/ggsidekick")

for( pkgIdx in 1:length(githubPackages) )
  while(!require(names(githubPackages)[pkgIdx], character.only = TRUE))
    devtools::install_github(githubPackages[pkgIdx])
