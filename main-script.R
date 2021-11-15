#######################################################
### Codes for analyzing the BCI 50-ha plot data     ### 
### and making figures for the BCI 100 years volume ###
### Author: Camille Piponiot, github.com/cpiponiot  ###
#######################################################

# list of pacakges needed to run this code
req_packages <- c("rdryad", "data.table", "ggplot2", "utils")

# packages that are not yet installed on the computer
ins_packages <-  req_packages[!(req_packages %in% rownames(installed.packages()))]

# install missing packages
if (length(ins_packages) > 0)
  install.packages(ins_packages)

# load all packages
lapply(req_packages, require, character.only = TRUE)

# download BCI data from the Dryad repository
# citation: Condit, Richard et al. (2019), Complete data from the Barro Colorado
# 50-ha plot: 423617 trees, 35 years, Dryad, Dataset,
# https://doi.org/10.15146/5xcp-0d46
dryad_data_path <- rdryad::dryad_download("10.15146/5xcp-0d46")

# unzip files
zip_files <- grep("\\.zip", dryad_data_path$`10.15146/5xcp-0d46`, value = TRUE)
zip_folders <- sapply(zip_files, function(dir) {
  name <- gsub("\\.zip", "", data.table::last(strsplit(dir, '/')[[1]]))
  utils::unzip(dir, exdir = name)
})


