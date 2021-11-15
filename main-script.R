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

# source internal functions 
# source("functions.R")

### Get BCI 50-ha plot data ####

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

# load species table
load(grep("spp", dryad_data_path$`10.15146/5xcp-0d46`, value = TRUE))

# load stem censuses
# 1. list census files
bci_stem <- list.files("bci.stem")

# 2. load census files as list
census_list <- lapply(bci_stem, function(name) {
  load(paste0("bci.stem/", name)); get(strsplit(name, "\\.rda")[[1]][1])
})

# 3. add census number information as list names
names(census_list) <- data.table::tstrsplit(bci_stem, "\\.")[[2]]

# 4. compile as one data table
df_census <- data.table::rbindlist(census_list, fill = TRUE, idcol = "censusID")









