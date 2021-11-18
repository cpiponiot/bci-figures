#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Codes for analyzing the BCI 50-ha plot data and #
# making figures for the BCI 100 years volume     #
# Author: Camille Piponiot, github.com/cpiponiot  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# list of package dependencies
req_packages <- c("rdryad", "data.table", "ggplot2", "utils", "BIOMASS")

# packages that are not yet installed on the computer
ins_packages <-  req_packages[!(req_packages %in% rownames(installed.packages()))]

# install missing packages
if (length(ins_packages) > 0) 
  install.packages(ins_packages)

# load packages
library(ggplot2)

# source internal functions 
source("functions.R")

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
df_stem <- data.table::rbindlist(census_list, fill = TRUE, idcol = "censusID")

# add species information
df_stem <- merge(df_stem, bci.spptable, by = "sp", all.x = TRUE) 

# get census year (median value of all years of measurement for one census)
df_stem[, year := as.numeric(data.table::tstrsplit(ExactDate, "-")[[1]])]
df_stem[, year := median(year, na.rm = TRUE), .(censusID)]

### estimate individual aboveground biomass using different methods ####

# no correction

# Chave et al 2014 allometric equation, no height information
df_stem[, chave14 := agb_bci(dbh = dbh/10, wd = wsg, method = "chave14")]

# Chave et al 2014 allometric equation, tree height from Martinez Cano et al., 2019
df_stem[, chave14_h := agb_bci(dbh = dbh/10, wd = wsg, method = "chave14", use_height_allom = TRUE)]

# Chave et al 2005 allometric equation, no height information
df_stem[, chave05 := agb_bci(dbh = dbh/10, wd = wsg, method = "chave05")]

# Chave et al 2005 allometric equation, tree height from Martinez Cano et al., 2019
df_stem[, chave05_h := agb_bci(dbh = dbh/10, wd = wsg, method = "chave05", use_height_allom = TRUE)]

# corr1: taper correction
# from Cushman et al., 2021, using WSG 
df_stem[, b := 0.151 - 0.025 * log(dbh/10) - 0.02 * log(hom) - 0.021 * log(wsg)]
df_stem[!is.na(hom), dbh_t := dbh * exp(b * (hom - 1.3))]
df_stem[, chave14_t := agb_bci(dbh = dbh_t/10, wd = wsg, method = "chave14", use_height_allom = TRUE)]

# corr2: interpolate missing DBHs
df_stem[, chave14_ti := agb_bci(dbh = dbh_ti/10, wd = wsg, method = "chave14", use_height_allom = TRUE)]


# corr3: replace DHB or AGB growth
missing_trees <- unique(df_stem$treeID[df_stem$DFstatus=="missing"&!is.na(df_stem$DFstatus)])
df_stem[df_stem$treeID %in% missing_trees, any(!is.na(dbh) & year < min(year[missing])), .(stemID)]

# plot-level AGB ####

# melt to long format for different methods
df_stem_melt <-
  data.table::melt(
    df_stem,
    id.vars = c("stemID", "treeID", "year", "quadrat", "DFstatus"), 
    measure.vars = grep("chave", colnames(df_stem)),
    variable.name = "allometry",
    value.name = "agb"
  )

# replace NAs with 0 in agb values
df_stem_melt[is.na(agb), agb := 0]

# aggregate agb values, in Mg/ha (divide by area = 50 ha)
df_plot <- df_stem_melt[DFstatus == "alive", .(agb = sum(agb) / 50), .(allometry, year)]


# corr4: kohyama correction 

#### Figure 1 - AGB estimation and potential sources of uncertainty ####

# create data table for illustrating different allometric equations at the
# individual tree level
dfallom <- data.table::data.table(expand.grid(dbh = 1:200, wd = c(0.4, 0.8)))

# Chave et al 2014 allometric equation, no height information
dfallom[, chave14 := agb_bci(dbh = dbh, wd = wd, method = "chave14")]

# Chave et al 2014 allometric equation, tree height from Martinez Cano et al., 2019
dfallom[, chave14_h := agb_bci(dbh = dbh, wd = wd, method = "chave14", use_height_allom = TRUE)]

# Chave et al 2005 allometric equation, no height information
dfallom[, chave05 := agb_bci(dbh = dbh, wd = wd, method = "chave05")]

# Chave et al 2005 allometric equation, tree height from Martinez Cano et al., 2019
dfallom[, chave05_h := agb_bci(dbh = dbh, wd = wd, method = "chave05", use_height_allom = TRUE)]

dfallom <- data.table::melt(dfallom, measure.vars = grep("chave", colnames(dfallom)), 
                            variable.name = "allometry", value.name = "agb")

ggplot(dfallom, aes(x = dbh, y = agb, colour = allometry, linetype = as.factor(wd))) +
  geom_line() +
  theme_classic()

ggplot(subset(df_plot, year >= 1985)) + 
  geom_line(aes(x = year, y = agb, color = allometry))+
  theme_classic()

