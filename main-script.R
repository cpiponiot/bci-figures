#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Codes for analyzing the BCI 50-ha plot data and #
# making figures for the BCI 100 years volume     #
# Author: Camille Piponiot, github.com/cpiponiot  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# list of package dependencies
req_packages <- c("rdryad", "data.table", "ggplot2", "utils", "truncnorm", "ggpubr")

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
df_stem[, census_year := median(year, na.rm = TRUE), .(censusID)]

# order df_stem by treeID, stemID and year
data.table::setorder(df_stem, treeID, stemID, year)

# dbh in cm 
df_stem[, dbh := dbh/10]

# remove measures under 1 cm 
df_stem <- subset(df_stem, dbh >= 1 | is.na(dbh))

### estimate individual aboveground biomass using different methods ####

# no correction

# Chave et al 2014 allometric equation, no height information
df_stem[, chave14 := agb_bci(dbh = dbh, wd = wsg, method = "chave14")]

# Chave et al 2014 allometric equation, tree height from Martinez Cano et al., 2019
df_stem[, chave14_h := agb_bci(dbh = dbh, wd = wsg, method = "chave14", use_height_allom = TRUE)]

# Chave et al 2005 allometric equation, no height information
df_stem[, chave05 := agb_bci(dbh = dbh, wd = wsg, method = "chave05")]

# Chave et al 2005 allometric equation, tree height from Martinez Cano et al., 2019
df_stem[, chave05_h := agb_bci(dbh = dbh, wd = wsg, method = "chave05", use_height_allom = TRUE)]

# corr1: taper correction
# from Cushman et al., 2021, using WSG 
df_stem[, b := 0.151 - 0.025 * log(dbh) - 0.02 * log(hom) - 0.021 * log(wsg)]
df_stem[!is.na(hom), dbh_t := dbh * exp(b * (hom - 1.3))]
df_stem[, chave14_t := agb_bci(dbh = dbh_t, wd = wsg)]

# corr2: interpolate missing DBHs
df_stem[, dbh_ti := interpolate_missing(dbh_t, year, DFstatus), .(stemID)]
df_stem[, chave14_ti := agb_bci(dbh = dbh_ti, wd = wsg)]

# corr3: replace DHB or AGB growth
# estimate dbh variation between two censuses, in cm/yr
df_stem[, Ddbh := c(diff(dbh_ti)/diff(year), NA), .(stemID)]
df_stem[, Dagb := c(diff(chave14)/diff(year), NA), .(stemID)]
# size groups / other grouping factors?
df_stem[, size := cut(dbh_ti, c(0.9, 10, 20, 30, 50, 100, 1000))]
## change abnormal dbh changes, grouping by size
df_stem[!is.na(Ddbh), `:=`(
  Ddbh_s = substitute_change(varD = Ddbh, cut = c(-0.5, 5)),
  Dagb_s = substitute_change(D = dbh_ti, varD = Ddbh, WD = wsg, 
                             value = "AGB", cut = c(-0.5, 5))),
  size]


# plot-level AGB and AWP ####

# melt to long format with a method column = allometries and methods used, and
# an 'agb' column = agb estimations under different methods
df_stem_melt <-
  data.table::melt(
    df_stem,
    id.vars = c("stemID", "treeID", "census_year", "quadrat", "DFstatus"), 
    measure.vars = grep("chave|Dagb", colnames(df_stem)),
    variable.name = "method",
    value.name = "agb"
  )

# replace NAs with 0 in agb values
df_stem_melt[is.na(agb), agb := 0]

# aggregate agb values, in Mg/ha (divide by area = 50 ha)
df_plot <- df_stem_melt[, .(value = sum(agb) / 50), .(method, year = census_year)]
df_plot[, variable := c("agb", "awp")[grepl("Dagb", method)+1]]

levels(df_plot$method) <- list("Chave 2014" = "chave14", 
                               "Chave 2014 + height allom" = "chave14_h", 
                               "Chave 2005" = "chave05", 
                               "Chave 2005 + height allom" = "chave05_h", 
                               "Taper correction" = "chave14_t", 
                               "Taper correction + missing stems" = "chave14_ti", 
                               "No correction" = "Dagb", 
                               "Substitution" = "Dagb_s")

# corr4: kohyama correction 
# df_plot[, kohyama_correction()

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
                            variable.name = "method", value.name = "agb")
levels(dfallom$method) <- list("Chave 2014" = "chave14", 
                               "Chave 2014 + height allom" = "chave14_h", 
                               "Chave 2005" = "chave05", 
                               "Chave 2005 + height allom" = "chave05_h")
fig1a_allom_ind <-
  ggplot(dfallom, aes(
    x = dbh,
    y = agb,
    colour = method,
    linetype = as.factor(wd)
  )) +
  geom_line() +
  theme_classic() + 
  expand_limits(x = 0, y = 0) +   
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Stem diameter (cm)",
       y = "Estimated aboveground biomass of trees (Mg)",
       colour = "Allometry used",
       lty = "Wood density") 
fig1_leg <- ggpubr::as_ggplot(ggpubr::get_legend(fig1a_allom_ind))

fig1b_allom_plot <- ggplot(subset(df_plot, year >= 1985 & !grepl("tion", method))) +
  geom_line(aes(x = year, y = value, color = method)) +
  theme_classic() +
  # expand_limits(y = 0) +   
  # scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Census year",
       y = "Estimated plot aboveground biomass (Mg/ha)",
       colour = "Allometry used") +
  theme(legend.position = "none")

ggpubr::ggarrange(fig1a_allom_ind + theme(legend.position = "none"), 
                  fig1b_allom_plot, fig1_leg, 
                  labels = c("a", "b", ""), ncol = 3, widths = c(2,2,1))
ggsave("fig1_allom.pdf", height = 4, width = 10)

dfcorr <- subset(df_plot, year >= 1985 &
                   (grepl("corr", method) | method == "Chave 2014") & 
                   variable == "agb")
levels(dfcorr$method)[levels(dfcorr$method) == "Chave 2014"] <- "No correction"

fig2b_corr_agb <-
  ggplot(dfcorr) +
  geom_line(aes(x = year, y = value, color = method)) +
  theme_classic() +
  labs(x = "Census year",
       y = "Estimated plot aboveground biomass (Mg/ha)",
       colour = "Correction applied") 

dfcorr <- subset(df_plot, year >= 1985 &
                   (grepl("corr", method) | method == "Chave 2014") & 
                   variable == "agb")
levels(dfcorr$method)[levels(dfcorr$method) == "Chave 2014"] <- "No correction"

fig2c_corr_awp <-
  ggplot(subset(df_plot, year >= 1985 & variable == "awp")) +
  geom_line(aes(x = year, y = value, color = method)) +
  theme_classic() +
  labs(x = "Census year",
       y = "Estimated plot aboveground woody \nproductivity (Mg/ha/yr)",
       colour = "Correction applied") 
