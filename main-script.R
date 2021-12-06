#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Codes for analyzing the BCI 50-ha plot data and #
# making figures for the BCI 100 years volume     #
# Author: Camille Piponiot, github.com/cpiponiot  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


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

# remove large strangler figs (> 50 cm DBH) from the following species: Ficus
# costaricana, Ficus obtusifolia, Ficus popenoei, and Ficus trigonata (see
# Rutishauser et al 2020)

str_figs <- subset(df_stem, (dbh > 50 & Genus == "Ficus" & Species %in% c("costaricana", "obtusifolia", "popenoei", "trigonata")))
length(unique(str_figs$quadrat))
length(unique(str_figs$treeID))
str_figs[, .(nquadrat = length(unique(quadrat)), nfigs = length(unique(treeID))), .(census_year)]

df_stem <- subset(df_stem, ! treeID %in% str_figs$treeID)

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
maxD <- ceiling(max(df_stem$dbh_ti, na.rm = TRUE))
df_stem[, size := cut(dbh_ti, c(1, 10, 20, 30, 50, 100, maxD), include.lowest = TRUE)]

# functional groups, based on Ruger et al 2020, Data S1
# download data
if (!file.exists("ruger_data_s1.xlsx"))
  utils::download.file(url = "https://www.science.org/doi/suppl/10.1126/science.aaz4797/suppl_file/aaz4797_ruger_data_s1.xlsx", 
                       destfile = "ruger_data_s1.xlsx")
# read downloaded data
ruger_data <- readxl::read_xlsx("ruger_data_s1.xlsx", sheet = "Data")
ruger_data$sp <- tolower(ruger_data$sp)

# get PFT definition from metadata
ruger_levels <- strsplit("1=slow, 2=fast, 3=LLP, 4=SLB, 5=intermediate", "=|, ")[[1]]
ruger_levels <- data.frame(matrix(ruger_levels, ncol = 2, byrow = TRUE))
colnames(ruger_levels) = c("PFT_2axes", "PFT")
ruger_data <- merge(ruger_data, ruger_levels)

# add PFTs to df_stem
df_stem <- merge(df_stem, ruger_data[, c("sp", "PFT")], all.x = TRUE)

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
    id.vars = c("stemID", "treeID", "census_year", "quadrat", "DFstatus", "size", "PFT"), 
    measure.vars = grep("chave|Dagb", colnames(df_stem)),
    variable.name = "method",
    value.name = "agb"
  )

# replace NAs with 0 in agb values
df_stem_melt[is.na(agb), agb := 0]

# aggregate agb values ####
# agb in Mg/ha (divide by area = 50 ha)
# 1. by year ####
df_plot <- df_stem_melt[, .(value = sum(agb) / 50), .(method, year = census_year)]
df_plot[, variable := c("agb", "awp")[grepl("Dagb", method)+1]]

levels(df_plot$method) <- list("Chave 2014" = "chave14", 
                               "Chave 2014 + height allom" = "chave14_h", 
                               "Chave 2005" = "chave05", 
                               "Chave 2005 + height allom" = "chave05_h", 
                               "Taper correction" = "chave14_t", 
                               "Taper correction\n+ missing stems" = "chave14_ti", 
                               "No correction" = "Dagb", 
                               "Substitution" = "Dagb_s")

# corr4: kohyama correction 
# > need to estimate mortality
dfK <-
  df_plot[year > 1982 &
            method %in% c("Substitution", "Taper correction\n+ missing stems")]
dfK <- data.table::dcast(dfK, year ~ variable)
dfK[order(year), `:=` (dT = c(diff(year), NA), 
                       awm = awp - c(diff(agb)/diff(year), NA))]
dfK[, awp_k := kohyama_correction(agb, awp, awm, dT)]
dfK <- data.table::data.table(method = "kohyama", year = dfK$year, 
                              value = dfK$awp_k, variable = "awp")
df_plot <- rbind(df_plot, dfK)

# 2. by quadrat ####
df_quadrat <- df_stem_melt[, .(value = sum(agb) / 0.2^2), 
                           .(method, year = census_year, quadrat)]
df_quadrat[, variable := c("agb", "awp")[grepl("Dagb", method)+1]]

# kohyama correction 
# > need to estimate mortality
df_quadrat <- df_quadrat[method %in% c("chave14_ti", "Dagb_s")]
df_quadrat <- data.table::dcast(df_quadrat, year + quadrat ~ variable)
df_quadrat[order(year), `:=` (dT = c(diff(year), NA), 
                              awm = awp - c(diff(agb)/diff(year), NA)), .(quadrat)]
df_quadrat[, awp_k := kohyama_correction(agb, awp, awm, dT)]
df_quadrat[, `:=`(dT = NULL, awp = awp_k, awm = NULL, awp_k = NULL)]

# mean across all years
df_quadrat <- df_quadrat[year > 1982 & year < 2015, .(
  agb = mean(agb), awp = mean(awp)
), .(quadrat)]

# add quadrat coordinates
df_quadrat$X <- as.numeric(substr(df_quadrat$quadrat, 1, 2))*20 + 10
df_quadrat$Y <- as.numeric(substr(df_quadrat$quadrat, 3, 4))*20 + 10
df_quadrat <- subset(df_quadrat, !is.na(quadrat) & quadrat != "")

# 3. by size class ####
df_size <- df_stem_melt[, .(value = sum(agb) / 50), 
                        .(method, year = census_year, size)]
df_size[, variable := c("agb", "awp")[grepl("Dagb", method)+1]]

# kohyama correction 
# > need to estimate mortality
df_size <- df_size[method %in% c("chave14_ti", "Dagb_s")]
df_size <- data.table::dcast(df_size, year + size ~ variable)
df_size[order(year), `:=` (
  dT = c(diff(year), NA), 
  awm = awp - c(diff(agb)/diff(year), NA)
), .(size)]
df_size[, awp_k := kohyama_correction(agb, awp, awm, dT)]
df_size[, `:=`(dT = NULL, awp = awp_k, awm = NULL, awp_k = NULL)]

# mean across all years
df_size <- df_size[!is.na(size) & year > 1982 & year < 2015, .(
  agb = mean(agb), awp = mean(awp)
), .(size)]

# 4. by functional group ####
# xxx create function for this
df_pft <- df_stem_melt[, .(value = sum(agb) / 50), 
                        .(method, year = census_year, PFT)]
df_pft[, variable := c("agb", "awp")[grepl("Dagb", method)+1]]

# kohyama correction 
# > xxx need to estimate mortality
df_pft <- df_pft[method %in% c("chave14_ti", "Dagb_s")]
df_pft <- data.table::dcast(df_pft, year + PFT ~ variable)
df_pft[order(year), `:=` (
  dT = c(diff(year), NA), 
  awm = awp - c(diff(agb)/diff(year), NA)
), .(PFT)]

df_pft[, awp_k := kohyama_correction(agb, awp, awm, dT)]
df_pft[, `:=`(dT = NULL, awp = awp_k, awm = NULL, awp_k = NULL)]

# mean across all years
df_pft <- df_pft[year > 1982 & year < 2015, .(
  agb = mean(agb), awp = mean(awp)
), .(PFT)]

df_929 <- subset(df_stem, stemID==929)
save(df_plot, df_size, df_pft, df_929, file = "data/data-main-script.rda")

# 
# # figure 4 ####
# fig4a_map_agb <- ggplot(df_quadrat, aes(x = X, y = Y, fill = agb)) + 
#   geom_raster() + 
#   labs(fill = "AGB\n(Mg/ha)") +
#   coord_equal() + 
#   expand_limits(x = c(0, 1000), y = c(0, 500)) +   
#   scale_x_continuous(expand = c(0, 0)) + 
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_gradientn(colours = rev(heat.colors(10))) + 
#   theme_bw()
# fig4b_map_awp <- ggplot(df_quadrat, aes(x = X, y = Y, fill = awp)) + 
#   geom_raster() + 
#   labs(fill = "AWP\n(Mg/ha/yr)") +
#   coord_equal() + 
#   expand_limits(x = c(0, 1000), y = c(0, 500)) +   
#   scale_x_continuous(expand = c(0, 0)) + 
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_gradientn(colours = rev(heat.colors(10))) + 
#   theme_bw()
# ggpubr::ggarrange(fig4a_map_agb, fig4b_map_awp, ncol = 1, labels = "auto")
# ggsave("fig4_map.pdf", height = 4, width = 4)
