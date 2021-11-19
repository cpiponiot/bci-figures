#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Functions used for the BCI 100 years volume     #
# Author: Camille Piponiot, github.com/cpiponiot  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# estimate agb values at BCI with allometric equations from Chave et al. 2014
# (method = "chave14") or Chave et al., 2005 (method = "chave05"); the height
# equation from Martinez Cano et al 2019 is used when use_height_allom = TRUE
agb_bci <- function(dbh,
                    wd,
                    method = "chave14",
                    use_height_allom = FALSE) {
  # with generic height allometry from Martinez Cano et al 2019
  if (use_height_allom) {
    # Chave et al 2005 - moist forests, with height (in Mg)
    h <- 58.0 * dbh ^ 0.73 / (21.8 + dbh ^ 0.73)
    if (method == "chave05") {
      agb <- 0.0509 * wd * dbh ^ 2 * h / 1000
    }
    if (method == "chave14") {
      # Chave et al. 2014, equation 4 with the BIOMASS package 
      agb <- (0.0673 * (wd * h * dbh^2)^0.976)/1000
    }
  } else {
    # without any height information
    # Chave et al 2005 - moist forests, without height (in Mg)
    if (method == "chave05") {
      agb <- wd * exp(-1.499 + 2.148 * log(dbh) + 
                        0.207 * log(dbh) ^ 2 - 0.0281 * log(dbh) ^ 3) /1000
    }
    if (method == "chave14") {
      # Chave et al. 2014, equation 7 with the BIOMASS package (transform into kg)
      E <- 0.05176398
      agb <- exp(-2.023977 - 0.89563505 * E + 0.92023559 * 
                   log(wd) + 2.79495823 * log(dbh) - 0.04606298 * (log(dbh)^2))/1000
    }
  }
  
  return(agb)
}


# function to interpolate missing dbh value ####

# param dbh: all dbhs measurements of one tree or stem
# param year: the years corresponding to those measurements
# param DFstatus: the status of the tree or stem for all those measurements

# returns a dbh vector with all missing values interpolated (only when there is
# at least one dbh measurement before and after the missing dbh)

interpolate_missing <- function(dbh, year, DFstatus) {
  if (any(DFstatus == "missing" & !is.na(DFstatus))) { 
    # order dbh and year vector in chronological order
    year_ord <- year[order(year)]
    dbh_ord <- dbh[order(year)]
    
    # years of all missing values 
    for (y0 in year[DFstatus == "missing" & !is.na(DFstatus)]) {
      # last (non NA) dbh measurement prior to the missing measurement
      dbh1 <- data.table::last(dbh_ord[!is.na(dbh_ord) & year_ord < y0])
      # first (non NA) dbh measurement after the missing measurement
      dbh2 <- data.table::first(dbh_ord[!is.na(dbh_ord) & year_ord > y0])
      
      # check that dbh1 and dbh2 have values (otherwise: no data to interpolate from)
      if (length(dbh1) > 0 & length(dbh2) > 0) {
        # last year with (non NA) dbh measurement prior to the missing measurement
        year1 <- max(year_ord[!is.na(dbh_ord) & year_ord < y0])
        # first year with (non NA) dbh measurement after the missing measurement
        year2 <- min(year_ord[!is.na(dbh_ord) & year_ord > y0])
        
        # slope of the regression between the measurements before and after the
        # missing value
        slope <- (dbh2 - dbh1) / (year2 - year1)
        # interpolation:
        dbh_ord[year_ord == y0] <- slope * (y0 - year1) + dbh1
      } 
    }
    
    return(dbh_ord[order(order(year))])
    
  } else {
    return(dbh)
  }
}

# substitute abnormal hcanges in DBH ####
## get symmetrical distribution of dbh change by transforming with modulus function
## then calculate mean dbh or agb change and backtransform

substitute_change = function(varD,
                             cut = c(-0.5, 5),
                             lambda = 0.5,
                             value = "D",
                             D = NULL,
                             WD = NULL) {
  
  keep_values = which(varD > cut[1] & varD < cut[2] & !is.na(varD))
  change_values = which((varD <= cut[1] | varD >= cut[2]) & !is.na(varD))
  transf_values = modulus(varD[keep_values], lambda)
  mu = mean(transf_values)
  sigma = sd(transf_values)
  
  if (value == "D") {
    diffModD = function(x) modulus(x, 1 / lambda) * dnorm(x, mu, sigma)
    varD[change_values] = integrate(diffModD,-Inf, Inf)$value
    return(varD)
  }
  
  if (value == "AGB") {
    dAGB = agb_bci(dbh = D + varD, wd = WD) - agb_bci(dbh = D, wd = WD)
    if (length(change_values) == 1) {
      dAGB[change_values] = ExpDiffAGB(
        d = D[change_values],
        wd = WD[change_values],
        lambda = lambda,
        mu = mu,
        sigma = sigma
      )
    } else if (length(change_values) > 1) {
      dAGB[change_values] = apply(cbind(D, WD)[change_values, ], 1, function(x) {
        ExpDiffAGB(
          d = x[1],
          wd = x[2],
          lambda = lambda,
          mu = mu,
          sigma = sigma
        )
      })
    }
    return(dAGB)
  }
}

modulus = function(d, lambda = 0.4) {
  return(sign(d) * abs(d) ^ lambda)
}

ExpDiffAGB = function(d, wd, lambda, mu, sigma) {
  # AGB > 0 => d + modulus(x, 1/lambda) > 0 => x > -(d^lambda)
  minVar = -d ^ lambda
  pdfdAGB = function(x) {
    dAGB = agb_bci(d + modulus(x, 1 / lambda), wd) - agb_bci(d, wd)
    dens = truncnorm::dtruncnorm(x, mean = mu, sd = sigma, a = minVar)
    return(dAGB * dens)
  }
  return(integrate(pdfdAGB, lower = minVar, upper = Inf)$value)
}

## modulus(Ddbh) not exactly normal: bimodal distribution (looks more like N(0,sd) + logN(mu, sd2))

# test
# x = fgeo_data[site=="bci"&dbhc>50]
# x[, dvar := c(diff(dbhc) / diff(year), NA)]
# # last stem measurement: dvar = NA
# x[, lastMeas := c(stemid[-1] != stemid[-nrow(x)], TRUE)]
# x[(lastMeas), dvar := NA]
# varD = x$dvar
# D = x$dbhc
# WD = x$wsg
# E = 1.5
# hom_change = x$dHOM
#


# Kohyama correction function ####
# instanteneous biomass (or other) fluxes as recommended by Kohyama 2019 (eq 1-2 in Table 1)

kohyama_correction <- function(stock, influx, outflux) {
  dt_sub = dt[variable %in% vars, .(value = sum(value * weight) / sum(weight),
                                    weight = sum(weight)),
              .(dT, year, variable, group, site)]
  dt_new = dcast(dt_sub, year + dT + site + weight + group ~ variable)
  dt_new$B0 = dt_new[, vars[1], with = FALSE]
  dt_new$gain = dt_new[, vars[2], with = FALSE]
  dt_new$loss = dt_new[, vars[3], with = FALSE]
  dt_new[, BS0 := B0 - loss * dT]
  dt_new[, BT := B0 + (gain - loss) * dT]
  dt_new[, vars[2] := (log(BT / BS0) * (BT - B0)) / (dT * log(BT / B0))]
  dt_new[, vars[3] := (log(B0 / BS0) * (BT - B0)) / (dT * log(BT / B0))]
  
  dt_new = melt(dt_new,
                id.vars = c("site", "weight", "group"),
                measure.vars = vars)
  dt_corr = rbind(dt[!variable %in% vars, colnames(dt_new), with = FALSE], dt_new)
  dt_corr = subset(dt_corr,!is.na(value))
  return(dt_corr)
}
