#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Functions used for the BCI 100 years volume     #
# Author: Camille Piponiot, github.com/cpiponiot  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


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

# substitute abnormal changes in DBH ####
## get symmetrical distribution of dbh change by transforming with modulus function
## then calculate mean dbh or agb change and backtransform

substitute_change = function(varD,
                             cut = c(-0.5, 5),
                             lambda = 0.5,
                             value = "D",
                             D = NULL,
                             WD = NULL, 
                             palms = NULL) {
  
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
    dAGB = agb_bci(dbh = D + varD, wd = WD, palms) - 
      agb_bci(dbh = D, wd = WD, palms)
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

ExpDiffAGB = function(d, wd, lambda, mu, sigma, palms = NULL) {
  # AGB > 0 => d + modulus(x, 1/lambda) > 0 => x > -(d^lambda)
  minVar = -d ^ lambda
  pdfdAGB = function(x) {
    dAGB = agb_bci(d + modulus(x, 1 / lambda), wd, palms) - agb_bci(d, wd, palms)
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

kohyama_correction <- function(stock, gain, loss, dT, output = "prod") {
  B0 <- stock
  BS0 <- B0 - loss * dT
  BT <- B0 + (gain - loss) * dT
  if (output == "prod") {
    outp <- (log(BT / BS0) * (BT - B0)) / (dT * log(BT / B0))
  } else if (output == "mort"){
    outp <- (log(B0 / BS0) * (BT - B0)) / (dT * log(BT / B0))
  } else 
    stop("Please provide either 'prod' (production) or 'mort' (mortality) as an output.")
  return(outp)
}

cumsum_naomit <- function (x) cumsum(ifelse(is.na(x), 0, x)) + x*0

resequence <- function (x, dx) {
  if (any(!is.na(x))) {
    # get the first non NA measurement (t0: position in the vector)
    t0 <- which(!is.na(x)) [1]
    # use only diff x measurements after this first measurement, and remove the
    # last one that is always NA (by definition)
    dx <- dx[t0:(length(dx)-1)]
    # new values: same until t0, then add the cumulative sum of dx to the first
    # non NA value
    x <- c(x[1:t0], x[t0] + cumsum_naomit(dx))
  } 
  return(x)
}

# area of the intersection of a circle and a rectangle ####
# adapted from https://stackoverflow.com/questions/622287/area-of-intersection-between-circle-and-rectangle
# integral calculated with www.wolframalpha.com

# param x0  lower
# param x1 

# returns the positive root of intersection of line y = h with circle centered at the origin and radius r
section <- function(h, r = 1) 
{
  # assert(r >= 0); # assume r is positive, leads to some simplifications in the formula below (can factor out r from the square root)
  return(ifelse(h < r, sqrt(r * r - h * h), 0)) # http:#www.wolframalpha.com/input/?i=r+*+sin%28acos%28x+%2F+r%29%29+%3D+h
}

# indefinite integral of circle segment
g <- function(x, h, r = 1) 
{
  return(.5 * (sqrt(1 - x * x / (r * r)) * x * r + r * r * asin(x / r) - 2 * h * x)) # http:#www.wolframalpha.com/input/?i=r+*+sin%28acos%28x+%2F+r%29%29+-+h
}

# area of intersection of an infinitely tall box with left edge at x0, right edge at x1, bottom edge at h and top edge at infinity, with circle centered at the origin with radius r
area_inf <- function(x0, x1, h, r) 
{
  if (x0 > x1)
    seqinr::swap(x0, x1) # this must be sorted otherwise we get negative area
  
  s <- section(h, r)
  return(g(max(-s, min(s, x1)), h, r) - g(max(-s, min(s, x0)), h, r)) # integrate the area
}

# area of the intersection of a finite box with a circle centered at the origin with radius r
area_fin <- function(x0, x1, y0, y1, r) 
{
  if (y0 > y1)
    seqinr::swap(y0, y1)
  
  if (y0 >= 0) {
    # y1 > y0 >=0
    return(area_inf(x0, x1, y0, r) - area_inf(x0, x1, y1, r)) # area of the lower box minus area of the higher box
    
  } else {
    if (y1 < 0) {
      # the box is completely under, just flip it above (opposite y coordinates)
      return(area_inf(x0, x1, -y1, r) - area_inf(x0, x1, -y0, r))
    } else # the box is both above and below, divide it to two boxes 
      return(area_inf(x0, x1, 0, r) - area_inf(x0, x1, -y0, r) + 
               area_inf(x0, x1, 0, r) - area_inf(x0, x1, y1, r))
  }  
}

# area of the intersection of a general box with a general circle
area_circle_rect <- function(x0, x1, y0, y1, cx, cy, r) 
{
  x0 <- x0 - cx; x1 <- x1 - cx;
  y0 <- y0 - cy; y1 <- y1 - cy;
  # get rid of the circle center
  
  return(area_fin(x0, x1, y0, y1, r))
}

### test
test = FALSE
if (test) {
  n=1000
  coords = data.table(id = 1:n,
                      cx = rnorm(n, 0, 10), 
                      cy = rnorm(n, 0, 10),
                      r = rlnorm(n, 0, 1), 
                      x0 = rnorm(n, 0, 10), 
                      y0 = rnorm(n, 0, 10), 
                      x1 = rnorm(n, 0, 10), 
                      y1 = rnorm(n, 0, 10))
  coords[x0 > x1, `:=`(x0=x1, x1=x0)]
  coords[y0 > y1, `:=`(y0=y1, y1=y0)]
  coords[, `:=`(area_i = area_circle_rect(x0, x1, y0, y1, cx, cy, r), 
                area_r = (x1-x0)*(y1-y0), area_c = r^2*pi), .(id)]
  
  coords[area_i<0]
  coords[area_i>area_r]
  coords[area_i>area_c]
  
  Id = 11
  cx = coords$cx[coords$id == Id]
  cy = coords$cy[coords$id == Id]
  r = coords$r[coords$id == Id]
  x0 = coords$x0[coords$id == Id]
  y0 = coords$y0[coords$id == Id]
  x1 = coords$x1[coords$id == Id]
  y1 = coords$y1[coords$id == Id]
  
  
  library(ggplot2)
  ggplot() +
    geom_rect(aes(xmin = x0, ymin = y0, xmax = x1, ymax = y1)) +
    ggforce::geom_circle(aes(x0 = cx, y0 = cy, r = r)) +
    # ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r)) +
    coord_equal() +
    geom_hline(yintercept = 0)
  
  ## create shapefiles 
  library(sf)
  library(sfheaders)
  
  st_area_inter <- function(x0, x1, y0, y1, cx, cy, r) {
    rect <- sf_polygon(data.frame(x = c(x0, x1, x1, x0), y = c(y0, y0, y1, y1)))
    center <- st_point(x = c(cx, cy))
    circle <- st_buffer(x = center, dist = r)
    inter <- st_intersection(rect, circle)
    return(st_area(inter))
  }
  
  coords[area_i>0, area_i_sf := st_area_inter(x0, x1, y0, y1, cx, cy, r), .(id)]
  
  summary(coords[area_i>0, area_i/area_i_sf])
}