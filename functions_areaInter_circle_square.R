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