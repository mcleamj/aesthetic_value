variablecol <- function(colvar, col, clim) {
  ncol <- length(col)
  colvar[colvar < min(clim)] <- NA
  colvar[colvar > max(clim)] <- NA
  rn <- clim[2] - clim[1]
  ifelse (rn != 0, Col <- col[1 + trunc((colvar - clim[1])/rn *
                                          (ncol - 1)+1e-15)], Col <- rep(col[1], ncol))               
  return(Col)
}