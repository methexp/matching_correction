# fitting using mptinr
library("MPTinR")
#check.mpt("pudel1.eqn")

# function to read data from an .MPT file into a vector for use with MPTinR
read.mdt <- function(file=file){
  #check which whitespace delimiter
  d <- read.delim(file=file, header = FALSE, sep = "", comment.char = "=", skip=1)
  return(d$V2[1:nrow(d)])
}

# function to create a restrictions file for use with MPTinR
make.restrictions <- function(file=file, restrictions=c("")){
  my.con <- file(file, open = "w")
  for(i in 1:length(restrictions)){
    writeLines(restrictions[i], con = my.con)
  }
  close(my.con)
}

# wrapper function for fit.mpt that takes .mdt files and restrictions in text lists as parameters
fit_mpt <- function(eqnfile, mdtfile, restrictions = c(""), ...){
  data <- read.mdt(file=mdtfile)
  if(restrictions[1] != ""){
    make.restrictions(file="restrictions.tmp", restrictions = restrictions)
    fit <- fit.mpt(data = data, model.filename = eqnfile, restrictions.filename = "restrictions.tmp", ...)
  } else {
    fit <- fit.mpt(data = data, model.filename = eqnfile, ...)
  }
  return(fit)
}

# output formatting function for Gsquare tests
apa.g2 <- function(fitmpt=list(), fitmpt_baseline=list(), df, q, p, delta=FALSE){
  # use values of df, p, q if no fitmpt object is provided
  if(length(fitmpt) == 0){  
    out <- "$"
    if(delta==TRUE) out <- paste0(out, "\\Delta ")
    out <- paste0(out, "G^2_{(df=")
    df_ <- df
    q_ <- q
    p_ <- p
  } else {
  # compute values for G², df, p
    if(length(fitmpt_baseline) == 0){  
      out <- "$G^2_{(df="
      df_ <- fitmpt$goodness.of.fit$df
      q_ <- fitmpt$goodness.of.fit$G.Squared
      p_ <- fitmpt$goodness.of.fit$p.value
    } else {
      out <- "$\\Delta G^2_{(df="
      df_ <- fitmpt$goodness.of.fit$df - fitmpt_baseline$goodness.of.fit$df
      q_ <- fitmpt$goodness.of.fit$G.Squared - fitmpt_baseline$goodness.of.fit$G.Squared
      p_ <- pchisq(q = q_, df = df_, lower.tail = FALSE)
    }
  }
  # format and return results
  out <- paste0(out, df_)
  out <- paste0(out, ")} = ")
  out <- paste0(out, papaja::printnum(q_, digits=2, gt1=TRUE))
  if(p_ < .0005) {
    out <- paste0(out, ", p ")
  } else {
    out <- paste0(out, ", p = ")
  }
  out <- paste0(out, papaja::printp(p_))
  out <- paste0(out, "$")
  return(out)      
}

make.resulttable <- function(fitmpt=list(), conditions=0, condnames, parorder, sep="_"){
  if(length(fitmpt) != 0){
    # get parameter estimates
    pars <- fitmpt$parameters
    total <- length(pars$estimates)
    percond <- total/conditions
    est <- low <- upp <- out <- matrix(NA, nrow=percond, ncol=conditions)
    parnames <- 1:percond
    for(i in 1:percond){
      parnames[i] <- strsplit(rownames(pars)[(i-1)*4+1], sep)[[1]][1]
      for(j in 1:conditions){
        est[i,j] <- pars$estimates[(i-1)*4+j]
        low[i,j] <- max(pars$lower.conf[(i-1)*4+j], 0)
        upp[i,j] <- min(pars$upper.conf[(i-1)*4+j], 1)
        out[i,j] <- paste0(round(est[i,j],2), " (", round(low[i,j],2), ", ", round(upp[i,j],2), ")")
        if(is.na(low[i,j])) out[i,j] <- paste0(round(est[i,j],2))
      }
    }
    out <- as.data.frame(out, stringsAsFactors=FALSE)
    names(out) <- condnames
    rownames(out) <- parnames
    out <- out[parorder,]
    return(out)
  }
}
