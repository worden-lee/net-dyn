# source("analyze-batch.r")

# what directory to look in?

load.batch.results <- function(path="batch-data")
{ # list the relevant files under that directory
  files <- dir(path,pattern="network.csv",recursive=TRUE, full.names=TRUE);
  all.rows <- NULL;
  for (f in files)
  { fd <- read.csv(f,header=TRUE);
    if (is.null(all.rows))
      all.rows <- fd
    else
      all.rows <- rbind(all.rows,fd)
  }
  all.rows
}

relative <- function(vec)
{ vec / (sum(vec)/length(vec)) }

stars <- function(p)
{ if(p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else if (p < 0.1) {
    return(".")
  } else {
    return(" ")
  }
}

doplot <- function(xv, yv, xlab, ylab="influence",
                   xlabelpos=0, ylabelpos=0, filename=NULL)
{
  plot(xv,yv,xlab=xlab,ylab=ylab)
  lm <- lm(yv ~ xv)
  lms <- summary(lm)
  abline(lm,col="red")
#  lines=sprintf("r-squared=%g\nslope=%g\nstd.err(slope)=%g\np=%g %s",
#    lms$r.squared, lms$coefficients[2,1],
#    lms$coefficients[2,2],
#    lms$coefficients[2,4], stars(lms$coefficients[2,4]))
# --
#  lines=sprintf("slope=%g\nstd.err(slope)=%g\nr-squared=%g",
#    lms$coefficients[2,1], lms$coefficients[2,2], lms$r.squared)
  lines=sprintf("coefficient=%g %s\nstd. error=%g\nr-squared=%g",
    lms$coefficients[2,1], stars(lms$coefficients[2,4]),
    lms$coefficients[2,2], lms$r.squared)
  text(xlabelpos,ylabelpos, lines, pos=4)
  if (!is.null(filename))
  { dev.copy2eps(file=filename)
    dev.off()
  }
}

