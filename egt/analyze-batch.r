# source("analyze-batch.r")

# what directory to look in?

load.batch.results <- function(path="batch-data",target="network.csv")
{ # list the relevant files under that directory
  cat("expanding list of files in (",path,")\n");
  files <- dir(path,pattern=target,recursive=TRUE, full.names=TRUE);
  cat("found ",length(files)," files\n");
  #if (length(files) > 1000)
  #{ files <- files[1:1000]; }
  all.rows <- NULL;
  for (f in files)
  { cat(f,"\n");
    fd <- read.csv(f,header=TRUE);
    rep <- sub("^.*out.([^/]*)/.*$","\\1",f);
    fd <- cbind(data.frame(rep=rep),fd);
    if (target == "nodes.csv")  # hell of inefficient
    { netf <- sub("nodes.csv","network.csv",f);
      net <- read.csv(netf,header=TRUE);
      net <- cbind(data.frame(rep=rep),net);
      fd <- merge(net,fd,by=c("rep"),suffixes=c(".net",".node"));
    }
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
                   xlabelpos=0, ylabelpos=0, lm=NULL, filename=NULL)
{
  plot(xv,yv,xlab=xlab,ylab=ylab)
  if (is.null(lm))
  { lm <- lm(yv ~ xv) }
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

