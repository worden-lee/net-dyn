# source("load-nodes.r")

relative <- function(vec)
{ vec / (sum(vec)/length(vec)) }

#nodes4 <- read.table("batch-data/network-influence-power/out/out/nodes.csv.~4~",sep=",",header=TRUE)
#nodes4$dp.dr.i. <- relative(nodes4$dp.dr.i.)

nodes.15 <- read.table("batch-data/network-influence-power/out.15/nodes.csv",sep=",",header=TRUE)
nodes.15$dp.dr.i. <- relative(nodes.15$dp.dr.i.)
nodes.26 <- read.table("batch-data/network-influence-power/out.26/nodes.csv",sep=",",header=TRUE)
nodes.26$dp.dr.i. <- relative(nodes.26$dp.dr.i.)
nodes.27 <- read.table("batch-data/network-influence-power/out.27/nodes.csv",sep=",",header=TRUE)
nodes.27$dp.dr.i. <- relative(nodes.27$dp.dr.i.)
nodes.30 <- read.table("batch-data/network-influence-power/out.30/nodes.csv",sep=",",header=TRUE)
nodes.30$dp.dr.i. <- relative(nodes.30$dp.dr.i.)
nodes.36 <- read.table("batch-data/network-influence-power/out.36/nodes.csv",sep=",",header=TRUE)
nodes.36$dp.dr.i. <- relative(nodes.36$dp.dr.i.)
nodes.45 <- read.table("batch-data/network-influence-power/out.45/nodes.csv",sep=",",header=TRUE)
nodes.45$dp.dr.i. <- relative(nodes.45$dp.dr.i.)
nodes.46 <- read.table("batch-data/network-influence-power/out.46/nodes.csv",sep=",",header=TRUE)
nodes.46$dp.dr.i. <- relative(nodes.46$dp.dr.i.)
nodes.47 <- read.table("batch-data/network-influence-power/out.47/nodes.csv",sep=",",header=TRUE)
nodes.47$dp.dr.i. <- relative(nodes.47$dp.dr.i.)
nodes.67 <- read.table("batch-data/network-influence-power/out.67/nodes.csv",sep=",",header=TRUE)
nodes.67$dp.dr.i. <- relative(nodes.67$dp.dr.i.)
nodes.71 <- read.table("batch-data/network-influence-power/out.71/nodes.csv",sep=",",header=TRUE)
nodes.71$dp.dr.i. <- relative(nodes.71$dp.dr.i.)
nodes.72 <- read.table("batch-data/network-influence-power/out.72/nodes.csv",sep=",",header=TRUE)
nodes.72$dp.dr.i. <- relative(nodes.72$dp.dr.i.)
nodes.79 <- read.table("batch-data/network-influence-power/out.79/nodes.csv",sep=",",header=TRUE)
nodes.79$dp.dr.i. <- relative(nodes.79$dp.dr.i.)

betweenness.15 <- read.table("batch-data/network-influence-power/out.15/betweenness.csv",sep=",",header=TRUE)
betweenness.26 <- read.table("batch-data/network-influence-power/out.26/betweenness.csv",sep=",",header=TRUE)
betweenness.27 <- read.table("batch-data/network-influence-power/out.27/betweenness.csv",sep=",",header=TRUE)
betweenness.30 <- read.table("batch-data/network-influence-power/out.30/betweenness.csv",sep=",",header=TRUE)
betweenness.36 <- read.table("batch-data/network-influence-power/out.36/betweenness.csv",sep=",",header=TRUE)
betweenness.45 <- read.table("batch-data/network-influence-power/out.45/betweenness.csv",sep=",",header=TRUE)
betweenness.46 <- read.table("batch-data/network-influence-power/out.46/betweenness.csv",sep=",",header=TRUE)
betweenness.47 <- read.table("batch-data/network-influence-power/out.47/betweenness.csv",sep=",",header=TRUE)
betweenness.67 <- read.table("batch-data/network-influence-power/out.67/betweenness.csv",sep=",",header=TRUE)
betweenness.71 <- read.table("batch-data/network-influence-power/out.71/betweenness.csv",sep=",",header=TRUE)
betweenness.72 <- read.table("batch-data/network-influence-power/out.72/betweenness.csv",sep=",",header=TRUE)
betweenness.79 <- read.table("batch-data/network-influence-power/out.79/betweenness.csv",sep=",",header=TRUE)
network.15 <- cbind(nodes.15,betweenness.15)
network.26 <- cbind(nodes.26,betweenness.26)
network.27 <- cbind(nodes.27,betweenness.27)
network.30 <- cbind(nodes.30,betweenness.30)
network.36 <- cbind(nodes.36,betweenness.36)
network.45 <- cbind(nodes.45,betweenness.45)
network.46 <- cbind(nodes.46,betweenness.46)
network.47 <- cbind(nodes.47,betweenness.47)
network.67 <- cbind(nodes.67,betweenness.67)
network.71 <- cbind(nodes.71,betweenness.71)
network.72 <- cbind(nodes.72,betweenness.72)
network.79 <- cbind(nodes.79,betweenness.79)

allnetworks <- rbind(network.15, network.26, network.27, network.30, network.36,
	network.45, network.46, network.47, network.67, network.71, 
	network.72, network.79)

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

doplot <- function(xv, yv, filename, xlab, xlabelpos, ylabelpos)
{
  plot(xv,yv,xlab=xlab,ylab="influence")
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
  dev.copy2eps(file=filename)
  dev.off()
}

doplot(network.26$in.degree,network.26$dp.dr.i.,
       "out-26.influence.v.indegree.eps","in-degree",0.28,1.8)

doplot(network.26$out.degree,network.26$dp.dr.i.,
	"out-26.influence.v.outdegree.eps","out-degree",2.5,0.2)

doplot(network.26$temperature,network.26$dp.dr.i.,
	"out-26.influence.v.temperature.eps","temperature",-0.01,1.8)

doplot(network.26$in.closeness,network.26$dp.dr.i.,
	"out-26.influence.v.incloseness.eps","in-closeness",0.27,1.8)

doplot(network.26$out.closeness,network.26$dp.dr.i.,
	"out-26.influence.v.outcloseness.eps","out-closeness",0.35,0.2)

doplot(network.26$betweenness,network.26$dp.dr.i.,
	"out-26.influence.v.betweenness.eps","betweenness",-5,1.8)

doplot(network.26$undirected.eigenvector.centrality,network.26$dp.dr.i.,
       "out-26.influence.v.eigenvector.centrality.eps",
       "eigenvector centrality (undirected)",0.045,1.8)
	
doplot(network.26$pagerank,network.26$dp.dr.i.,
	"out-26.influence.v.pagerank.eps","pagerank",0.002,1.8)


doplot(allnetworks$in.degree,allnetworks$dp.dr.i.,
	"all.influence.v.indegree.eps","in-degree",0.28,1.8)

doplot(allnetworks$out.degree,allnetworks$dp.dr.i.,
	"all.influence.v.outdegree.eps","out-degree",2.1,1.8)

doplot(allnetworks$temperature,allnetworks$dp.dr.i.,
	"all.influence.v.temperature.eps","temperature",-0.01,1.8)

doplot(allnetworks$in.closeness,allnetworks$dp.dr.i.,
	"all.influence.v.incloseness.eps","in-closeness",0.27,1.8)

doplot(allnetworks$out.closeness,allnetworks$dp.dr.i.,
	"all.influence.v.outcloseness.eps","out-closeness",0.62,1.8)

doplot(allnetworks$betweenness,allnetworks$dp.dr.i.,
	"all.influence.v.betweenness.eps","betweenness",-5,1.8)

doplot(allnetworks$undirected.eigenvector.centrality,allnetworks$dp.dr.i.,
       "all.influence.v.eigenvector.centrality.eps",
       "eigenvector centrality (undirected)",0.02,1.8)

doplot(allnetworks$pagerank,allnetworks$dp.dr.i.,
       "all.influence.v.pagerank.eps","pagerank",0.002,1.8)

doplot(allnetworks$constraint,allnetworks$dp.dr.i.,
       "all.influence.v.constraint.eps","Burt's constraint",0.2,1.8)

