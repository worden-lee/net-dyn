# source("load-nodes.r")

nodes15 <- read.table("batch-data/network-influence-power/out.15/nodes.csv",sep=",",header=TRUE)
nodes26 <- read.table("batch-data/network-influence-power/out.26/nodes.csv",sep=",",header=TRUE)
nodes27 <- read.table("batch-data/network-influence-power/out.27/nodes.csv",sep=",",header=TRUE)
nodes30 <- read.table("batch-data/network-influence-power/out.30/nodes.csv",sep=",",header=TRUE)
nodes36 <- read.table("batch-data/network-influence-power/out.36/nodes.csv",sep=",",header=TRUE)
nodes45 <- read.table("batch-data/network-influence-power/out.45/nodes.csv",sep=",",header=TRUE)
nodes46 <- read.table("batch-data/network-influence-power/out.46/nodes.csv",sep=",",header=TRUE)
nodes47 <- read.table("batch-data/network-influence-power/out.47/nodes.csv",sep=",",header=TRUE)
nodes67 <- read.table("batch-data/network-influence-power/out.67/nodes.csv",sep=",",header=TRUE)
nodes71 <- read.table("batch-data/network-influence-power/out.71/nodes.csv",sep=",",header=TRUE)
nodes72 <- read.table("batch-data/network-influence-power/out.72/nodes.csv",sep=",",header=TRUE)
nodes79 <- read.table("batch-data/network-influence-power/out.79/nodes.csv",sep=",",header=TRUE)

betweenness15 <- read.table("batch-data/network-influence-power/out.15/betweenness.csv",sep=",",header=TRUE)
betweenness26 <- read.table("batch-data/network-influence-power/out.26/betweenness.csv",sep=",",header=TRUE)
betweenness27 <- read.table("batch-data/network-influence-power/out.27/betweenness.csv",sep=",",header=TRUE)
betweenness30 <- read.table("batch-data/network-influence-power/out.30/betweenness.csv",sep=",",header=TRUE)
betweenness36 <- read.table("batch-data/network-influence-power/out.36/betweenness.csv",sep=",",header=TRUE)
betweenness45 <- read.table("batch-data/network-influence-power/out.45/betweenness.csv",sep=",",header=TRUE)
betweenness46 <- read.table("batch-data/network-influence-power/out.46/betweenness.csv",sep=",",header=TRUE)
betweenness47 <- read.table("batch-data/network-influence-power/out.47/betweenness.csv",sep=",",header=TRUE)
betweenness67 <- read.table("batch-data/network-influence-power/out.67/betweenness.csv",sep=",",header=TRUE)
betweenness71 <- read.table("batch-data/network-influence-power/out.71/betweenness.csv",sep=",",header=TRUE)
betweenness72 <- read.table("batch-data/network-influence-power/out.72/betweenness.csv",sep=",",header=TRUE)
betweenness79 <- read.table("batch-data/network-influence-power/out.79/betweenness.csv",sep=",",header=TRUE)
network15 <- cbind(nodes15,betweenness15)
network26 <- cbind(nodes26,betweenness26)
network27 <- cbind(nodes27,betweenness27)
network30 <- cbind(nodes30,betweenness30)
network36 <- cbind(nodes36,betweenness36)
network45 <- cbind(nodes45,betweenness45)
network46 <- cbind(nodes46,betweenness46)
network47 <- cbind(nodes47,betweenness47)
network67 <- cbind(nodes67,betweenness67)
network71 <- cbind(nodes71,betweenness71)
network72 <- cbind(nodes72,betweenness72)
network79 <- cbind(nodes79,betweenness79)

allnetworks <- rbind(network15, network26, network27, network30, network36,
	network45, network46, network47, network67, network71, 
	network72, network79)

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

doplot <- function(xv, yv, filename, xlabelpos, ylabelpos)
{
  plot(xv,yv)
  lm <- lm(yv ~ xv)
  lms <- summary(lm)
  abline(lm,col="blue")
  text(xlabelpos,ylabelpos,
	sprintf("r-squared=%g\nslope=%g\nstd.err(slope)=%g\np=%g %s",
	lms$r.squared, lms$coefficients[2,1], lms$coefficients[2,2],
	lms$coefficients[2,4], stars(lms$coefficients[2,4])),
	pos=4)
#  dev.copy2eps(file=filename)
#  dev.off()
}

doplot(network26$in.degree,network26$dp.dr.i.,
	"out-26.influence.v.indegree.eps",0.28,0.033)

doplot(network26$out.degree,network26$dp.dr.i.,
	"out-26.influence.v.outdegree.eps",2.5,0.0045)

doplot(network26$temperature,network26$dp.dr.i.,
	"out-26.influence.v.temperature.eps",0,0.033)

doplot(network26$in.closeness,network26$dp.dr.i.,
	"out-26.influence.v.incloseness.eps",0.28,0.033)

doplot(network26$out.closeness,network26$dp.dr.i.,
	"out-26.influence.v.outcloseness.eps",0.35,0.0045)

doplot(network26$betweenness,network26$dp.dr.i.,
	"out-26.influence.v.betweenness.eps",0.28,0.033)

doplot(network26$undirected.eigenvector.centrality,network26$dp.dr.i.,
	"out-26.influence.v.betweenness.eps",0.045,0.033)
	
doplot(network26$pagerank,network26$dp.dr.i.,
	"out-26.influence.v.pagerank.eps",0.002,0.033)


doplot(allnetworks$in.degree,allnetworks$dp.dr.i.,
	"all.influence.v.indegree.eps",0.28,0.033)

doplot(allnetworks$out.degree,allnetworks$dp.dr.i.,
	"all.influence.v.outdegree.eps",2.1,0.033)

doplot(allnetworks$temperature,allnetworks$dp.dr.i.,
	"all.influence.v.temperature.eps",0.28,0.033)

doplot(allnetworks$in.closeness,allnetworks$dp.dr.i.,
	"all.influence.v.incloseness.eps",0.27,0.033)

doplot(allnetworks$out.closeness,allnetworks$dp.dr.i.,
	"all.influence.v.outcloseness.eps",0.35,0.033)

doplot(allnetworks$betweenness,allnetworks$dp.dr.i.,
	"all.influence.v.betweenness.eps",0,0.033)

doplot(allnetworks$undirected.eigenvector.centrality,allnetworks$dp.dr.i.,
	"all.influence.v.betweenness.eps",0.02,0.033)

doplot(allnetworks$pagerank,allnetworks$dp.dr.i.,
	"all.influence.v.pagerank.eps",0.002,0.033)

