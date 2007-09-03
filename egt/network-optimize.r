# trying some stuff out, learning to write R code -lw

# code for adjusting a network to optimize a given indicator

# this is not super documented but here's how you might run it:
#     start R
# $ R
#     load in this file
# > source("network-optimize.r")
#     create an empty network with 5 nodes
# > n <- network.initialize(5)
#     fix it up until its density is about 0.3 (i.e. close to 3/10 of
#     the possible edges are present)
#     this will display a twitchy animation of the network changing
# > n <- improve.network(n,density.indicator(0.3))
#     now start looking for ways to improve the network's propensity
#     to amplify selection, i.e. improve fixation probability of a fitter
#     node occupant.  This will animate the process of fixation of the
#     mutant (plotted in green, I think) as well as changing network.
#     it will run for a long time and you'll probably give up waiting.
# > n <- improve.network(n,selection.indicator)

# we assume a network is directed, allow networks with
# loops, but do not allow parallel edges

library("network")

# =================================================================
# default parameter values, you can change them
network.evolution.peak.verifications = 5;
network.evolution.fixation.fitness   = 1.1;
network.evolution.fixation.nsamples  = 300;
network.evolution.fixation.maxsteps  = 10000;
network.evolution.fixation.animate   = TRUE;

# note, for 5 nodes the Moran process fixation probability
# for fitness 1.1 is 0.2398

# note also, N = var / (accuracy * mean)^2
# for p = 0.25, accuracy = 0.1, we find N = 300
# (Bernoulli distribution, mean = p, var = p(1-p))

# =================================================================
# here, code for optimizing networks for a given indicator
#  you need an indicator function to indicate what to optimize,
#  these are found further below

# modify a network by turning a given edge on or off
# returns: modified network
flipedge <- function(net,pair)
{
  i <- pair[1];
  j <- pair[2];
  net[i,j] <- 1 - net[i,j];
  return(net);
}

# what's the indicator if a given edge is flipped
# returns: value of indicator function
trymove <- function(net,indicator,pair)
{ indicator(flipedge(net,pair));
}

# what edge makes the best improvement to the indicator when flipped
# returns: a pair, or an empty list when there's no good move
bestmove.global <- function(net,indicator)
{
  best <- indicator(net);
cat("given network: ",best,"\n");
  bestmove <- c();
  vs <- network.size(net);
  ilist <- sample(1:vs); # random permutation
  jlist <- sample(1:vs);
  for (i in ilist)
    for (j in jlist)
    { # try flipping edge ij on/off
      last <- trymove(net,indicator,c(i,j));
cat("toggle ",c(i,j),": ",last,"\n");
      if (last > best)
      { best <- last;
        bestmove <- c(i,j);
      }
    }
cat("best = ",best,": toggle ",c(i,j),"\n");
  return(bestmove);
}

# look for an edge that improves the indicator when flipped, return the
# first one that does
bestmove.greedy <- function(net,indicator)
{
  given <- indicator(net);
cat("given network: ",given,"\n");
  vs <- network.size(net);
  ilist <- sample(1:vs); # random permutation
  jlist <- sample(1:vs);
  for (i in ilist)
    for (j in jlist)
      if (j != i)
      { # try flipping edge ij on/off
        last <- trymove(net,indicator,c(i,j));
cat("toggle ",c(i,j),": ",last,"\n");
        if (last > given)
        { return(c(i,j)); }
      }
  return(c());
}
  
# modify a network by flipping whatever edge improves the indicator the most
# returns: modified network
improve.network.once <- function(net,indicator)
{ edge <- bestmove.greedy(net,indicator);
  if (length(edge) == 2)
  { return(flipedge(net,edge)); }
  else
  { cat("no improvement found\n");
    return(net);
  }
}

# repeatedly improve network until a local optimum (relative to flipping
# a single edge) is reached.  plot each step on the screen.
# returns: locally optimal network
improve.network <-
  function(net, indicator,
           n.verifications=network.evolution.peak.verifications)
{ verified <- 0;
  plot(net,displaylabels=TRUE);
  while(1)
  { net1 <- improve.network.once(net,indicator);
    plot(net1,displaylabels=TRUE);
    if (identical(net,net1))
    { cat("networks are the same:\n");
      print(as.matrix.network.adjacency(net));
      print(as.matrix.network.adjacency(net1));
      verified <- verified + 1;
      if (verified == n.verifications)
      { break; }
    }
    else
    { verified <- 0; }
    net <- net1;
# keep a record in case of interruption
    .last.network <<- net;
    #Sys.sleep(1);
  }
  return(net);
}

# =================================================================
# indicator functions

# factory for network density indicator functions
# returns: a function
# density(d) is a function
# density(d)(net) is how close the density of net is to d
density.indicator <- function(d)
{ return(function(net) { 1 - abs(network.density(net) - d); });
}

# the serious indicator function is the fixation probability.
# it requires a couple of subroutines.

# this picks from a list of weights
# returns: index into the list
# i.e. given list [3,5,1,2], it will return 1 with probability 3/11,
# 2 with probability 5/11, etc. (the sum of the weights being 11).
pick.weighted <- function(weights)
{
#cat("pick.weighted from ",weights,": ");
  cum.weight <- runif(1,0,sum(weights));
  for (index in 1:length(weights))
  { cum.weight <- cum.weight - weights[index];
    if (cum.weight < 0)
    { #cat(index,"\n");
      return(index);
    }
  }
}
  
# this does one reproduction on a network
# state is a list: fitness at each network vertex
# returns: new state
reproduction.step <- function(net,state)
{
  if(network.edgecount(net)==0) { return(state); }
  while(1)
  { # first choose a parent node randomly, proportionately to fitness
    parent <- pick.weighted(state);
    # now choose child down an arc from the parent
    out.edges <- get.edges(net,parent,neighborhood="out");
    if (length(out.edges)>0) {break;}
  }
  # #@($(*!! arrows point from $outl to $inl!!
  # note also, n[1,2] points from 1 to 2. i.e. $outl is 1
  child.candidates <- unlist(lapply(out.edges,function(ed){ed$inl}));
#cat("choose child from ",child.candidates,"\n");
  out.weights <-
    unlist(lapply(child.candidates,function(ch){net[parent,ch]}));
#cat("weights: ",out.weights,"\n");
  child.index <- pick.weighted(out.weights);
  child <- child.candidates[child.index];
#cat("n[",parent,",",child,"] = ",n[parent,child],"\n");
  # replace child with parent type, i.e. fitness
  state[child] <- state[parent];
  return(state);
}

# this drops a mutant in and sees if it fixates.
# with animation
# returns: which fitness fixated.
#   the mean fitness, if no fixation.
montecarlofixation.once <-
  function(net, mutant.fitness,
           animate=network.evolution.fixation.animate,
           max.steps=network.evolution.fixation.maxsteps)
{ if(network.edgecount(net)==0) { return(1); }
  state <- rep(1,network.size(net)); # start with all 1's
  # pick a node and mutate it
  mutant.index <- pick.weighted(state);
  state[mutant.index] <- mutant.fitness;
  if(animate && !identical(.cache.network,net))
  { # cache network layout for plotting
    .cache.network <<- net;
    .cache.network.coords <<-
      network.layout.fruchtermanreingold(as.matrix.network.adjacency(net),
                                         NULL);
  }
  step <- 1;
  while(length(unique(state)) > 1 && step <= max.steps)
  { state <- reproduction.step(net,state);
    if (animate)
    { #Sys.sleep(0.5);
      plot(net, vertex.col=ifelse(state<=1,2,3),
           coord=.cache.network.coords, jitter=FALSE);
    }
    step <- step + 1;
  }
#  cat("after ",step," steps: ",state,": mean = ",mean(state),"\n");
  return(mean(state));
}
.cache.network <<- 0;

# this tests for fixation a bunch of times and estimates
# fixation probability.
# returns: fixation probability
#  adjusted to include 'partial fixation'
montecarlofixation.batch <-
  function(net, repetitions=network.evolution.fixation.nsamples,
           mutant.fitness=network.evolution.fixation.fitness)
{ total.fitness <- 0;
  fitnesses <-
    replicate(repetitions,montecarlofixation.once(net,mutant.fitness));
##   for (rep in 1:repetitions)
##   { total.fitness <-
##       total.fitness + montecarlofixation.once(net,mutant.fitness)
##   }
  mean.fitness = mean(fitnesses);
#cat(fitnesses,": mean = ",mean.fitness,"\n");
  # keep cached results over multiple batches
  
  return( (mean.fitness - 1)/(mutant.fitness - 1) );
}

# indicator function reporting fixation probability of a favorable
#  mutant
# note, this indicator won't get you anywhere if your graph isn't
#  at least connected and not multirooted. might require more than that.
selection.indicator <- function(net)
{ return(montecarlofixation.batch(net)); }

# for testing
cacher <- function()
{ new.samples <- runif(100,0,1);
  print(new.samples);
  .cached.mean <- .cached.nsamples * .cached.mean + 100 * mean(new.samples);
  .cached.nsamples <- .cached.nsamples + 100;
  .cached.mean <- .cached.mean / .cached.nsamples;
  cat("long-term average =",.cached.mean,"\n");
}
