#!/usr/bin/perl
use File::Path;

my $reps = 100;

# parameters are compiled into the program, so this is awkward for now.
# set this to match what's compiled in, we can improve it later.
my $outdir = "batch-data/final-networks/12/choose-parent-drift";

if (!-e $outdir)
{ mkpath($outdir) or die "couldn't create $outdir";
  symlink("$outdir/network","$outdir/network.dot")
    or die "couldn't create 'network' symbolic link";
}

for my $i (1 .. $reps)
{ my $comm = "./network-optimize";
  system($comm) == 0 or die "error running $comm";
  $comm = "cp --force --backup=numbered --suffix=.dot dot/network.dot $outdir/network";
  system($comm) == 0 or die "error running $comm";
}
