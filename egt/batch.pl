#!/usr/bin/perl
use File::Path;

my $reps = 100;

my $pwd = `pwd`;
chomp($pwd);
# parameters are compiled into the program, so this is awkward for now.
# set this to match what's compiled in, we can improve it later.
my $outdir = "$pwd/batch-data/final-networks/12/choose-parent";

if (!-e $outdir)
{ mkpath($outdir) or die "couldn't create $outdir";
#   symlink("$outdir/network","$outdir/network.dot")
#     or die "couldn't create 'network' symbolic link";
}

for my $i (1 .. $reps)
{ if (!-e "out") { mkdir("out") or die "couldn't mkdir out"; }
  my $comm = "./network-optimize";
  system($comm) == 0 or die "error running $comm";
  $comm = "cp -r --force --backup=numbered out/ $outdir/out";
  system($comm) == 0 or die "error running $comm";
}
