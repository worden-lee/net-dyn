#!/usr/bin/perl
use File::Path;

my $reps = 10000;

my $pwd = `pwd`;
chomp($pwd);
# parameters are compiled into the program, so this is awkward for now.
# set this to match what's compiled in, we can improve it later.
my $experiment = "network-power-sample";
my $outdir = "$pwd/batch-data/$experiment";

if (!-e $outdir)
{ mkpath($outdir) or die "couldn't create $outdir";
#   symlink("$outdir/network","$outdir/network.dot")
#     or die "couldn't create 'network' symbolic link";
}

for my $i (1 .. $reps)
{ my($dest) = "$outdir/out.$i";
  next if (-e $dest);
  if (!-e "out") { mkdir("out") or die "couldn't mkdir out"; }
  system("rm -rf out/*");
  my $comm = "./network-optimize -f settings/$experiment.settings";
  system($comm) == 0 or die "error running $comm";
  $comm = "cp -r --force --backup=numbered out/ $dest";
  system($comm) == 0 or die "error running $comm";
}
