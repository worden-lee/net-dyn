#!/usr/bin/perl
use File::Path;

my $reps = 10000;
my @exprange = (1.0,3.0);
my $expstep = 0.1;

my @explist =
 map { $_ * $expstep } (($exprange[0]/$expstep) .. ($exprange[1]/$expstep));

my $pwd = `pwd`;
chomp($pwd);

my $experiment = "network-power-sample";
my $outdir = "$pwd/batch-data/$experiment";

if (!-e $outdir)
{ mkpath($outdir) or die "couldn't create $outdir";
#   symlink("$outdir/network","$outdir/network.dot")
#     or die "couldn't create 'network' symbolic link";
}

print "@explist\n";

for my $i (1 .. $reps)
{ for my $inexp (@explist)
  { for my $outexp (@explist)
    { my @extra_args = ("pl_in_exp=$inexp","pl_out_exp=$outexp");
      my($catdir) = "$outdir/".join(".",@extra_args);
      my($dest) = "$catdir/out.$i";
      next if (-e $dest);
      if (!-e "out") { mkdir("out") or die "couldn't mkdir out"; }
      system("rm -rf out/*");
      my $comm = "./network-sample -f settings/$experiment.settings ".
        join("", map {" --$_" } @extra_args);
      print "$comm\n";
      system($comm) == 0 or die "error running $comm";
      if (!-e $catdir) { mkpath($catdir) or die "couldn't create $catdir"; }
      $comm = "cp -r --force --backup=numbered out/ $dest";
      print "$comm\n";
      system($comm) == 0 or die "error running $comm";
    }
  }
}
