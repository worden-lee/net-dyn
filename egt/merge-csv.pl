#!/usr/bin/perl
use warnings;
use strict;
use IO::File;

my %options =
( path => "batch-power",
  allnet => "",
  allnodes => ""
);

while (@ARGV)
{ my $arg = shift @ARGV;
  # else
  $options{path} = $arg;
}
if ($options{allnet} eq "")
{ $options{allnet} = "$options{path}/allnet.csv";
  $options{allnodes} = "$options{path}/allnodes.csv";
}

my $allnet = new IO::File;
my $allnodes = new IO::File;
$::wrote_headers = 0;
open($allnet,">$options{allnet}")
  or die "couldn't open $options{allnet} for writing: $!";
open($allnodes,">$options{allnodes}")
  or die "couldn't open $options{allnodes} for writing: $!";
recurse_dir($options{path}, $allnet, $allnodes);
close $allnet;
close $allnodes;
exit;

# ===== end of main code, beginning of subroutines =====

# why reinvent the wheel? I don't know.
# this reads a csv file, returns a hash object with 3 things in it:
#  "names" => list of column names.
#  "indices" => hash mapping column names to column numbers.
#  "rows" => list of rows.  each entry is a list of values.
sub parse_csv
{ my($file) = @_;
  my(%data);
  print "$file\n";
  open(CSV,"$file")
    or die "couldn't open $file for reading: $!";
  while(<CSV>)
  { chomp;
    if ($. == 1) # header line
    { s/^\#\s*//;
      my @names = split(/,/);
      @{$data{names}} = @names;
      @{$data{indices}}{@names} = ($[ .. $#names);
    }
    else
    { my @vals = split(/,/);
      push @{$data{rows}}, \@vals;
    }
  }
  close CSV;
  return \%data;
}

# this reads a .settings file and returns the same structure as parse_csv
sub parse_settings
{ my($file) = @_;
  my(%data);
#   print "$file\n";
  open(SETTINGS,"$file")
    or die "couldn't open $file for reading: $!";
  $data{rows} = [ [] ];
  while(<SETTINGS>)
  { chomp;
    next if (/^#/);
    # hopefully no need to worry about 'include' lines
    if (/\S\s+\S/)
    { my($key,$value) = split(/\s+/);
      push @{$data{names}}, $key;
      $data{indices}{$key} = scalar(@{$data{rows}->[0]});
      push @{$data{rows}[0]}, $value;
    }
  }
  close SETTINGS;

#   print "# ",join(",",@{$data->{names}}),"\n";
#   for my $l (@{$data->{rows}})
#   { print join(",",@{$l}),"\n"; }

  return \%data;
}

sub get_column
{ my($data,$name) = @_;
  my $index = $data->{indices}{$name};
  return map { $_->[$index] } @{$data->{rows}};
}

sub get_row_column
{ my($data,$row,$col) = @_;
  my $col_index = $data->{indices}{$col};
  my $val = $data->{rows}[$row][$col_index];
#   print "get_row_column $row $col($col_index) - $val\n";
  return $val;
}

# pass a 'data' structure and a column, as refs, and a column name
sub add_column
{ my($data,$name,$column) = @_;
  push @{$data->{names}}, $name;
  for my $i ($[ .. (scalar(@{$data->{rows}}) + $[ - 1))
  { push @{$data->{rows}[$i]}, $column->[$i]; }
  $data->{indices}{$name} = (scalar(@{$data->{rows}[0]}) + $[ - 1);
#   print "is ",$column->[0]," equal to ",
#     $data->{rows}[0][$data->{indices}{$name}],"?\n";
}

sub process_files
{ my($dir, $allnet, $allnodes) = @_;
  my($net) = parse_csv("$dir/network.csv");
  my($nodes) = parse_csv("$dir/nodes.csv");
  my($settings) = parse_settings("$dir/all.settings");
  # disambiguate column name
  map {s/^p$/p.graph/} @{$net->{names}};
  $net->{indices}{"p.graph"} = $net->{indices}{"p"};
  map {s/^p$/p.vertex/} @{$nodes->{names}};
  $nodes->{indices}{"p.vertex"} = $nodes->{indices}{"p"};
  # create rep column
  my $rep;
  if ("$dir/" =~ m|out\.([^/]*)/|) { $rep = $1; }
  add_column($net,"rep",[$rep]);
  # create s column
#   print "mutantFitness index = ",
#     $settings->{indices}{"mutantFitness"}, "\n";
  my $mf = get_row_column($settings,0,"mutantFitness");
#   print "mutantFitness = ",$mf,"\n";
  my $rf = get_row_column($settings,0,"residentFitness");
  my $s = ($mf-$rf)/$rf;
  add_column($net,"s",[$s]);
  # now a nice calculation
  my $rule = get_row_column($settings,0,"update_rule");
  if ($rule eq "EGT")
  { my $mu_m10 = get_row_column($net,0,"mu.-1.0");
    my $N = get_row_column($net,0,"nv");
    my $denom = 1/(1 - exp(-$s*$N));
    my @pcol;
    my $sum_p = 0;
    my $xi_index = $nodes->{indices}{"in degree"};
    for my $vx (@{$nodes->{rows}})
    { my $xi = $vx->[$xi_index];
      my $pfix = (1 - exp(-$s/($mu_m10*$xi)))*$denom;
      push @pcol, $pfix;
      $sum_p = $sum_p + $pfix;
    }
    add_column($nodes, "analytic.p.vertex", \@pcol);
    add_column($net, "analytic.p.graph", [ $sum_p/$N ]);
  }
  elsif ($rule eq "VM")
  { my $mu_01 = get_row_column($net,0,"mu.0.1");
    my $s_adj = $s*$mu_01/get_row_column($net,0,"mu.0.2");
    my $N = get_row_column($net,0,"nv");
    my $denom = 1/(1 - exp(-$s_adj*$N*$mu_01));
    my @pcol;
    my $sum_p = 0;
    my $xo_index = $nodes->{indices}{"out degree"};
    for my $vx (@{$nodes->{rows}})
    { my $xo = $vx->[$xo_index];
      my $pfix = (1 - exp(-$s_adj*$xo))*$denom;
      push @pcol, $pfix;
      $sum_p = $sum_p + $pfix;
    }
    add_column($nodes, "analytic.p.vertex", \@pcol);
    add_column($net, "analytic.p.graph", [ $sum_p/$N ]);
  }

  # test the data handling
#   print "# ",join(",",@{$net->{names}}),"\n";
#   for my $l (@{$net->{rows}})
#   { print join(",",@{$l}),"\n"; }
#   # print a column
#   print "# p columns:\n";
#   for my $l (@{$net->{rows}})
#   { print @{$l}[$net->{indices}{"p.graph"}]," ",
#       @{$l}[$net->{indices}{"analytic.p.graph"}],"\n";
#   }
#   for my $l (@{$nodes->{rows}})
#   { print @{$l}[$nodes->{indices}{"p.vertex"}]," ",
#       @{$l}[$nodes->{indices}{"analytic.p.vertex"}],"\n";
#   }
  # write em out
  my $netrow; # we're going to join the unique row of "net" to
              # all the rows of "nodes"
  my $local_allnet = "$dir/network.extended.csv";
  open(LOCAL_ALLNET,">$local_allnet")
    or die "couldn't open $local_allnet for writing: $!";
  print LOCAL_ALLNET "# ",join(",",@{$net->{names}}),"\n";
  for my $l (@{$net->{rows}})
  { print LOCAL_ALLNET join(",",@{$l}),"\n";
    $netrow = $l;
  }
  close LOCAL_ALLNET;
  my $local_allnodes = "$dir/nodes.extended.csv";
  open(LOCAL_ALLNODES,">$local_allnodes")
    or die "couldn't open $local_allnodes for writing: $!";
  print LOCAL_ALLNODES
    "# ",join(",",@{$nodes->{names}},@{$net->{names}}),"\n";
  for my $l (@{$nodes->{rows}})
  { print LOCAL_ALLNODES
      join(",",@{$l},@{$netrow}),"\n";
  }
  close LOCAL_ALLNODES;
  if (!$::wrote_headers)
  { print $allnet "# ",join(",",@{$net->{names}}),"\n";
    print $allnodes
      "# ",join(",",@{$nodes->{names}},@{$net->{names}}),"\n";
    $::wrote_headers = 1;
  }
  for my $l (@{$net->{rows}})
  { print $allnet join(",",@{$l}),"\n"; }
  for my $l (@{$nodes->{rows}})
  { print $allnodes
      join(",",@{$l},@{$netrow}),"\n"; }
}

sub recurse_dir
{ my($dir, $allnet, $allnodes) = @_;
  print "examining $dir\n";
  if (-f "$dir/nodes.csv")
  { process_files($dir, $allnet, $allnodes);
  }
  else
  { opendir(DIR,$dir) or die "couldn't open $dir for recursion: $!";
    my(@subdirs) =
      grep { $_ ne "." and $_ ne ".." and -d "$dir/$_" } readdir(DIR);
    foreach (@subdirs)
    { recurse_dir("$dir/$_",$allnet,$allnodes); }
    closedir(DIR);
  }
}

