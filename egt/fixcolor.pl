#!/usr/bin/perl

# comm. line options
my($dotfile,$csvfile,$outfile) = ("","","");

while (@ARGV)
{
  my $arg = shift ARGV;
  print "arg: $arg\n";
  if ($arg eq "-d")
  { $dotfile = shift ARGV;
    $dbg && print "dotfile $dotfile\n";
  }
  elsif ($arg eq "-c")
  { $csvfile = shift ARGV;
    $dbg && print "csvfile $csvfile\n";
  }
  elsif ($arg eq "-o")
  { $outfile = shift ARGV;
    $dbg && print "output to $outfile\n";
  }
  elsif ($arg eq "--debug")
  { $dbg = 1;
    print "debug output on\n";
  }
  else
  { if ($dotfile eq "")
    { $dotfile = $arg;
      $dbg && print "dotfile $dotfile\n";
    }
    elsif ($csvfile eq "")
    { $csvfile = $arg;
      $dbg && print "csvfile $csvfile\n";
    }
  }
}

open DOT, "<$dotfile"
  or die "couldn't open $dotfile";

open CSV, "<$csvfile"
  or die "couldn't open $csvfile";

open OUT, ">$outfile"
  or die "couldn't open $outfile for writing";

my $csline = <CSV>;
$dbg && print "read first line: $csline";
my(@names) = split(/,/,$csline);
my($iix);
for $i ($[ .. $#names)
{ if ($names[$i] =~ /^dp.dr.i.$/i or $names[$i] =~ /^influence$/i)
  { $iix = $i;
    last;
  }
}
my(%influence);
my($maxinf) = 0;
my($counter) = 0;
while(<CSV>)
{ my(@csline) = split(/,/);
  $influence{$counter} = $csline[$iix];
  $dbg && print "influence{$counter} = $csline[$iix]\n";
  if ($maxinf < $csline[$iix])
  { $maxinf = $csline[$iix]; }
  ++$counter;
}

# normalize the influences
for my $k (keys(%influence))
{ $influence{$k} /= $maxinf;
  $dbg && print "normalize influence{$k} = $influence{$k}\n";
}

while(<DOT>)
{ my $line = $_;
  chomp $line;
  if ($line =~ /^\s*"(\S+)"(\s*\[(.*,)*)(color=[^,\]]*)(.*)$/
      || $line =~ /^\s*(\S+)(\s*\[(.*,)*)(color=[^,\]]*)(.*)/)
  { $dbg && print "caught line [$1 $2 $3 $4 $5] <$line>: ";
    my($key, $color) = ($1, $4);
    $dbg && print "[key=$key] ";
    if (defined $influence{$key})
    { my($number) = 255*(1-$influence{$key});
      my($whx) = sprintf("%02x",$number);
      $dbg && print "fixing [$key] ";
      $line = "$1$2 style=filled,fillcolor=\"#ff${whx}${whx}\"$5";
    }
    $dbg && print "<$line>\n";
  }
  print OUT "$line\n";
}
