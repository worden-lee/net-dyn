#!/usr/bin/perl
while(<>)
{
  my($hx);
  if (/color=\"#(..)0000\"/)
  { $hx = $1;
    my($number) = hex($hx) / 256;
    $number *= $number;
#    my($nhx) = sprintf("%02x",256*$number);
    my($whx) = sprintf("%02x",256*(1 - $number));
    s/${hx}0000/ff$whx$whx/;
  }
  print;
}
