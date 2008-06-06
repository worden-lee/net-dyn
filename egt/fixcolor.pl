#!/usr/bin/perl
while(<>)
{
  my($hx);
  if (/color=\"#(..)0000\"/)
  { $hx = $1;
    my($number) = hex($hx) / 256;
    $number *= $number;
#    my($nhx) = sprintf("%02x",256*$number);
    my($wh) = 256*(1 - $number);
    if ($wh > 255) { $wh = 255; }
    my($whx) = sprintf("%02x",$wh);
    s/${hx}0000/ff$whx$whx/;
    s/\bcolor=/style=filled, fillcolor=/;
  }
  print;
}
