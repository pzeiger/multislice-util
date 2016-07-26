#!/usr/bin/perl

open inp, "<$ARGV[0]";
open out, ">>$ARGV[1]";
#open out, ">$ARGV[1]";

$N="$ARGV[2]";

print "$N\n";

#for($i=1;$i<8;$i++) {
#  $l=<inp1>;
#}

#print out "\n";
#print out "  F\n";

for($i=1;$i<=$N;$i++) {

  $l = <inp>;
  @tar = split /\s+/, $l;

#  print "$tar[1] "; print "$tar[2] "; print "$tar[3] "; print "$tar[4]\n";

  if($tar[2]<0) {$tar[2]+=1;}
  if($tar[3]<0) {$tar[3]+=1;}
  if($tar[4]<0) {$tar[4]+=1;}

  if($tar[2]>1) {$tar[2]-=1;}
  if($tar[3]>1) {$tar[3]-=1;}
  if($tar[4]>1) {$tar[4]-=1;}


  printf out "%2i %1.10f %1.10f %1.10f\n", $tar[1], $tar[2], $tar[3], $tar[4];
}
