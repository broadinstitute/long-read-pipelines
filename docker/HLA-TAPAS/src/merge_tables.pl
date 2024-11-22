#!/usr/bin/perl
#

use strict;

my $Table1 = $ARGV[0];
my $Table2 = $ARGV[1];
my $IndexStr = $ARGV[2];

if (scalar(@ARGV) != 3) {
    print "usage: %>merge_tables.pl datafile_1 datafile_2 index_string\n";
    print "prints all rows of datafile_2, followed by corresponding rows (if available) from datafile_1\n";
    exit();
}   


my @headers1 = ();
my %data1 = ();
open(T1,$Table1);

my $linecount = 0;
my $IndexCol = -1;

while(my $c = <T1>) {
  $c=~s/\s+$//;
  $c=~s/^\s+//;
  my @line = split /\s+/, $c;

  if ( $linecount == 0 ) {
    for (my $i=0; $i<=$#line; $i++) {
      if ( $line[$i] eq $IndexStr ) {
        if ( $IndexCol >= 0 ) {
          die "Duplicate column label $line[$i] in $Table1 - exiting.\n";
        }
        $IndexCol = $i;
      } 
      $headers1[$i] = $line[$i];
    }
  }
  else {
    if ( $IndexCol == -1 ) {
      die "Did not find label $IndexStr in $Table1 - exiting.\n";
    }

    for (my $i=0; $i<=$#line; $i++) {
      if ( $line[$IndexCol] ne "NA" ) {
        $data1{$line[$IndexCol]}{$headers1[$i]} = $line[$i];
      }
    }
  } 

  $linecount++;
}
close(T1);

#print STDERR "read $linecount lines from $Table1\n";

open(T2,$Table2);

$linecount = 0;
$IndexCol = -1;

while(my $c = <T2>) {
  $c=~s/\s+$//;
  $c=~s/^\s+//;
  my @line = split /\s+/, $c;

  if ( $linecount == 0 ) {
    for (my $i=0; $i<=$#line; $i++) {
      if ( $line[$i] eq $IndexStr ) {
        if ( $IndexCol >= 0 ) {
          die "Duplicate column label $line[$i] in $Table2 - exiting.\n";
        }
        $IndexCol = $i;
      }
      print "$line[$i]\t";
    }
    foreach my $header (@headers1) {
      if ( $header ne $IndexStr ) {
        print "$header\t";
      }
    }
    print "\n";
  }
  else {
    if ( $IndexCol == -1 ) {
      die "Did not find label $IndexStr in $Table2 - exiting.\n";
    }

    for (my $i=0; $i<=$#line; $i++) {
      print "$line[$i]\t";
    }
    foreach my $header (@headers1) {
      if ( $header ne $IndexStr ) {
        if ( exists($data1{$line[$IndexCol]}{$header}) ) { 
          print "$data1{$line[$IndexCol]}{$header}\t";
        } else {
          print "NA ";
        }
      }
    }
    print "\n";
  }

  $linecount++;
}
close(T2);

#print STDERR "read $linecount lines from $Table2\n";


