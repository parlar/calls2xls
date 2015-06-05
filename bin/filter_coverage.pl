#!/usr/bin/perl -w
#
use List::Util qw(max min);
die "filter_coverage.pl <bed>\n" if (!(@ARGV));
die "filter_coverage.pl <bed>\n" if ( $#ARGV != 0 );

chomp($d1 = $ARGV[0]);

open(FH0, "<$d1");
chomp(@FQ = <FH0>);
close(FH0);  

foreach $row (@FQ) {
    @tmp = split(/\t/, $row);
    if($tmp[3] < 30) {
	print "$row\n";
    }
}
