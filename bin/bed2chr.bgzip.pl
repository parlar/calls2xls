#!/usr/bin/perl -w
#

use DBI;

die "bed2chr-starch.pl <in.bed> <out.path> \n" if (!(@ARGV));
die "bed2chr-starch.pl <in.bed> <out.path> \n" if ( $#ARGV != 1 );

chomp($bed = $ARGV[0]);
chomp($path = $ARGV[1]);

if(!-d $path){
   system("mkdir -p $path");
} 

open(FH, "<$bed");
while($line=<FH>){
    @brow = split(/\t/, $line);
    $outfile = $path . "/" . $brow[0] . ".gz";
#    print "$outfile\n";
    if(!-e $outfile){
	if(fileno(FH2)){
	    close(FH2);
	}
	open(FH2, "| sort-bed - | bgzip -c > $outfile") or die "Couldn't open: $!";
    }
    else{
	if(fileno(FH2)){
	    print FH2 $line;
	}
	else{
	    print "$bed $path\n";
	    die;
	}
    }
}
close(FH2);
close(FH);

    
