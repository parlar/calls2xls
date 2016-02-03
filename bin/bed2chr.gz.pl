#!/usr/bin/perl -w
#

die "bed2chr.gz.pl <in.bed> <out.path> \n" if (!(@ARGV));
die "bed2chr.gz.pl <in.bed> <out.path> \n" if ( $#ARGV != 1 );

chomp($bed = $ARGV[0]);
chomp($path = $ARGV[1]);

if(!-d $path){
   system("mkdir -p $path");
} 

open(FH, "<$bed");
   
while($line=<FH>){
	@brow = split(/\t/, $line);
	$addpath = $path . "/" . $brow[5] . "/" . $brow[6] . "/" . $brow[7];

	if(!-d $addpath){
	  system("mkdir -p $addpath");
	}
	$outfile = $addpath . "/" . $brow[0] . ".gz";
#    print "$outfile\n";
	if(!-e $outfile){
	  if(fileno(FH2)){
		 close(FH2);
	  }
	  open(FH2, "| sort-bed - | pigz -c > $outfile") or die "Couldn't open: $!";
   }
   else{
	  if(fileno(FH2)){
		 print FH2 $line;
	  }
	  else{
		 print "$outfile exists, filehandle not open, $bed $path\n";
		 last;
	  }
	}
}
close(FH2);
close(FH);

    
