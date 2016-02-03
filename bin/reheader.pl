#!/usr/bin/perl -w
#

die "reheader.pl <infile.vcf.gz> <outfile.vcf.gz>\n" if (!(@ARGV));
die "reheader.pl <infile.vcf.gz> <outfile.vcf.gz>\n" if ( $#ARGV != 1 );

chomp($infile = $ARGV[0]);
chomp($outfile = $ARGV[1]);

$head = $infile;
$head =~ s/.vcf.gz$/.nhead.txt/g;

@header = ();
open FH, "gunzip -c $infile | ";
while ($line = <FH>) {
    if($line =~ /^#/){
		push @header, $line;
    }
    else{
		$info = (split(/\s/, $line))[7];
		@entries = split(/;/, $info);
		foreach $row (@entries){
			if($row =~ /=/){
				$key = (split(/=/, $row))[0];
				$keys{$key} = 1;
			}
		}
    }
}
close(FH);
open FH, ">$head";
print FH "$header[0]";
foreach $key ( sort { $b cmp $a } keys(%keys)) {
    print FH "##INFO=<ID=" . $key . ",Number=A,Type=String,Description=\"sample:rundate:calldate=foundAlleles/totalAlleles.\">\n";
}
print FH "$header[-1]";

close(FH);
#print "bcftools reheader -h $head -o $outfile $infile\n";
system("bcftools reheader -h $head -o $outfile $infile");
system("tabix -fp vcf $outfile")
