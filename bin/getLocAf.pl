#!/usr/bin/perl -w
#
#use HTML::Extract;

die "getLocAf.pl <vcf> \n" if (!(@ARGV));
die "getLocAf.pl <vcf>\n" if ( $#ARGV != 0 );

chomp($vcf = $ARGV[0]);

$covdb = "/home/data_in/databases/covVarDB/covdb.bed.gz";
$vardb = "/home/data_in/databases/covVarDB/vardb.vcf.gz";


open (FH, "$vcf");
while ($line = <FH>) {
    if ($line =~ /^#/ &&  $line =~ /\S/){
	print "$line";
    }
    if ($line !~ /^#/ &&  $line =~ /\S/){
	chomp($line);
#	print "---\n$line";
	@vcfRow = split(/\t/, $line);
	
	$position = $vcfRow[0] . ":" . $vcfRow[1] . "-" . $vcfRow[1];
	
	unless (open(GFH, "tabix $covdb $position 2>&1 |")) {
	    die "Can't spawn external command!";
	}
	chomp(@covout = <GFH>);
	close(GFH);

	$nototalleles = 0;
	$nofoundalleles = 0;
	
	%coveredSamples = ();

	foreach $row(@covout){
	    @cov = split(/\t/, $row);
	    if($cov[3] >= 30) {
		if($cov[0] eq "chrM"){
		    $nototalleles += 1;
		}
		elsif($cov[0] eq "chrY" && $cov[8] eq "male"){
		    $nototalleles += 1;
		}
		elsif($cov[0] eq "chrX" && $cov[8] eq "female"){
		    $nototalleles += 2;
		}
		elsif($cov[0] eq "chrX" && $cov[8] eq "male"){
		    $nototalleles += 1;
		}
		else{
		    $nototalleles += 2;
		}
		$id = $cov[5] . ":" . $cov[6] . ":" . $cov[7];
		$coveredSamples{$id} = 1;
	    }
	}
	
	unless (open(GFH2, "tabix $vardb $position 2>&1 |")) {
	    die "Can't spawn external command!";
	}
	chomp(@varout = <GFH2>);
	
	close(GFH2);
	
	if(exists($varout[0]) && exists($covout[0])){
	    foreach $row(@varout) {
		@varsplit = split(/\t/, $row);
		if((($vcfRow[0] eq $varsplit[0] && $vcfRow[1] eq $vcfRow[1]) && $vcfRow[3] eq $varsplit[3]) && $vcfRow[4] eq $varsplit[4]){

		    @samples = split(/;/, $varsplit[7]);

		    foreach $row2(@samples){
			($key, $value) = split(/=/, $row2);
			if(exists($coveredSamples{$key})){
			    @tmp = split(/[\/\|]/, $value);
			    $nofoundalleles += $tmp[0];
			}
		    }
		}
	    }
	}
	$freq = sprintf("%.5f", $nofoundalleles / $nototalleles);

	$freqstr = ";LOC_AF=" . $freq . ";LOC_AC=" . $nofoundalleles . ";LOC_AN=" . $nototalleles;
	$vcfRow[7] .= $freqstr;
	$out = join("\t", @vcfRow);
	print "$out\n";
    }
}
