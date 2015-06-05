#!/usr/bin/perl -w
#
#use HTML::Extract;

die "getfreq.pl <vcf> \n" if (!(@ARGV));
die "getfreq.pl <vcf>\n" if ( $#ARGV != 0 );

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
	@tmp = split(/\t/, $line);
	
	$ssvcf = $tmp[0] . ":" . $tmp[1] . "-" . $tmp[1];
	$ssbed = $tmp[0] . ":" . $tmp[1] . "-" . $tmp[1];
	
	unless (open(GFH, "tabix $covdb $ssbed 2>&1 |")) {
	    die "Can't spawn external command!";
	}
	chomp(@covout = <GFH>);
	
	close(GFH);
	
	unless (open(GFH2, "tabix $vardb $ssvcf 2>&1 |")) {
	    die "Can't spawn external command!";
	}
	chomp(@varout = <GFH2>);
	
	close(GFH2);  
	@freqA=(0,0);
	if(exists($varout[0]) && exists($covout[0])){
	    
#	    print "---\n";
#	    print "$ssvcf\n";
#	    
#	    print "$line\n";

	    foreach $row(@varout) {
		@varsplit = split(/\t/, $row);
		if((($tmp[0] eq $varsplit[0] && $tmp[1] eq $varsplit[1]) && $tmp[3] eq $varsplit[3]) && $tmp[4] eq $varsplit[4]){
#		    $foundInDb = 1;
		    %chash = ();
		    %vhash = ();
		    foreach $row(@covout){
			@cov = split(/\t/, $row);
			if($cov[3] >= 30) {
			    $chash{$cov[5]} = $cov[8];
			}
		    }
		    @gts = split(/;/, $varsplit[7]);
		    foreach $row(@gts){
			@gt = split(/[:=]/, $row);
			if($chash{$gt[0]}){
			    $vhash{$gt[0]} = $gt[3];
			}
		    }
		    
		    foreach $g (sort { $a cmp $b } keys(%vhash)) {
			$vhash{$g} =~ s/\|/\//g;
			if($chash{$g} eq "male" && $tmp[0] =~ /chrX/){
			    $vhash{$g} = "1/1";
			}
			elsif($chash{$g} eq "male" && $tmp[0] =~ /chrY/){
			    $vhash{$g} = "1/1";
			}
			elsif($chash{$g} eq "female" && $tmp[0] =~ /chrY/){
			    delete $vhash{$g};
			}
			elsif($tmp[0] =~ /chrM/){
			    $vhash{$g} = "1/1";
			}
			elsif($vhash{$g} =~ /(\d)/ && $vhash{$g} !~ /\//){
			    $vhash{$g} = $1 . "/" . $1;
			}
		    }
		    $no_alleles = 0;
		    $tot_alleles = 0;
		    foreach $g (sort { $a cmp $b } keys(%vhash)) {
			if(exists($chash{$g})){
			    $vhash{$g} =~ /(\d)[\/\|](\d)/;			
			    $no_alleles += $1;
			    $tot_alleles += $2;
			}
		    }
#		    print "$no_alleles / $tot_alleles\n";
		    if($tot_alleles != 0){
			$freqA[0] = sprintf("%.2f", $no_alleles / $tot_alleles);
			$freqA[1] = $no_alleles;
		    }
		}
	    }
	}
	$freqstr = ";LOCFREQ=" . $freqA[0] . ";LOCALS=" . $freqA[1];
	$tmp[7] .= $freqstr;
	$out = join("\t", @tmp);
	print "$out\n";
    }
}
