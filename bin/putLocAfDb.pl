#!/usr/bin/perl -w 
#
#use HTML::Extract;

use YAML::XS 'LoadFile';
use Sort::Key::Natural qw(natsort);


die "putLocAfDb.pl <miseq folder> \n" if (!(@ARGV));
die "putLocAfDb.pl <miseq folder \n" if ( $#ARGV != 0 );

chomp($miseqroot = $ARGV[0]);

open (FH, "</home/data_in/calls2xls/calls2xls.cfg");
while ( $line = <FH> ) {
    if($line !~ /#/ && $line =~ /\S/){
	@tmp = split(/=/, $line);
	foreach $row (@tmp) {
	    $row =~ s/\s//g;
	}
	$c{$tmp[0]}=$tmp[1];
	if($tmp[0] eq "info_dp"){
	    @info = split(/,/,$tmp[1]);
	    foreach $i (@info){
		$c{$tmp[0]}{$i} = 1;
	    }
	}
	if($tmp[0] eq "info_ao"){
	    @info = split(/,/,$tmp[1]);
	    foreach $i (@info){
		$c{$tmp[0]}{$i} = 1;
	    }
	}
    }
}

$rundate = (split /_/, $miseqroot)[0];

$finalpath = $miseqroot . "/calling/final/";

opendir(DIR, $finalpath);
@files = grep { /^project_/ } readdir(DIR);
closedir(DIR);

$projectpath = $finalpath . $files[0];
$summary_file = $projectpath . "/project-summary.yaml";

$config = LoadFile( $summary_file );
$dt = $config -> {date};
$dt =~ /\d\d(\d\d)-(\d\d)-(\d\d)\s/;
$calldate = $1 . $2 . $3;

$ensvcfgz = $projectpath . "/batch1-ensemble.vcf.gz";
$outvar = $c{'tmpdir'} . "/tmp.vcf";
$outnorm = $outvar . ".n.vcf";
$vardb_nogz = $c{'tmpdir'} . "/tmp_vardb.vcf";
$outdbbed = $c{'tmpdir'} . "/tmp_projectSampleDbCovs.s.bed.gz";
$projectsamplecovs = $c{'tmpdir'} . "/tmp_projectSampleDbCovs.bed";

open FH, "gunzip -c $ensvcfgz | ";
open(FH2, ">$outvar");

while ($line = <FH>) {
    if($line =~ /^#CHROM/) {
	chomp($line);
	@h = split(/\t/, $line);
	for($i = 9; $i <= $#h; $i++){
	    if($h[$i] =~ /\S/){
		$s2c{$h[$i]}=$i;
	    }
	}
	@header = splice(@h,0,8);
	$hline = join("\t", @header);
	print FH2 "$hline\n";
    }
    elsif($line =~ /^#/){
	print FH2 "$line";
    }
    else{
	chomp($line);
	@r= split(/\t/, $line);
	@sampINFO = ();
	foreach $s (sort { $a cmp $b } keys(%s2c)) {
	    @tmp = split(/:/, $r[$s2c{$s}]);
	    $alleles = 0;
	    while ($tmp[0] =~ /(\d+)/g) {
		$alleles += $1;
	    }
	    if($tmp[0] =~ /[\|\/]/){
		$tot = 2;
	    }
	    else{
		$tot = 1;
	    }
	    if($alleles > 0) {
		push @sampINFO, $s . ":$rundate" . ":$calldate" . "=$alleles" . "/" . $tot;
	    }
	}
	if(exists ($sampINFO[0])){
	    $r[7] = join(';', @sampINFO);
	    @l = splice(@r,0,8);
	    $ol = join("\t", @l);
	    print FH2 "$ol\n";
	}
    }
}
close(FH);
close(FH2);

system("vt normalize $outvar -r $c{'genome'} -o $outnorm"); 

open FH, "gunzip -c $c{'vardb'} | ";
while ($line = <FH>) {
    if($line !~ /^#/ && $line =~ /\S/){
	chomp($line);
	@tmp = split(/\t/, $line);
	@infos = split(/;/, (split(/\t/, $line))[7]);
	foreach $info (@infos) {
	    $line =~ s/\s//g;
	    ($ikey, $val) = split(/=/, $info);
	    $dat{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}{$ikey} = $val;
	}
    }
}

close(FH);

open(FH, "<$outnorm");
while ($line = <FH>) {
    if($line !~ /^#/ && $line =~ /\S/){
	chomp($line);
	@tmp = split(/\t/, $line);
	@infos = split(/;/, (split(/\t/, $line))[7]);
	foreach $info (@infos) {
	    $line =~ s/\s//g;
	    ($ikey, $val) = split(/=/, $info);
	    $dat{$tmp[0]}{$tmp[1]}{$tmp[3]}{$tmp[4]}{$ikey} = $val;
	}
    }
}
close(FH);
open(OUT, ">$vardb_nogz");
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
foreach $chr (natsort (keys(%dat))) {
    foreach $p (sort { $a <=> $b } keys(%{$dat{$chr}})) {    
        foreach $ref (sort { $a cmp $b } keys(%{$dat{$chr}{$p}})) {    
            foreach $alt (sort { $a cmp $b } keys(%{$dat{$chr}{$p}{$ref}})) {
                @infos = ();
                foreach $val (sort { $a cmp $b } keys(%{$dat{$chr}{$p}{$ref}{$alt}})) {
                    push @infos, $val . "=" . $dat{$chr}{$p}{$ref}{$alt}{$val};
                }
                $row = $chr . "\t" .  $p . "\t" . "." . "\t" . $ref . "\t" . $alt . "\t" . "100" . "\t" . "PASS" . "\t" . join(';', @infos);
#		print "$vardb_nogz\n";
                print OUT "$row\n";
            }
        }
    }
}
close(OUT);

system("bgzip $vardb_nogz");
$vardb_gz = $vardb_nogz . ".gz";

@sample_folders = <$finalpath/*>;

foreach $folder (@sample_folders){
    if (-d $folder && $folder !~ /project/) {
	print "$folder\n";
	$ccov = $folder . "/sample.coverage.bed";
	if (-e $ccov) {
	    system("zgrep -hv ^\"#\" $ccov >> $projectsamplecovs");
	}
    }
}

print "$c{'covdb'} $projectsamplecovs\n";

system("zgrep -hv ^\"#\" $c{'covdb'} >> $projectsamplecovs");

system("cat $projectsamplecovs | LC_ALL=C sort -S 30G -k1,1 -k2,2n | LC_ALL=C uniq | bgzip > $outdbbed");

system("mv $outdbbed $c{'covdb'}");
system("tabix -p bed $c{'covdb'}");
system("mv $vardb_gz $c{'vardb'}");
system("tabix -p vcf $c{'vardb'}");

