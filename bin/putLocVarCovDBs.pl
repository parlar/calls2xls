#!/usr/bin/perl -w 
#
#use HTML::Extract;

use YAML::XS 'LoadFile';
use Sort::Key::Natural qw(natsort);
use Config::General;
use File::HomeDir;

die "putLocVarCovDBs.pl <miseq folder> \n" if (!(@ARGV));
die "putLocVarCovDBs.pl <miseq folder> \n" if ( $#ARGV != 0 );

chomp($miseqroot = $ARGV[0]);

## Read configs
$conf = Config::General->new("$ENV{'CALLS2XLS'}/calls2xls.cfg");
%c = $conf->getall;

$rundate    = (split /_/, $miseqroot)[0];
$finalpath  = $miseqroot . "/calling/final/";
$home       = File::HomeDir->my_home;


#print "$finalpath\n";

opendir(DIR, $finalpath);
chomp(@filesfolders = grep(!/^\.$|^\.\.$/,readdir(DIR)));
closedir(DIR);

foreach $row (@filesfolders){
    if (-d $finalpath . "/" . $row){
        if ($row =~ /^project_/) {
            $projectf = $row;
        }
        else{
            push @sample_folders, $row;
        }
    }
}

$projectpath = $finalpath . $projectf;
$summary_file = $projectpath . "/project-summary.yaml";

$config = LoadFile( $summary_file );
$dt = $config -> {date};
$dt =~ /\d\d(\d\d)-(\d\d)-(\d\d)\s/;
$calldate = $1 . $2 . $3;

$ensvcfgz = $projectpath . "/batch1-ensemble.vcf.gz";
$outvar = $c{'tmpdir'} . "/tmp.dc.n.vcf.gz";
$outvar_mod = $c{'tmpdir'} . "/tmp.dc.n.mod.vcf.gz";
$reheaded = $c{'tmpdir'} . "/tmp.dc.n.mod.rh.vcf.gz";

#read and restructure input vcf

system("vt decompose -s $ensvcfgz | vt normalize -r $c{'genome'} - | vt sort - | vt uniq - | bgzip -c >$outvar"); 

open FH2, '| pbgzip -c > ' . $outvar_mod;
open FH, "pbgzip -d -c $outvar | ";
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
			$tmp[0] =~ s/\./0/g;
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

print "$outvar_mod $reheaded\n";


system("tabix -f -p vcf $outvar_mod");
#print "tabix -fp vcf $reheaded\n";
system("reheader.pl $outvar_mod $reheaded");
#print "tabix -fp vcf $reheaded\n";
system("tabix -f -p vcf $reheaded");

#print "tabix -f -p vcf $reheaded\n";

## Merge vcfs

$combined_file = $reheaded . ".comb.vcf.gz";


$vardb = $c{'vardb'};
$vardb =~ s/~/$home/g;


print "$vardb\n";

if(-e $vardb){
#    print "vardb exists\n";
	system("bcftools merge -m none --output-type z $c{'vardb'} $reheaded >$combined_file");
    print "bcftools merge -m none --output-type z $c{'vardb'} $reheaded >$combined_file\n";
}
else{
	system("cp $reheaded $combined_file");
}

#copy, replace vardb 

system("cp $c{'vardb'} $c{'vardb'}.old");
system("cp $combined_file $c{'vardb'}");
system("tabix -f -p vcf $c{'vardb'}");


#copy bedgraph files to coverage db

foreach $folder (@sample_folders){
    print "$folder\n";
    if (-d $finalpath . "/" . $folder) {
		print "is folder $folder\n";
		$ccov = $finalpath . "/" . $folder . "/" . "sample.coverage.bed";
		$covd_in = $c{'_cov'} . "/covdata";  
		if (-e $ccov) {
			system("bed2chr.gz.pl $ccov $covd_in");
		}
    }
}

unless (open(FH0, "find $c('_cov') -name \"*.bed.gz\" 2>&1 |")) {
    die "Can't spawn external command!";
}
chomp(@COVS = <FH0>);
unless (close(FH0)) {  
    die "External command failed: $?";
}

foreach $cov (@COVS){
    print "$cov\n";
}

=pod

print "$c{'covdb'} $projectsamplecovs\n";

system("zgrep -hv ^\"#\" $c{'covdb'} >> $projectsamplecovs");

system("cat $projectsamplecovs | LC_ALL=C sort -S 30G -k1,1 -k2,2n | LC_ALL=C uniq | bgzip > $outdbbed");

system("mv $outdbbed $c{'covdb'}");
system("tabix -p bed $c{'covdb'}");
system("mv $combined_file $c{'vardb'}");
system("tabix -p vcf $c{'vardb'}");

#system("makeFreqDB.pl $c{'vardb'} $c{'covdb'} $c{'freqdb'});

