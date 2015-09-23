#!/usr/bin/perl -w
#
use warnings;

die "read_bcbio_yaml.pl <yaml file> <sample>\n" if (!(@ARGV));
die "read_bcbio_yaml.pl <yaml file> <sample>\n" if ( $#ARGV != 1 );    

chomp($file = $ARGV[0]);
chomp($sample = $ARGV[1]);

open(FH, "<$file");
chomp(@YML=<FH>);
close(FH);

$keep1 = 0;
@SS = ();
$cs = "EMTY";

foreach $row (@YML) {
    if($row =~ /^\s+-\s+/){
	if(exists($SS[0]) && $cs ne "EMTY"){
	    for ($i = 0;  $i <= $#SS; $i++){
		push @{$configs{$cs}}, $SS[$i];
	    }
	}
	@SS = ();
	push @SS, $row;
	$cs = "EMTY";
	$keep1 = 1;
	next;
    }
    if($row =~/^\w/) {
	if(exists($SS[0]) && $cs ne "EMTY"){
	    for ($i = 0;  $i <= $#SS; $i++){
		push @{$configs{$cs}}, $SS[$i];
	    }
	}
	@SS = ();
	$cs = "EMTY";
	$keep1 = 0;
    }
    if($keep1) {
	push @SS, $row;
    }   
    if($row =~/^\s+description:\s(\S+)/) {
	$cs = $1;
    }
}

if(exists($SS[0]) && $cs ne "EMTY"){
    for ($i = 0;  $i <= $#SS; $i++){
	push @{$configs{$cs}}, $SS[$i];
    }
}

foreach $row (@{$configs{$sample}}){
    print "$row\n";
}
