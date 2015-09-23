#!/usr/bin/perl -w
#
use warnings;
use Getopt::Long;
use YAML::XS 'LoadFile';

GetOptions (\%o, "summary_file=s", "output_type=s", "samples=s");

sub helpmess0{
    print "usage: read_bcbio_summary_yaml.pl [options]\n";
    print "  --summary_file=file          (required)\n";     
    print "             Pipeline summary output file\n";
    print "  --output_type=mapdata|date   (required)\n"; 
    print "             What to output\n";
    print "  --sample=sample_number(s)    (requred for mapdata)\n";
    print "             Ex. 15-1234,15-1235\n";
}

if(!(exists($o{'summary_file'}) && exists($o{'output_type'}))){
    &helpmess0();
    exit 0;
}

if(exists($o{'summary_file'}) && exists($o{'output_type'})){
    if($o{'output_type'} !~ /date|mapdata/ && !exists($o{'samples'})){
	&helpmess0();
	print "hopp\n";
	exit 0;
    }
}

$config = LoadFile( $o{'summary_file'} );

if($o{'output_type'} eq "date"){
    $dt = $config -> {date};
    $dt =~ /\d\d(\d\d)-(\d\d)-(\d\d)\s/;
    $calldate = $1 . $2 . $3;
    print "$calldate\n";
}
elsif($o{'output_type'} eq "mapdata"){
    @samples = split(/,/, $o{'samples'});
    for($i = 0;; $i++){
	if(exists($config->{samples}->[$i])){
	    foreach $sample (@samples) {
		if($config->{samples}->[$i]->{description} eq $sample){
		    print "$sample:\n";
		    foreach $key (sort keys(%{$config->{samples}->[$i]->{summary}->{metrics}})){
			print "  $key : $config->{samples}->[$i]->{summary}->{metrics}->{$key}\n";
		    }
		}
	    }
	}
	else{
	    last;
	}
    }
}

