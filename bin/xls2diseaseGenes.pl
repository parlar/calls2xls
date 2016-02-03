#!/usr/bin/perl -w
#use strict;
use Spreadsheet::ParseExcel;
use Encode qw(encode decode);
use Config::General;

die "xls2disease_genes.pl <xls> <date>\n" if (!(@ARGV));
die "xls2disease_genes.pl <xls> <date>\n" if ( $#ARGV != 1 );

chomp($xls = $ARGV[0]);
chomp($date = $ARGV[1]);

$enc = 'utf-8';


## Read configs
my $conf = Config::General->new("$ENV{'CALLS2XLS'}/calls2xls.cfg");
my %c = $conf->getall;


my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->parse($xls);

if ( !defined $workbook ) {
    die $parser->error(), ".\n";
}

$calls2xls_root = $ENV{'CALLS2XLS'};

## Read configs
open (FH, "<$calls2xls_root/calls2xls.cfg");
while ( $line = <FH> ) {
    if($line !~ /#/ && $line =~ /\S/){
        @tmp = split(/\s=\s/, $line);
        foreach $row (@tmp) {
            $row =~ s/\s//g;
        }
        if($tmp[0] eq "info_dp"){
            @info = split(/,/,$tmp[1]);
            foreach $i (@info){
                $c{$tmp[0]}{$i} = 1;
            }
        }
        elsif($tmp[0] eq "info_ao"){
            @info = split(/,/,$tmp[1]);
            foreach $i (@info){
                $c{$tmp[0]}{$i} = 1;
            }
        }
        else{
            $c{$tmp[0]}=$tmp[1];
        }
    }
}

#$date = ;

$output = $calls2xls_root . "/" . $c{'diseaseGeneAssocPath'} . "/" . "gl." . $date . ".txt"; 

open(FH, ">$output");

for my $worksheet ( $workbook->worksheets() ) {
    $name = $worksheet->get_name();
    $name_uni = encode($enc, $name);
    if($name_uni !~ /versikt/ && $name_uni !~ /flikar/){
	print FH ">$name_uni\n";    
	my ( $row_min, $row_max ) = $worksheet->row_range();
	for my $row ( 1 .. $row_max ) {
	    @genedat = ();
	    for my $col ( 0 .. 2 ) {
		my $cell = $worksheet->get_cell( $row, $col );
		if($cell){
		    $cv = $cell->value();
		    $cv =~ s/\s//g;
		    push @genedat, $cv;
		}
		else{
		    push @genedat, "NA";
		}
	    }
	    $genestr = join("\t", @genedat);
	    if( $genestr !~ /^NA\t/){
		print FH "$genestr\n";
	    }
	}
    }
}

close(FH);

system("mkConfigOptionsFile.pl");
