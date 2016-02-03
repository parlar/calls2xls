#!/usr/bin/perl -w
#use strict;
use Spreadsheet::ParseExcel;
use Encode qw(encode decode);
use Config::General;
use Getopt::Long;
use Time::localtime;

GetOptions (\%o, "help");

#$enc = 'utf-8';

if(exists($o{'help'})){
    print <<"::_USAGE0_BLOCK_END_::";
        
  usage: mkConfigOptionsFile.pl [options]
    options:         
      -h, --help    Display this help message
      
  The script auto-generates a config options file in folder 
  defined in calls2xls.cfg
      
::_USAGE0_BLOCK_END_::
   
    exit 1
}

## Read configs
my $conf = Config::General->new("$ENV{'CALLS2XLS'}/calls2xls.cfg");
my %c = $conf->getall;


$nowdate = &get_nowtime();

$designpath = "$ENV{'CALLS2XLS'}/" . $c{'designAndGeneListsRootPath'};
@designfiles = &get_files($designpath);

$diseasepath = "$ENV{'CALLS2XLS'}/" . $c{'diseaseGeneAssocPath'};
@diseasefiles = &get_files($diseasepath);

#print "@designfiles\n";
#print "@diseasefiles\n";


$outfile = $c{'miseqOptionsFilePath'} . "/" . "config." . $nowdate . ".txt"; 

open(FH, ">$outfile");

print FH "\nSampleSheet.csv template:\n";
print FH "   Design\$Gender\$Disease_hypothesis\$GeneListFile\n\n";
print FH "Design file:\n";

foreach $file (@designfiles){
    print FH "   $file\n";
}

print FH "Gender:\n";
print FH "   male/female\n";

print FH "\nGeneListFiles:\n";
foreach $file (@diseasefiles){
    print FH "   $file\n";
    print FH "      ";

    open(FH2, "<$diseasepath/$file");
    chomp(@ddata = <FH2>);
    close(FH2);

    foreach $row (@ddata){
	if($row =~ /^>(\S+)/){
	    print FH "$1 "
	}
    }
    print FH "\n\n";
}

close(FH);

sub get_nowtime(){
	my $tm = localtime;
	return sprintf("%04d.%02d.%02d-%02d.%02d.%02d", $tm->year+1900, ($tm->mon)+1, $tm->mday, $tm->hour, $tm->min, $tm->sec );
}

sub get_files(){
	my $finalpath = shift;
	my @files;
	opendir(DIR, $finalpath);
	@files = grep {/bed$|txt$/} readdir(DIR);
	close(DIR);
	return sort {$b cmp $a} @files;
}
