#!/usr/bin/perl -w
#

use Cwd;

die "make_fastq_links.pl <folder> \n" if (!(@ARGV));
die "make_fastq_links.pl <folder> \n" if ( $#ARGV != 0 );

open (FH, "</home/data_in/calls2xls/calls2xls.cfg");

while ( $line = <FH> ) {
    if($line !~ /#/ && $line =~ /\S/){
        @tmp = split(/=/, $line);
        foreach $row (@tmp) {
            $row =~ s/\s//g;
        }
        $c{$tmp[0]}=$tmp[1];
    }
}
close(FH);



chomp($root = $ARGV[0]);

$search = "*." . "fastq.gz";

unless (open(FH0, "find $root -name \"$search\" 2>&1 |")) {
    die "Can't spawn external command!";
}

chomp(@FQ = <FH0>);
unless (close(FH0)) {  
    die "External command failed: $?";
}

if($root !~ /\/$/){
    $root = $root . "/";
}
$currdir = getcwd;
$newdir = $root . "fastqs";
$newlndir = $root . "fastqs-lns";

system("mkdir $newdir");
system("mkdir $newlndir");

foreach $file (@FQ) {
    @tmp = split(/\//, $file);
    @tmp2 = split(/_/, $tmp[-1]);
    $longfile = $currdir . "/" . $file;
    $sampleFiles{$tmp2[0]}{$longfile} = 1;
}

open(FH1, ">commands1.txt");

foreach $sample (sort { $a cmp $b } keys(%sampleFiles)) {
    unless ($sample =~ /Undetermined/){
	@tmpFiles = ();
	$f1 = "";
	$f2 = "";
	foreach $file (sort { $a cmp $b } keys(%{$sampleFiles{$sample}})) {
	    push @tmpFiles, $file;
	}
	foreach $file (@tmpFiles) {
	    if ($file =~ /_R1_/) {
		$f1 = $file;
	    }
	    if ($file =~ /_R2_/) {
		$f2 = $file;
	    }
	}
	print FH1 "$c{'trim_galore'} --length 44 --paired --clip_R2 5 --clip_R1 5 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT $f1 $f2\n";
    }
}

close(FH1);

$currdir = getcwd;
print "$currdir $newdir\n";
system("mv commands1.txt $newdir");
chdir $newdir;
system("parallel -j 16 -- < commands1.txt");
system("rm commands1.txt");

@fqs = <*.fq.gz>;
$fasdir = getcwd;
foreach $file (@fqs) {
    $short = $file;
    $short =~ s/_S\d+.*_val_/_/;
    $pathshort = "../fastqs-lns/" . $short;
    $target = $fasdir . "/" . $file;
    print "ln -s $target $pathshort \n";
    system("ln -s $target $pathshort");
}
chdir $currdir;
