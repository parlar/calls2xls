#!/usr/bin/perl -w
#

use Config::General;
use Cwd;
use MCE::Map
      max_workers =>16,             ## Default 'auto'
      chunk_size => 1              ## Default 'auto'
;

die "makeFastqLinks.pl <runfolder> \n" if (!(@ARGV));
die "makeFastqLinks.pl <runfolder> \n" if ( $#ARGV != 0 );

$conf = Config::General->new("$ENV{'CALLS2XLS'}/calls2xls.cfg");
%c = $conf->getall;

sub run_trimgalore{
    ($f1, $f2) = split(/;/, $_);
    system("$ENV{'CALLS2XLS'}/$c{'trim_galore'} --length 40 --paired --clip_R2 5 --clip_R1 5 --output_dir $trimOutPath -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT $f1 $f2");
    return(1);
}

#my $currdir = getcwd;

chomp($runfolder = $ARGV[0]);
$runfolderDataPath = $runfolder . "/" . $c{'runFolderDataPath'};

$trimOutPath = $runfolder . "/fastqs";
$trimOutLnPath = $runfolder . "/fastqs-lns";

system("mkdir -p $trimOutPath");
system("mkdir -p $trimOutLnPath");

print "$runfolderDataPath/*.fastq.gz\n";

chomp(@GZ = <$runfolderDataPath/*.fastq.gz>);

foreach $file (@GZ) {
    if($file =~ /.gz$/){
	$fastqfile = (split(/\//, $file))[-1];
	$sample = (split(/_/, $fastqfile))[0];
	$sampleFiles{$sample}{$file} = 1;
    }
}

@cmds = ();

foreach $sample (sort { $a cmp $b } keys(%sampleFiles)) {
    unless ($sample =~ /Undetermined/){
	$f1 = "";
	$f2 = "";
	foreach $fastqfile (sort { $a cmp $b } keys(%{$sampleFiles{$sample}})) {
	    print "$fastqfile\n";
	    if ($fastqfile =~ /_R1_/) {
		$f1 = $fastqfile;
	    }
	    if ($fastqfile =~ /_R2_/) {
		$f2 = $fastqfile;
	    }
	}
	push @cmds, $f1 . ";" . $f2
    }
}
foreach $row (@cmds){
    print "$row\n";
}

@out = mce_map { run_trimgalore } @cmds;
print "@out\n";
$workdir = getcwd;

@trimmedFqs = <$trimOutPath/*.fq.gz>;

foreach $longfile (@trimmedFqs) {
    $file = (split(/\//, $longfile))[-1];
    $sample = (split(/_/, $file))[0];
    if($file =~ /^\S+.*_val(_[12].fq.gz)$/){
	print "$sample $1 \n";
	$linkname = $workdir . "/" . $trimOutLnPath . "/" . $sample . $1;
	$target = $workdir . "/" . $longfile;
	print "ln -sr $target $linkname \n";
	system("ln -sr $target $linkname");
    }
}
