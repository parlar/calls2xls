#!/usr/bin/perl -w
#


die "remEOFbgzip.pl <file> \n" if (!(@ARGV));
die "remEOFbgzip.pl <file> \n" if ( $#ARGV != 0 );

$refEOF = "1f 8b 08 04 00 00 00 00 00 ff 06 00 42 43 02 00 1b 00 03 00 00 00 00 00 00 00 00 00 ";

chomp($file = $ARGV[0]);

$fsize = (stat($file))[7];

open FILE, $file or die "Error message here: $!";

binmode FILE;

my $data; # read buffer

$ok = 0;

$EOF = "";
seek(FILE,-28,2);
while (read FILE, $data, 1) {
    $EOF .= unpack("H*",$data) . " ";
}
if($refEOF eq $EOF){
    print "$file ok EOF";
    $ok = 1;
}
else{
    print "$file incorrect EOF - will not be truncated\n";
}

print "\n";
close(FILE);

if($ok){
    $newsize = $fsize - 28;
    truncate($file, $newsize);
}
