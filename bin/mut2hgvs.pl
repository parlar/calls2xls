#!/usr/bin/perl -w
#
use strict;
#use Time::localtime;
#use List::Util qw(max min);
#use Config::General;
use XML::LibXML;
#use Digest::MD5  qw(md5_hex);
#use Sort::Key::Natural;

die "read_mut.pl <pf> \n" if (!(@ARGV));
die "read_mut.pl <pf> \n" if ( $#ARGV != 0 );

chomp(my $pf = $ARGV[0]);
chomp(my @files = <*.$pf>);

#my $nowtime = get_nowtime();
#my $md5 = md5_hex($nowtime);
my %alldat;
foreach my $file (@files){
	my $parser = XML::LibXML->new();
	my $doc    = $parser->parse_file($file);

	#print "$mut\n";
	foreach my $mutation_node ($doc->findnodes('/Mutations/Mutation')) {
		my ($id) = $mutation_node->getAttribute('id');
		my ($assembly) = $mutation_node->getAttribute('refAssembly');
		my ($chr) = $mutation_node->getAttribute('chr');
		my ($variant_node) = $mutation_node->findnodes('./Variant');
		my ($type) = $variant_node->getAttribute('type');
		my ($nomen_node) = $variant_node->findnodes('./Nomenclature');
		my ($refseq) = $nomen_node->getAttribute('refSeq');
		my ($cnomen_node) = $nomen_node->findnodes('./cNomen');
		my ($hgvs) = $cnomen_node->getAttribute('val');
		$assembly =~ s/\s//g;

		$alldat{$assembly}{$id}{'refseq'} = $refseq;
		$alldat{$assembly}{$id}{'hgvs'} = $hgvs;
#		print "$id\t$assembly\t$chr\t$refseq\t$hgvs\n";
	}
}

foreach my $assembly (sort { $a cmp $b } keys (%alldat)) {
	open (FH, ">/home/data_in/tmp/$assembly.txt");
	foreach my $id (sort { $a cmp $b } keys (%{$alldat{$assembly}})) {
		print FH $alldat{$assembly}{$id}{'refseq'} . ":" . $alldat{$assembly}{$id}{'hgvs'} . "\n";
	}
	close(FH);
}
		
	
	
