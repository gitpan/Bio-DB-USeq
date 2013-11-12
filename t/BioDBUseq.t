#!/usr/bin/perl

# Test script for Bio::DB::USeq

use strict;
use Test;
use FindBin '$Bin';

BEGIN {
	plan tests => 20,
}

use lib "$Bin/../lib";
use Bio::DB::USeq qw(binMean binStdev);

my $file = "$Bin/data/sample.useq";

# initialize and open file
my $useq = Bio::DB::USeq->new($file);
ok($useq);
die "unable to open useq file!\n" unless $useq;

# check metadata parameters of useq file
ok($useq->type, 'region');

# check seq_ids
my @seq_ids = $useq->seq_ids;
ok(scalar @seq_ids, 2);
ok($useq->length('chrI'), 230042);

# check segment
my $segment = $useq->segment(
	-seq_id     => 'chrI',
	-start      => 166301,
	-end        => 166500,
);
ok($segment);
ok($segment->start, 166301);

# check scores
my @scores = $segment->scores;
# print "These are ", scalar(@scores), " segment scores: ", join(' ', @scores), "\n";
ok(scalar @scores, 31);
ok(sprintf("%.2f", $scores[0]), 7.91); # round to 2 decimal places to avoid any
ok(sprintf("%.2f", $scores[2]), 8.02); # potential floating point errors

# check features via an iterator
my $i = $segment->get_seq_stream;
ok($i);
my $f = $i->next_seq;
ok($f);
# print "f end is ", $f->end, "\n";
ok($f->end, 166302);

# check wiggle function via bins and across a member junction
my ($wiggle) = $useq->features(
	-seq_id     => 'chrII',
	-start      => 455001,
	-end        => 475000,
	-type       => 'wiggle:1000',
);
ok($wiggle);
my @wig_bins = $wiggle->wiggle;
# print "There are ", scalar(@wig_bins), " wig bins: ", join(' ', @wig_bins), "\n";
ok(scalar(@wig_bins), 1000);
ok(sprintf("%.2f", $wig_bins[9]), 1.94);

# check chromosome stats
my $chr_mean = $useq->chr_mean('chrI');
# print "The chromosome mean for chrI is $chr_mean\n";
ok(sprintf("%.2f", $chr_mean), 1.88);

# check statistical function
my ($summary) = $useq->features(
	-seq_id     => 'chrII',
	-start      => 1,
	-end        => 500000,
	-type       => 'summary',
);
ok($summary);
my $stats = $summary->statistical_summary(10);
ok(scalar @$stats, 10);
my $second_mean  = binMean( $stats->[1] );
my $second_stdev = binStdev( $stats->[1] );
# print "Second interval mean is $second_mean\n";
# print "Second interval stdev is $second_stdev\n"; 
ok(sprintf("%.2f", $second_mean), '0.60');
ok(sprintf("%.2f", $second_stdev), 0.62);


exit 0;
