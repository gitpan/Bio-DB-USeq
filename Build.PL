#!/usr/bin/env perl

use strict;
use warnings;
use Module::Build;

my $build = Module::Build->new(
	build_class       => 'Module::Build',
	module_name       => 'Bio::DB::USeq',
	license           => 'perl',
	dist_version_from => 'lib/Bio/DB/USeq.pm',
	dist_author       => 'Timothy Parnell <parnell.tj@gmail.com>',
	requires          => {
		'Bio::Root::Version'  => '1.0069001',
		'Archive::Zip'        => '1.30',
	},
);

$build->create_build_script;
exit 0;
