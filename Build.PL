#!/usr/bin/perl

use strict;
use warnings;
use Module::Build;

my $build = Module::Build->new(
	build_class       => 'Module::Build',
	module_name       => 'Bio::DB::USeq',
	license           => 'perl',
	dist_version_from => 'lib/Bio/DB/USeq.pm',
	dist_author       => 'Timothy Parnell <parnell.tj@gmail.com>',
	configure_requires => {
		'Module::Build'       => 0,
	},
	requires          => {
		'perl'                => '5.10.0',
		'Bio::Root::Version'  => '1.0069001',
		'Archive::Zip'        => '1.30',
	},
);

$build->create_build_script;
exit 0;
