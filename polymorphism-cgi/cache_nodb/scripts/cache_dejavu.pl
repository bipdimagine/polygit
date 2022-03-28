#!/usr/bin/perl
use strict;
use FindBin qw($RealBin);
use Getopt::Long;
use Data::Dumper;
use Parallel::ForkManager;
use lib "$RealBin/../../../GenBo/lib/obj-nodb/";
use GBuffer;

my $fork = 1;
my ($project_name, $chr_name);
GetOptions(
	'project=s' => \$project_name,
	'chr=s'     => \$chr_name,
);
die("\n\nERROR: -project option not defined. Die...\n\n") unless ($project_name);
die("\n\nERROR: -chr option not defined. Die...\n\n") unless ($chr_name);



my $hh;
my $buffer = new GBuffer;
$buffer->vmtouch(1);
my $project = $buffer->newProjectCache( -name => $project_name );
my $chr = $project->getChromosome($chr_name);
my $nb_var_chr = $chr->countThisVariants($chr->getVariantsVector());
if ($nb_var_chr > 0) {
	foreach my $var (@{$chr->getListVarObjects($chr->getVariantsVector())}) {
		$hh->{$var->id()} = $var->{heho_string};
	}
	if ($hh) {
		my $nb_var_hh = scalar (keys %{$hh});
		if ($nb_var_hh < $nb_var_chr) {
			die("\n\nERROR: cache_dejavu died for $project_name in CHR$chr_name. (Nb Var OBS: $nb_var_hh - Exp: $nb_var_chr). Die...\n\n");
		}
		else {
			#store( $hh, $project->lmdb_cache_dir."/$chr_name.dv.freeze" );
			warn Dumper $hh;
			warn $project->lmdb_cache_dir."/$chr_name.dv.freeze";
			warn "Nb Var OBS: $nb_var_hh - Exp: $nb_var_chr"; 
		}
	}
	else {
		die("\n\nERROR: cache_dejavu died for $project_name in CHR$chr_name. (Nb Var OBS: 0 - Exp: $nb_var_chr). Die...\n\n");
	}
}
