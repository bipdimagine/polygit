package VcfMerge;

use strict;
use Moo;

use warnings;
use Carp;
use FindBin qw($Bin);
use Vcf;
use Data::Dumper;


has tabix => (
	is		=> 'ro',
	required=> 1,
);

has list_VCF => (
	is		=> 'ro',
	required=> 1,
);

has VCF_out => (
	is		=> 'ro',
	required=> 1,
);

has remove_duplicates => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> undef,
);

has ref_for_missing => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> undef,
);

has trim_ALTs => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> 0,
);

has vcf_header => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> undef,
);

has collapse => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> undef,
);

has regions => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> undef,
);

has silent => (
	is		=> 'ro',
	lazy	=> 1,
	default	=> undef,
);


sub parse_params {
	my $self = shift;
	my $opts = 
    { 
        args => [$0, @ARGV],
        joiner => 
        {
            'I16' => \&joiner_dp4,
            'DP4' => \&joiner_dp4,
            'MQ0'  => \&joiner_sum,
            'DP'  => \&joiner_sum,
        },
        trim_redundant_ALTs => 0,
        collapse_any => 1,
    };
    $$opts{tabix_path} = $self->tabix();
    $$opts{files} = $self->list_VCF();
    $$opts{outfile} = $self->VCF_out();
    $$opts{rm_dups} = 1 if ($self->remove_duplicates());
    $$opts{ref_for_missing} = $self->ref_for_missing() if ($self->ref_for_missing());
    $$opts{trim_redundant_ALTs} = 1 if ($self->trim_ALTs());
    $$opts{silent_dups} = 1 if ($self->silent());
    $$opts{vcf_header} = 1 if ($self->vcf_header());
    $$opts{regions_list} = 1 if ($self->regions());
    if ($self->collapse()) {
    	my $c = $self->collapse();
		$$opts{collapse_any} = 0;
		if ( $c eq 'snps' ) { $$opts{collapse_snps} = 1; }
		elsif ( $c eq 'indels' ) { $$opts{collapse_indels} = 1; }
		elsif ( $c eq 'both' ) {
			$$opts{collapse_snps} = 1;
			$$opts{collapse_indels} = 1;
		}
		elsif ( $c eq 'any' ) { $$opts{collapse_any} = 1; }
		elsif ( $c eq 'none' ) {
			$$opts{collapse_any} = 0;
			$$opts{collapse_snps} = 0;
			$$opts{collapse_indels} = 0;
		}
		else { error("Expected one of <snps|indels|both|any> with -c"); }
    }
    return $opts;
}


#sub parse_params
#{
#    $0 =~ s{^.+/}{}; $0 .= "($Vcf::VERSION)";
#    my $opts = 
#    { 
#        args => [$0, @ARGV],
#        joiner => 
#        {
#            'I16' => \&joiner_dp4,
#            'DP4' => \&joiner_dp4,
#            'MQ0'  => \&joiner_sum,
#            'DP'  => \&joiner_sum,
#        },
#        trim_redundant_ALTs => 0,
#        collapse_any => 1,
#    };
#    while (my $arg=shift(@ARGV))
#    {
#        if ( $arg eq '-d' || $arg eq '--remove-duplicates' ) { $$opts{rm_dups}=1; next; }
#        if ( $arg eq '-R' || $arg eq '--ref-for-missing' ) { $$opts{ref_for_missing}=shift(@ARGV); next; }
#        if ( $arg eq '-t' || $arg eq '--trim-ALTs' ) { $$opts{trim_redundant_ALTs}=1; next; }
#        if ( $arg eq '-s' || $arg eq '--silent' ) { $$opts{silent_dups}=1; next; }
#        if ( $arg eq '-H' || $arg eq '--vcf-header' ) { $$opts{vcf_header}=shift(@ARGV); next; }
#        if ( $arg eq '-T' || $arg eq '--tabix_path' ) { $$opts{tabix_path}=shift(@ARGV); next; }
#        if ( $arg eq '-c' || $arg eq '--collapse' ) 
#        {
#            $$opts{collapse_any} = 0;
#            my $c = shift(@ARGV);
#            if ( $c eq 'snps' ) { $$opts{collapse_snps}=1; }
#            elsif ( $c eq 'indels' ) { $$opts{collapse_indels}=1; }
#            elsif ( $c eq 'both' ) { $$opts{collapse_snps}=1; $$opts{collapse_indels}=1; }
#            elsif ( $c eq 'any' ) { $$opts{collapse_any}=1; }
#            elsif ( $c eq 'none' ) { $$opts{collapse_any}=0; $$opts{collapse_snps}=0; $$opts{collapse_indels}=0; }
#            else { error("Expected one of <snps|indels|both|any> with -c"); }
#            next; 
#        }
#        if ( $arg eq '-r' || $arg eq '--regions' ) { $$opts{regions_list}=shift(@ARGV); next; }
#        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
#        if ( -e $arg ) { push @{$$opts{files}},$arg; next; }
#        error("Unknown parameter or non-existent file \"$arg\". Run -? for help.\n");
#    }
#    if ( !exists($$opts{files}) ) { error() }
#    return $opts;
#}

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { croak join('',@msg); }
    die
        "About: Merges VCF files by position, creating multi-sample VCFs from fewer-sample VCFs.\n",
        "   The tool requires bgzipped and tabix indexed VCF files on input. (E.g. bgzip file.vcf; tabix -p vcf file.vcf.gz)\n",
        "   If you need to concatenate VCFs (e.g. files split by chromosome), look at vcf-concat instead.\n",
        "Usage: vcf-merge [OPTIONS] file1.vcf file2.vcf.gz ... > out.vcf\n",
        "Options:\n",
        "   -c, --collapse <snps|indels|both|any|none>      treat as identical sites with differing alleles [any]\n",
        "   -d, --remove-duplicates                         If there should be two consecutive rows with the same chr:pos, print only the first one.\n",
        "   -H, --vcf-header <file>                         Use the provided VCF header\n",
        "   -h, -?, --help                                  This help message.\n",
        "   -r, --regions <list|file>                       Do only the given regions (comma-separated list or one region per line in a file).\n",
        "   -R, --ref-for-missing <string>                  Use the REF allele instead of the default missing genotype. Because it is not obvious\n",
        "                                                       what ploidy should be used, a user-defined string is used instead (e.g. 0/0).\n",
        "   -s, --silent                                    Try to be a bit more silent, no warnings about duplicate lines.\n",
        "   -t, --trim-ALTs                                 If set, redundant ALTs will be removed\n",
        "\n";
}


# Returns the common prefix of the files.
sub common_prefix
{
    my ($files) = @_;
    my @paths;
    my $len = -1;
    for my $file (@$files)
    {
        my @path = split(m{/+},$file);
        if ( $len<0 || $len>scalar @path ) { $len=scalar @path; }
        push @paths, \@path;
    }
    my @common;
    for (my $i=0; $i<$len; $i++)
    {
        my $identical=1;
        for (my $ifile=1; $ifile<scalar @paths; $ifile++)
        {
            if ( $paths[$ifile]->[$i] ne $paths[0]->[$i] ) { $identical=0; last; }
        }
        if ( !$identical ) { last; }
        push @common, $paths[0]->[$i];
    }
    return join('/+',@common);
}


sub read_region_list
{
    my ($opts) = @_;

    my @regions = ();
    if ( exists($$opts{regions_list}) ) 
    { 
        if ( -e $$opts{regions_list} )
        {
            open(my $rgs,'<',$$opts{regions_list}) or error("$$opts{regions_list}: $!");
            while (my $line=<$rgs>)
            {
                chomp($line);
                push @regions, $line;
            }
            close($rgs);
        }
        else
        {
            @regions = split(/,/,$$opts{regions_list}); 
        }
    }
    return (@regions);
}

sub check_AGtags_definition
{
    my ($vcf) = @_;
    if ( $$vcf{version} >= 4.1 ) { return; }

    # Whatever is the value set to, the user takes the responsibility for the merging strategy used
    if ( exists($ENV{DONT_FIX_VCF40_AG_TAGS}) ) { return; }

    my @tags;
    if ( exists($$vcf{header}{INFO}{PL}) && $$vcf{header}{INFO}{PL}{Number} != -1  ) { push @tags, 'PL'; }
    if ( exists($$vcf{header}{INFO}{GL}) && $$vcf{header}{INFO}{GL}{Number} != -1  ) { push @tags, 'GL'; }
    if ( exists($$vcf{header}{INFO}{AC}) && $$vcf{header}{INFO}{AC}{Number} != -1  ) { push @tags, 'AC'; }
    if ( exists($$vcf{header}{INFO}{AF}) && $$vcf{header}{INFO}{AF}{Number} != -1  ) { push @tags, 'AF'; }

    if ( !@tags ) { return; }

    $ENV{DONT_FIX_VCF40_AG_TAGS} = 1;
    my $tags = join(',',@tags);
    print STDERR 
        "Warning: The $tags tag(s) will not be merged correctly for multiallelic sites.\n",
        "   To be handled correctly, please redefine with Number=. or set the environment\n",
        "   variable DONT_FIX_VCF40_AG_TAGS=0.\n";
}

sub init_cols
{
    my ($opts,$vcf_out) = @_;

    my $prefix;
    my @regions = read_region_list($opts);
    my @vcfs;
    my @cols;
    my %has_chrom;
    my %col_names;
    my $icol = 9;
    my $ncols_total = 0;

    if ( !$$opts{has_col_names} ) { $prefix = common_prefix($$opts{files}); }

    # Go through all files and read header, obtain list of chromosomes. The file names will be used for columns, unless
    #   they were read from the header.
    for my $file (@{$$opts{files}})
    {
        my $vcf = Vcf->new( tabix=>$$opts{tabix_path}, file=>$file);
        $$vcf{line_buffer} = [];
        $vcf->parse_header();
        check_AGtags_definition($vcf);
        $vcf->close();
        push @vcfs, $vcf;

        # Precompute the weighting factor for the QUAL column
        my $ncols = scalar @{$$vcf{columns}} - 9;
        if ( $ncols<=0 ) { $ncols = 1; }
        $$vcf{qual_weight} = 1.0*$ncols;
        $ncols_total += $ncols;

        # Update the list of known chromosomes
        if ( !exists($$opts{regions_list}) )
        {
            my $chrms = $vcf->get_chromosomes();
            for my $chr (@$chrms)
            {
                if ( exists($has_chrom{$chr}) ) { next; }
                $has_chrom{$chr} = 1;
                push @regions, $chr;
            }
        }

        my $col_prefix = '';
        if ( !$$opts{has_col_names} )
        {
            # Make the column names nice - strip common prefix and the suffix .vcf.gz
            $col_prefix = $file; 
            $col_prefix =~ s{^/*$prefix/*}{};
            $col_prefix =~ s/\.gz$//i;
            $col_prefix =~ s/\.vcf$//i;
            $col_prefix .= '_';
        }

        if ( !exists($$vcf{columns}) ) { error("No header present? $file\n"); }

        # Create good names for the columns in the merged vcf file
        my @vcf_cols = @{$$vcf{columns}};
        $$vcf{__col_names} = [];
        for my $col (@vcf_cols[9..$#vcf_cols])
        {
            my $col_name = $col;
            if ( $$opts{has_col_names} ) 
            {
                if ( $icol >= @{$$vcf_out{columns}} ) { error("Fewer columns in the header than in the VCF files total.\n"); }
                $col_name = $$vcf_out{columns}[$icol];
                $icol++;

                if ( exists($col_names{$col_name}) ) { error("The column names not unique in the header: $col_name\n"); }
            }
            else
            {
                if ( exists($col_names{$col_name}) ) { $col_name = $col_prefix.$col; }
                if ( exists($col_names{$col_name}) ) { warn("FIXME: the column name [$col_name] not unique.\n"); }
            }
            warn("Using column name '$col_name' for $file:$col\n");
            $col_names{$col_name} = 1;

            push @cols, $col_name;
            push @{$$vcf{__col_names}}, $col_name;
        }
    }

    if ( $$opts{has_col_names} && $icol!=@{$$vcf_out{columns}} ) { error("More columns in the header than in the VCF files total.\n"); }

    # QUAL weighting
    for my $vcf (@vcfs)
    {
        $$vcf{qual_weight} /= $ncols_total;
    }

    $$opts{vcfs} = \@vcfs;
    $$opts{cols} = \@cols;
    $$opts{regions} = \@regions;
}


sub merge_vcf_files
{
	my $self = shift;
    my $opts = $self->parse_params();
    
    open (OUT, ">".$$opts{outfile});
    
    # Create output VCF
    my $vcf_out;
    if ( $$opts{vcf_header} )
    {
        $vcf_out = Vcf->new( tabix=>$$opts{tabix_path}, file=>$$opts{vcf_header});
        $vcf_out->parse_header();
        if ( $$vcf_out{columns} && @{$$vcf_out{columns}} ) { $$opts{has_col_names}=1; }
    }
    else
    {
        $vcf_out = Vcf->new( tabix=>$$opts{tabix_path} );
    }

    $$vcf_out{trim_redundant_ALTs} = $$opts{trim_redundant_ALTs};
    
    init_cols($opts,$vcf_out);
    my @regions = @{$$opts{regions}};
    my @cols    = @{$$opts{cols}};
    my @vcfs    = @{$$opts{vcfs}};

    #warn Dumper $opts; die;
    
    # Get the header of the output VCF ready
    $vcf_out->add_columns(@cols);
    if ( !$$vcf_out{has_header} )
    {
        for my $vcf (@vcfs)
        {
            # To get the missig fields filled by the default values
            for my $hline (@{$$vcf{header_lines}})
            {
                if ( $$hline{key} eq 'fileformat' ) { next; }
                $vcf_out->add_header_line($hline,silent=>1);
            }
        }
    }

    # List source files
    my $source;
    for (my $i=0; $i<@vcfs; $i++)
    {
        if ( $i ) { $source .= ','; }
        $source .= "$i:$vcfs[$i]{file}";
    }
    $vcf_out->add_header_line({key=>'source',value=>join(' ',@{$$opts{args}})},append=>'timestamp');
    $vcf_out->add_header_line({key=>'sourceFiles',value=>$source},append=>'timestamp');
    $vcf_out->add_header_line({key=>'INFO',ID=>'SF',Number=>-1,Type=>'String',Description=>'Source File (index to sourceFiles, f when filtered)'});

    my $have_samples = @{$$vcf_out{columns}}>9 ? 1 : 0; 

    $vcf_out->recalc_ac_an($have_samples ? 2 : 0);
    $vcf_out->add_header_line({key=>'INFO',ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'});
    $vcf_out->add_header_line({key=>'INFO',ID=>'AN',Number=>1,Type=>'Integer',Description=>'Total number of alleles in called genotypes'});
    print OUT $vcf_out->format_header();

    # Go through all VCF files simultaneously and output each line, one region at a time.
    for my $region (@regions)
    {
        # Open files
        for my $vcf (@vcfs) 
        { 
            delete($$vcf{done});
            $vcf->open(region=>$region); 
        }

        while ( my $pos=advance_position($opts,\@vcfs) )
        {
            my %out;
            $out{POS}    = $pos;
            $out{ID}     = '.';
            $out{ALT}    = [];
            $out{FORMAT} = [];
            my %format;
            my %info;
            my @src_files;
            my %filters;
            my (@quals,@qual_weights,$qual_weights_sum,%ac,$an);

            my %ref_alt_map = ();
            # Find out the REFs and ALTs: in VCFv4.0, the REFs can differ and ALTs must be converted
            for my $vcf (@vcfs)
            {
                my $line = $$vcf{last_line};
                if ( !$line ) { next; }
                if ( !exists($out{CHROM}) ) { $out{CHROM} = $$line{CHROM}; }
                my $ref = $$line{REF};
                for my $alt (@{$$line{ALT}}) { $ref_alt_map{$ref}{$alt}=$alt; }
            }
            # Do the REF,ALT conversion only when necessary
            my $new_ref; 
            if ( scalar keys %ref_alt_map > 1 ) 
            { 
                $new_ref = $vcf_out->fill_ref_alt_mapping(\%ref_alt_map); 
                if ( !defined $new_ref ) { error("Failed on line $out{CHROM}:$out{POS}\n"); }
            }
            if ( !$have_samples or !$$opts{trim_redundant_ALTs} )
            {
                # Do not loose information from the ALT column when samples are not present
                my %alts;
                for my $vcf (@vcfs)
                {
                    my $line = $$vcf{last_line};
                    if ( !$line ) { next; }
                    my $ref = $$line{REF};
                    for my $alt (@{$$line{ALT}}) { $alts{$ref_alt_map{$ref}{$alt}}=1; }
                    delete($alts{'.'});
                    $out{ALT} = [ keys %alts ];
                }
            }
            for (my $ivcf=0; $ivcf<@vcfs; $ivcf++)
            {
                my $vcf  = $vcfs[$ivcf];
                my $line = $$vcf{last_line};

                # If this file does not have a record for this position, then for all its columns output undef gtype
                if ( !$line )
                {
                    for (my $i=0; $i<@{$$vcf{__col_names}}; $i++)
                    {
                        my $name = $$vcf{__col_names}->[$i];
                        $out{gtypes}{$name}{GT} = exists($$opts{ref_for_missing}) ? $$opts{ref_for_missing} : $$vcf_out{defaults}{GT};
                    }
                    next;
                }

                # Check if the site has been filtered
                if ( scalar @{$$line{FILTER}}>1 or ($$line{FILTER}[0] ne $$vcf{filter_passed} && $$line{FILTER}[0] ne $$vcf{defaults}{default}) )
                {
                    push @src_files,$ivcf.'f';
                }
                else
                {
                    push @src_files,$ivcf;
                }

                # Collect information for the FILTER field
                for my $flt (@{$$line{FILTER}})
                {
                    if ( $flt eq $$vcf{filter_passed} )
                    {
                        $filters{$$vcf_out{filter_passed}} = 1;
                    }
                    elsif ( $flt ne $$vcf{defaults}{default} )
                    {
                        $filters{$flt} = 1;
                    }
                }

                # Collect information for the QUAL field
                if ( $$line{QUAL} ne $$vcf{defaults}{QUAL} && $$line{QUAL} ne $$vcf{defaults}{default} && $$line{QUAL}>0 )
                {
                    push @quals,$$line{QUAL};
                    push @qual_weights,$$vcf{qual_weight};
                    $qual_weights_sum += $$vcf{qual_weight};
                }

                if ( $$line{ID} ne '.' && $out{ID} eq '.' ) { $out{ID}=$$line{ID}; }

                # Remember the FORMAT fields
                for my $field (@{$$line{FORMAT}}) { $format{$field} = 1; }

                # VCF without genotypes: calculate AC,AN if present
                if ( !$have_samples )
                {
                    if ( exists($$line{INFO}{AN}) ) { $an += $$line{INFO}{AN}; }
                    if ( exists($$line{INFO}{AC}) ) 
                    { 
                        my (@acs) = split(/,/,$$line{INFO}{AC});
                        for (my $i=0; $i<@acs; $i++) 
                        { 
                            my $alt = $ref_alt_map{$$line{REF}}{$$line{ALT}[$i]};
                            $ac{$alt} += $acs[$i];
                        }
                    }
                }

                # Join the INFO field
                for my $inf (keys %{$$line{INFO}}) 
                {
                    # When conflicting INFO fields are present, use the first one, unless a joining method exists
                    if ( exists($info{$inf}) )
                    {
                        if ( exists($$opts{joiner}{$inf}) )
                        {
                            &{$$opts{joiner}{$inf}}(\$info{$inf},$$line{INFO}{$inf});
                        }
                        next;
                    }
                    $info{$inf} = $$line{INFO}{$inf};
                }

                my $ref = $$line{REF};

                # The ALT column may change after the merge, take care of ALT dependent tags such as GL.
                if ( $have_samples )
                {
                    if ( defined $new_ref )
                    {
                        $vcf->parse_AGtags($line,\%ref_alt_map,$$line{REF});
                    }
                    else
                    {
                        $vcf->parse_AGtags($line);
                    }
                }

                # Now fill in the genotype information for each column
                for (my $i=0; $i<@{$$vcf{__col_names}}; $i++)
                {
                    my $ori_name = $$vcf{columns}->[$i+9];
                    my $out_name = $$vcf{__col_names}->[$i];

                    $out{gtypes}{$out_name} = $$line{gtypes}{$ori_name};

                    # This is to convert 0/1 to G/C
                    my ($alleles,$seps,$is_phased,$is_empty) = $vcf->parse_haplotype($line,$ori_name);
                    if ( defined $new_ref )
                    {
                        my @als;
                        for my $al (@$alleles) 
                        { 
                            push @als, exists($ref_alt_map{$ref}{$al}) ? $ref_alt_map{$ref}{$al} : '.';
                        }
                        $out{gtypes}{$out_name}{GT} = $vcf->format_haplotype(\@als,$seps);
                    }
                    else
                    {
                        $out{gtypes}{$out_name}{GT} = $vcf->format_haplotype($alleles,$seps);
                    }
                }
                $out{REF} = defined $new_ref ? $new_ref : $ref;
            }

            $out{INFO} = { %info }; 
            $out{INFO}{SF} = join(',',@src_files);

            # Output the QUAL information
            my $qual;
            for (my $i=0; $i<@quals; $i++)
            {
                $qual += $quals[$i] * $qual_weights[$i] * (1.0 / $qual_weights_sum);
            }
            $out{QUAL} = defined $qual ? sprintf("%.2f",$qual) : $$vcf_out{defaults}{QUAL};

            # Output the FILTER information: remove PASS or missing value if some other information
            #   is present.
            delete($filters{$$vcf_out{defaults}{default}});
            if ( exists($filters{$$vcf_out{filter_passed}}) && scalar keys %filters > 1 )
            {
                delete($filters{$$vcf_out{filter_passed}});
            }
            $out{FILTER} = [ keys %filters ];
            if ( !@{$out{FILTER}} ) { push @{$out{FILTER}},$$vcf_out{defaults}{default}; }

            # The GT field must come as first
            delete($format{GT});
            $out{FORMAT} = ['GT'];
            for my $key (keys %format) { push @{$out{FORMAT}},$key; }

            if ( $have_samples )
            {
                $vcf_out->format_genotype_strings(\%out);
            }
            else
            {
                if ( defined $an ) { $out{INFO}{AN}=$an; }
                if ( scalar keys %ac ) 
                { 
                    my @acs;
                    for my $alt (@{$out{ALT}})
                    {
                        # Some of the files may not have AC, the AC count can be undefined in such a case.
                        push @acs, exists($ac{$alt}) ? $ac{$alt} : 0; 
                    }
                    $out{INFO}{AC} = join(',',@acs);
                }
            }
            print OUT $vcf_out->format_line(\%out);
        }
    }
    for my $vcf (@vcfs) 
    { 
        $vcf->close() or error("close failed: $$vcf{file}\n"); 
    }
    close (OUT);
}


sub advance_position
{
    my ($opts,$vcfs) = @_;

    my $min_pos;
    for my $vcf (@$vcfs)
    {
        fill_buffer($opts,$vcf) unless $$vcf{done};
        if ( @{$$vcf{line_buffer}} && (!defined $min_pos or $min_pos>$$vcf{line_buffer}[0]{POS}) ) { $min_pos = $$vcf{line_buffer}[0]{POS}; }
    }
    if ( !defined $min_pos ) { return undef; }

    my ($first,$has_snp,$has_indel);
    for my $vcf (@$vcfs)
    {
        delete($$vcf{last_line});

        if ( @{$$vcf{line_buffer}} && $min_pos ne $$vcf{line_buffer}[0]{POS} ) { next; }
        if ( !defined $first ) 
        { 
            $$vcf{last_line} = shift @{$$vcf{line_buffer}};
            $first = $$vcf{last_line};
            next;
        }
        my $irec;
        for (my $i=0; $i<@{$$vcf{line_buffer}}; $i++)
        {
            my $line = $$vcf{line_buffer}[$i];
            if ( $$line{POS} ne $$first{POS} ) { last; }

            if ( $$opts{collapse_any} ) { $irec=$i; last; }  # checking position only
            if ( $$opts{collapse_snps} && $$first{variant_type}&1 && $$line{variant_type}&1 ) { $irec=$i; last; }
            if ( $$opts{collapse_indels} && $$first{variant_type}&2 && $$line{variant_type}&2 ) { $irec=$i; last; }

            if ( $$vcf{line_buffer}[$i]{REF} ne $$first{REF} ) { next; } # refs do not match
            for my $al1 (@{$$line{ALT}})
            {
                for my $al2 (@{$$first{ALT}})
                {
                    if ( $al1 eq $al2 ) { $irec=$i; last; }
                }
                if ( defined $irec ) { last; }
            }
            if ( defined $irec ) { last; }
        }
        if ( defined $irec )
        {
            $$vcf{last_line} = splice(@{$$vcf{line_buffer}},$irec,1);
        }
    }
    return $min_pos;
}

sub fill_buffer
{
    my ($opts,$vcf) = @_;
    if ( @{$$vcf{line_buffer}} && $$vcf{line_buffer}[0]{POS}!=$$vcf{line_buffer}[-1]{POS} ) { return; }
    while ( 1 )
    {
        my $line = $vcf->next_data_hash();
        if ( !$line ) { $$vcf{done} = 1; return; }
        if ( !$$opts{collapse_any} ) 
        {
            for my $al (@{$$line{ALT}})
            {
                my ($type,$len,$ht) = $vcf->event_type($$line{REF},$al);
                if ( $type eq 's' or $type eq 'r' ) { $$line{variant_type} |= 1; }
                if ( $type eq 'i' or $type eq 'o' ) { $$line{variant_type} |= 2; }
            }
        }
        push @{$$vcf{line_buffer}}, $line;
        if ( $$vcf{line_buffer}[0]{POS} != $$vcf{line_buffer}[-1]{POS} ) { return; }
    }
}


# Field joiner methods
sub joiner_sum
{
    my ($ori,$new) = @_;
    $$ori += $new;
}

sub joiner_dp4
{
    my ($ori,$new) = @_;
    my @vals1 = split(/,/,$$ori);
    my @vals2 = split(/,/,$new);
    if ( @vals1 != @vals2 ) { error("Cannot join: $$ori vs $new\n"); }
    for (my $i=0; $i<@vals1; $i++)
    {
        $vals1[$i] += $vals2[$i];
    }
    $$ori = join(',',@vals1);
}


1;