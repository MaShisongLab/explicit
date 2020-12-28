use strict;
my $pvalue_cutoff = 1E-9;  			# pValue cut-off for TF-target gene pair.  
my $qvalue_cutoff = 1E-5;  			# qvalue cut-off for a TF to be considered as predictor TF for a module.
my $total_gene = 29182; 			# the total number of non-TF target genes.
my $cluster_file = "modules_to_analyze.txt"; 	# The file with cluster information. The first column is gene id, the second column is gene name.
my $target_gene_file = "./data/At.target.genes.txt";
my $gene_symbl_file = "./data/At.gene.symbls.txt";
my $edge_file = "./data/At.SigEdges.txt"; 		# with singinificant interacting TF-target gene pair.
my $output_file = "results.regulator.tfs.txt";

sub LnGamma;
sub pValue; 	#( int in_c, int in_g, int c_size, int total)
sub HyperG;	# hypergeomteric distribution
sub BH;		# adjustment for multiple testing vi BH procedure.

my %e;		# to store predictor TFs info for target genes.
my %tfc;	# number of target genes for a TF within the genome.
my %seen;
my $l;
my @t;
my %symbl;
my @c;
my $gene;
my $n;
my $ci;
my $o;

unless( -e $edge_file ){
	print "The file 'At.SigEdges.txt' not found. It should be saved in the 'data' sub-directory. And the location of the script 'getArabidopsisRegulatorTFs.pl' should not be changed.\nExit in 5 second.\n";
	sleep 5;
	exit 1;
}

unless( -e $gene_symbl_file ){
	print "The file 'At.gene.symbls.txt' not found. It should be saved in the 'data' sub-directory. And the location of the script 'getArabidopsisRegulatorTFs.pl' should not be changed.\nExit in 5 seconds.\n";
	sleep 5;
	exit 1;
}

unless( -e $target_gene_file ){
	print "Target gene file 'At.target.genes.txt' not found. It should be saved in the 'data' sub-directory. And the location of the script 'getArabidopsisRegulatorTFs.pl' should not be changed.\nExit in 5 seconds.\n";
	sleep 5;
	exit 1;
}

unless( -e $cluster_file ){
	print "The file 'modules_to_analyze.txt' not found. It should be saved in the same directory as the script 'getArabidopsisRegulatorTFs.pl'.\nExit in 5 seconds.\n";
	sleep 5;
	exit 1;
}


our @Logfact;
for ( my $i = 0; $i < 1000000; $i++) { $Logfact[$i] = LnGamma ($i + 1) ; } # Log factorial values for 0 to 1000000


open IN, $edge_file ;
while ( $l = <IN>){
	chomp $l;
	@t = split("\t",$l);
	my $pc = int($t[4] * 10000 + 0.5) / 10000;  #deleted
	my $coef = int($t[2] * 10000 + 0.5) / 10000;
	if ( $t[3] <= $pvalue_cutoff ){
		$e{$t[0]} .= "$t[1] $coef $t[3] $pc\t";
		$tfc{$t[1]} ++;
	}
}

close IN;

open IN, $target_gene_file;
while ( $l = <IN>){
	chomp $l;
	$seen{$l} = 1;
}
close IN;

open IN, $gene_symbl_file;
while ( $l = <IN>){
	chomp $l;
	@t = split("\t", $l);
	$symbl{$t[0]} = $t[1];
}
close IN;

my %g;		# store the genes for clusters.
my %count;	# gene count for clusters.
my @t;

open IN, $cluster_file;
while ( $l = <IN>){
	chomp $l;
	@t = split("\t", $l);
	if ( $seen{$t[0]}){				# only consider target genes
		unless ( $seen{$l} ){			# remove duplicated gene-cluser line
			$g{$t[1]} .= "$t[0]\t";
			$count{$t[1]} ++;
			unless( $seen{$t[1]}){
				push( @c, $t[1]);
				$seen{$t[1]} = 1;
			}
			$seen{$l} = 1;
		}
	}
}
close IN;

unless( scalar( keys %count) >= 1){
	print "No clusters found. Please double-check the module information, especially the gene name.\nExit in 5 seconds.\n";
	sleep 5;
	exit 1;
}

my %ept;
my @ept;
my @tf;
my $r;
my $symbl_now;
my $k;

open OUT, ">$output_file";
print OUT "Module\tRank\tTF\tSymbl\tModuleSize(without TF genes)\tCount\tCountInGenome\tGenomeSize\tpValueEnrichment\tpValue(bh-adjusted)\tFraction\tmean_beta\tTargetAGI\tTargetSymbl\tbeta\tTargetpValue( -log10 p )\n";
foreach $k ( @c ){
	$l = $g{$k};
	$l =~ s/\t$//;
	@t = split("\t",$l);
	my %beta = %ept;
	my %num = %ept;
	my %pvs = %ept;
	my %b_s = %ept;
	my %syb = %ept;
	my %idg = %ept;
	my @line = @ept;
	my @idx = @ept;
	my @table = @ept;
	foreach $gene (@t){
		@tf = split("\t",$e{$gene});

		if (length( $symbl{$gene} )> 2 ){
			$symbl_now = $symbl{$gene};
		}else{
			$symbl_now = $gene;
		}

		foreach $r (@tf){
			my @m = split(" ", $r);
			$beta{$m[0]} += $m[1];
			$num{$m[0]} ++;
			$pvs{$m[0]} .= "$m[2] ";
			$syb{$m[0]} .= "$symbl_now#####";
			$idg{$m[0]} .= "$gene ";
			my $b_now = int( $m[1] * 1000 + 0.5) / 1000;
			$b_s{$m[0]} .= "$b_now ";
		}
	}
	foreach $n (keys %beta){
		my $n_c = $num{$n};
		my $n_b = sprintf "%0.4f", $beta{$n} / $num{$n};	# average beta
		
		my $ps = $pvs{$n};			
		$ps =~ s/" $"//;
		my @pp = split(" ",$ps);		# array of pvalues

		$b_s{$n} =~ s/" $"//;		
		my @c_b_s = split(" ", $b_s{$n});	# array of beta

		$syb{$n} =~ s/"#####$"//;		# array of gene symbl
		my @c_syb = split("#####", $syb{$n});

		$idg{$n} =~ s/" $"//;			# array of gene
		my @c_id = split(" ", $idg{$n});

		my @si = sort { $pp[$a] <=> $pp[$b] } 0..$#pp;

		my @logp = @ept;
		my @top_gene_symbl = @ept;
		my @top_gene_id = @ept;
		my @top_gene_b = @ept;

		my $j = scalar( @pp );
		for ( my $i = 0; $i < $j; $i++){
			$ci = $si[$i];
			if ( $pp[$ci] <= 0){
				$logp[$i] = "Inf";
			}else{
				$logp[$i] = int( - log( $pp[$ci])/log (10) * 10 ) / 10;
			}
			$top_gene_symbl[$i] = $c_syb[$ci];
			$top_gene_b[$i] = $c_b_s[$ci];
			$top_gene_id[$i] = $c_id[$ci];

		}
		my $p5 = join("/", @logp);
		my $p6 = join("/", @top_gene_symbl);
		my $p7 = join("/", @top_gene_b);
		#$p8 = join("/", @top_gene_pc);
		my $p9 = join("/", @top_gene_id);

		my $percentage = int($n_c / $count{$k}*1000 + 0.5) / 1000 ;
		my $pvalue_enrichment = pValue($n_c, $tfc{$n}, $count{$k}, $total_gene); #( int in_c, int in_g, int c_size, int total)
		my $ll ="$n\t $symbl{$n}\t$count{$k}\t$n_c\t$tfc{$n}\t$total_gene\t$pvalue_enrichment\t0\t$percentage\t$n_b\t\"$p9\"\t\"$p6\"\t\"$p7\"\t\" $p5\"";
		if ( $n_c >= 1){
			push @table, {
				'line' => $ll,
				'count' => $n_c,
				'pvalue' => $pvalue_enrichment,
				'sum' => $beta{$n},
				'percentage' => $percentage
			};
		}
	}
	@table = sort {
	#	$b->{'sum'} <=> $a->{'sum'} ||
	#	$b->{'count'} <=> $a->{'count'} ||
		$a->{'pvalue'} <=> $b->{'pvalue'}
		} @table;
		
	my @raw_p = @ept;
	foreach $o ( @table) {
		push (@raw_p, $o->{'pvalue'});
	}
	my @bhpvalue = BH ( \@raw_p );

	my $i = 1;
	my $j = 0;
	foreach $o ( @table ){
		if ( $bhpvalue[$j] <= $qvalue_cutoff ){
			my $ll = $o->{'line'};
			@t = split("\t", $ll);
			$t[7] = $bhpvalue[$j];
			$ll = join("\t", @t);
			print OUT "$k\t$i\t$ll\n";
			$i++;
		}
		$j++;
	}
}
close OUT;

#######################################
sub LnGamma
{
	my $x = $_[0];

	if ($x < 12 ){
		my $multiple = 1;
		for ( my $i = 1; $i < $x; $i++){ $multiple *= $i; }
		return log($multiple);
	} else {
		my @c = (1.0/12.0,
			-1.0/360.0,
			1.0/1260.0,
			-1.0/1680.0,
			1.0/1188.0,
			-691.0/360360.0,
			1.0/156.0,
			-3617.0/122400.0,
			43867.0/244188.0,
			-174611.0/125400.0
		);
		my $x_square = 1.0 / ( $x * $x );
		my $factor = 1.0;
		my $sum = $c[0];
		for (my $i = 1; $i < 10; $i++)
		{
			$factor *= $x_square;
			$sum += $c[$i] * $factor;
		}
		$sum = $sum / $x;
		my $Pi = 3.141592653589793238462643383279502884197169;
		my $lnGamma = ( $x - 0.5)*log( $x ) - $x + 0.5 * log( 2 * $Pi ) + $sum ;
		return $lnGamma;
	}
}



######################################
sub pValue #( int in_c, int in_g, int c_size, int total)
{
	my $in_c = $_[0];	# number of promoters in the cluster containing the motif
	my $in_g = $_[1];	# number of promoters in the genome containing the motif
	my $c_size = $_[2];	# the size of the cluster
	my $total = $_[3];	# the total number of promoter is the genome
	
	my $p_1side = 0;
	
	my $max = $c_size;
	if ( $in_g < $c_size ){ $max = $in_g;}

	my $p_1side = HyperG( $in_c, $in_g, $c_size, $total );

	for ( my $i = $in_c + 1; $i <= $max; $i++){
		$p_1side += HyperG( $i, $in_g, $c_size, $total);
	}

	return $p_1side;
}

######################################
sub HyperG # ( int in_c, int in_g, int c_size, int total)
{
	my ($in_c, $in_g, $c_size, $total) = @_;
	my $p1 = $Logfact[$total - $in_g] + $Logfact[$in_g] + $Logfact[$c_size] + $Logfact[$total - $c_size];
	my $p2 = $Logfact[$in_c] + $Logfact[ $in_g - $in_c] + $Logfact[ $in_c + $total - $in_g - $c_size ] + $Logfact[ $c_size - $in_c] + $Logfact[ $total ];
	my $p3 = exp( $p1 - $p2 );
	return $p3;
}


######################################
sub BH # ( \@idx )
{
	my $p_array_ref = shift;
	my @p_value = @{$p_array_ref};
	my @sorted_pvalue = sort { $b <=> $a } @p_value;
	my @bh_pvalue = @sorted_pvalue;
	my $total_num = scalar ( @bh_pvalue);
	my $min_p = 1;
	my $i;
	if ( $total_num > 1 ){
		for ( $i = 0; $i < $total_num; $i++){
			$bh_pvalue[$i] = $sorted_pvalue[$i] * $total_num / ($total_num - $i) ;
			if ( $bh_pvalue[$i] > $min_p){
				$bh_pvalue[$i] = $min_p;
			}else{
				$min_p = $bh_pvalue[$i];
			}
		}
	}
	@bh_pvalue = reverse (@bh_pvalue);
	return @bh_pvalue;
}


