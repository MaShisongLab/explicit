use strict;
my $color_file = "./data/colors.txt";
my $regulator_file = "results.regulator.tfs.txt";
my $out_file = "chord.lists.txt";
my $qvalue_cutoff = 1E-05;
my $tf_num = 50;
my $tg_num = 15;

unless ( -e $color_file ){
	print "The file 'colors.txt' was not found with the data folder. EXIT in 3 second.\n";
	sleep 3;
	exit 1;
}

unless ( -e $regulator_file ){
	print "The file 'results.regulator.tfs.txt' was not found. EXIT in 3 second.\n";
	sleep 3;
	exit 1;
}

our $lower_bound = (- log($qvalue_cutoff) / log(10));

my $l;
my @t;
my @rand_colors;
my @colors;

open IN, $color_file;
while ( $l = <IN>){
	chomp $l;
	@t = split("\t", $l);
	if ( $t[1] eq "random_color"){
		push @rand_colors, $t[0];
	}elsif( $t[1] eq "darkblue2darkred"){
		push @colors, $t[0];
	}
}
close IN;

my $m = $ARGV[0];
chomp $m;

if ( length( $m ) < 1 ){
	print "\nNo module name given in the command line.\nPlease input module name and press enter:";
	$m = <STDIN>;
	chomp $m;
}

if ( $ARGV[1] >= 1){
	$tf_num = $ARGV[1];
}

if ( $ARGV[2] >= 1){
	$tg_num = $ARGV[2];
}

my @tf;	
my %pv;		# pvalue for tf enrichment
my @tg;		# target gene
my @b;		# regression coefficient beta
my @b_pv;	# pvalue for beta
my $i;
my %count;
my %edge;
my %sum;
my %beta_v;
my @gene;

open IN, $regulator_file;
$l = <IN>;
while ( $l = <IN>){
	chomp $l;
	$l =~ s/ //g;
	$l =~ s/"//g;
	@t = split("\t", $l);
	if ( $t[0] eq $m & $t[1] <= $tf_num){
		push @tf, $t[3];
		$pv{$t[3]} = $t[9];
		@tg = split("\/", $t[13]);
		@b = split("\/", $t[14]);
		@b_pv = split("\/",$t[15]);
		my $n = scalar @tg;
		for( $i = 0; $i < $n; $i++){
			$count{$tg[$i]} ++;
			$edge{$tg[$i]} .= "$t[3] $tg[$i] $b[$i] $b_pv[$i]\t";
			$sum{$tg[$i]} += abs($b[$i]);
			$beta_v{"$t[3] $tg[$i]"} = $b[$i];
		}
	}
}
close IN;

my $k;

foreach $k (keys %count){
	push @gene, {'gene' => $k, 'count' => $count{$k}, 'sum' => $sum{$k} };
}

my @sorted_gene = sort {
	$b->{'count'} <=> $a->{'count'} ||
	$b->{'sum'} <=> $a->{'sum'}
} @gene;

if ( scalar(@gene) < 1 ){
	print "\n\nModule $m not found in the file$regulator_file.\nEXIT in 3 seconds.\n";
	sleep 3;
	exit 1;
}

if ( scalar @sorted_gene < $tg_num){
	$tg_num = scalar @sorted_gene;
}

my @q;
my %seen;

for ( $i = 0; $i < $tg_num; $i++ ){
	$k = $sorted_gene[$i];
	$k = $k->{'gene'};
	$l = $edge{$k};
	@t = split("\t", $l);
	foreach( @t ){
		if ( length( $_ ) > 3 ){
			@q = split(" ", $_);
			$seen{$q[0]} = 1;
		}
	}
}

my $j = 0;
my $pvalue;
my %node_col;
my %id;

foreach (@tf){
	$id{$_} = $j;
	$pvalue = $pv{$_};
	my $col = getColor(- log($pvalue) / log(10));
	$node_col{$_} = $col;
	$j++;
}

open OUT, ">$out_file";

for ( $i = 0; $i < $tg_num; $i++ ){
	$k = $sorted_gene[$i];
	$k = $k->{'gene'};
	$l = $edge{$k};
	@t = split("\t", $l);
	foreach( @t ){
		if ( length( $_ ) > 3 ){
			my $p = $_;
			$p =~ s/ /\t/g;
			my @q = split("\t", $p);
			my $beta_pv = $q[3];
			my $beta_color = getColor ( $beta_pv );
			my $ppvv = abs($q[2]);
			print OUT "$q[0]\t$q[1]\t$ppvv\t$node_col{$q[0]}\t$rand_colors[$id{$q[0]}]\t$beta_color\n"
		}
	}
}

close OUT;


########################

sub getColor
{
	my $cp = $_[0];
	my $t;

	my $n = 11;
	if ( $cp eq "Inf" ){
		$t = 200 - $n;
	}else{
		$t =  $cp;
		if ( $t > (200 - $n)){
			$t = 200 - $n;
		}
	}

	my $q = $t / ($t + $n);
	my $q2 =  $lower_bound / ( $lower_bound + $n );
	my $q3 = ( 200 - $n) / 200;
	my $q4 = ( $q3 - $q2 ) / 150;

	$q = int( ($q - $q2) / $q4);

	if ( $q < 0 ) {$q = 0;}
	if ( $q > 150) {$q = 150;}

	return ( $colors[$q] );
}


