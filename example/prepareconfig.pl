#!/usr/bin/perl
#input files needed: trzval.txt and triangles.txt

#this code generates
#the necessary input files for 
#modified WRF nests on subsets
#of total numbr of processes 

use strict;

# Declare the subroutines
sub trim($);

#argument passed #processors

my $maxproc = $ARGV[0];
my $opt = $ARGV[1];

our (@decomp);
our ($cmd, $maxdom, $e_we, $e_sn, $iter, $xy);
our $nlinp = 'namelist.input';

our $outfile = 'nestconfig';

our (@left, @right);
our (@e_we_all, @e_sn_all);
our ($i, $j, $k, $l, $t, $size, $proc, $data, $char, $family, $ar, $iter);
our (@lines, @array, @sum, @triangles, @vertices, @nestconfig, @prediction_mid, @prediction_1024);
my ($x, $y, $x_1, $y_1, $x_2, $y_2, $x_3, $y_3, $lambda_1, $lambda_2, $lambda_3, $k_1, $k_2, $k_3);
my (@nest_sizes, @input_nest_sizes, @aspectratio, @input_aspectratio);
my ($id, $start_rank, $nest_x, $nest_y, $xgreatery, $leftset);
my ($shorter, $longer, @dim, $divisor, $slabs);

my %HoA = ();

my @prcr = qw(32 64 128 196 256 324 512 784 900 1024);

my @e_we = qw(415 487 259 403 283 322 277 310 178 184 187 133 94);
my @e_sn = qw(445 334 505 310 394 268 286 202 307 232 163 100 124);

my $mid = int($#prcr+1)/2;

foreach $i (0 .. $#e_we) {
	$size = $e_we[$i] * $e_sn[$i]; 
	$nest_sizes[$i] = $size;
	$ar = $e_we[$i] / $e_sn[$i];
	$aspectratio[$i] = $ar;
}

#get the nest config from namelist.input
&get_dom();

@e_we_all = split(/,/,$e_we);
@e_sn_all = split(/,/,$e_sn);
$e_we_all[0] = trim($e_we_all[0]);
$e_sn_all[0] = trim($e_sn_all[0]);

for $iter (1 .. $maxdom-1)
{
	$nestconfig[$iter-1][0] = $e_we_all[$iter];				#e_we
	$nestconfig[$iter-1][1] = $e_sn_all[$iter];				#e_sn
	$nestconfig[$iter-1][2] = $e_we_all[$iter] * $e_sn_all[$iter] ;		#xy
	$nestconfig[$iter-1][3] = $e_we_all[$iter] / $e_sn_all[$iter] ;		#x/y
}

open (FH, "< trzval.txt") || die "trzval.txt open failed: $!\n";#exec times for 13 training set data for procs 32 64 128 196 256 324 512 784 900 1024
@lines = <FH>;
close(FH);

#build table of numbers
for($i=0; $i<=$#lines; $i++)
{
	chomp($lines[$i]);
	$lines[$i] = trim($lines[$i]);
	@array = split(/ /, $lines[$i]);
	
	foreach $j (0 .. $#prcr)
	{
		push @{ $HoA{$i} }, $array[$j];			#push to the nest number i
	}

}

open (FH, "< triangles.txt") || die "triangles.txt open failed: $!\n";		#delaunay triangulation from qhull.py
@lines = <FH>;
close(FH);

for($i=1; $i<=$#lines; $i++)
{
	chomp($lines[$i]);
	$j = $i-1;
	$triangles[$j] = $lines[$i];
}

#    Point p lies inside the triangle if and only if 0 < \lambda_i < 1 \;\forall\; i \text{ in } 1,2,3.
#    Otherwise, p lies on the edge or corner of the triangle if 0 \leq \lambda_i \leq 1 \;\forall\; i \text{ in } 1,2,3.
#    Otherwise, p lies outside the triangle. 

foreach $i (1 .. $maxdom-1) {

 	$x = $nestconfig[$i-1][3]; $y = $nestconfig[$i-1][2]/100000.0;

	foreach $j (0 .. $#triangles) {

		@vertices = split(/ /, $triangles[$j]);
		$k_1 = $vertices[0];
		$x_1 = $aspectratio[$k_1]; $y_1 = $nest_sizes[$k_1]/100000.0;
		$k_2 = $vertices[1];
		$x_2 = $aspectratio[$k_2]; $y_2 = $nest_sizes[$k_2]/100000.0;
		$k_3 = $vertices[2];
		$x_3 = $aspectratio[$k_3]; $y_3 = $nest_sizes[$k_3]/100000.0;

		$lambda_1 = (($y_2-$y_3)*($x-$x_3) + ($x_3-$x_2)*($y-$y_3))/(($y_2-$y_3)*($x_1-$x_3) + ($x_3-$x_2)*($y_1-$y_3));
		$lambda_2 = (($y_3-$y_1)*($x-$x_3) + ($x_1-$x_3)*($y-$y_3))/(($y_2-$y_3)*($x_1-$x_3) + ($x_3-$x_2)*($y_1-$y_3));
		$lambda_3 = 1 - $lambda_1 - $lambda_2;

		if(($lambda_1 > 0 and $lambda_1 < 1) and ($lambda_2 > 0 and $lambda_2 < 1) and ($lambda_3 > 0 and $lambda_3 < 1) )
		{
			$l = $#prcr;
			$prediction_1024[$i-1][0] = $lambda_1*$HoA{$k_1}[$l] + $lambda_2*$HoA{$k_2}[$l] + $lambda_3*$HoA{$k_3}[$l];
			$l = $mid;
			$prediction_mid[$i-1][0] = $lambda_1*$HoA{$k_1}[$l] + $lambda_2*$HoA{$k_2}[$l] + $lambda_3*$HoA{$k_3}[$l];
		}
	
	}

}

$sum[0] = 0; $sum[1] = 0;
foreach $i (1 .. $maxdom-1) {
	$sum[0] = $sum[0] + $prediction_1024[$i-1][0];		#1024
	$sum[1] = $sum[1] + $prediction_mid[$i-1][0];		#324	#$prediction_mid[$i-1][0] --> execution time prediction on 324
}

&mpaspect($maxproc);

if($decomp[1] >= 2 * $decomp[0]) {
	$xgreatery = 1;
	$shorter = $decomp[0];
	$longer = $decomp[1];
	@dim = (4, 3);			#y, x
}
else {
	$xgreatery = 0;
	$shorter = $decomp[1];
	$longer = $decomp[0];
	@dim = (3, 4);			#x, y
}

#stage 1 to get share of processors
foreach $i (1 .. $maxdom-1) {

	$prediction_mid[$i-1][1] = $prediction_mid[$i-1][0]/$sum[1];			#ratio for nest $i
	$prediction_mid[$i-1][2] = $prediction_mid[$i-1][1]*$maxproc;			#approximate share of processors
	$prediction_mid[$i-1][$dim[0]] = $prediction_mid[$i-1][2]/$shorter;		#get ntasks_x or ntasks_y -- removefrag
	$prediction_mid[$i-1][$dim[1]] = $shorter;					#get ntasks_x or ntasks_y
	$prediction_mid[$i-1][5] = int($prediction_mid[$i-1][2]/$shorter);		#now get the int value
	$prediction_mid[$i-1][6] = (int(($prediction_mid[$i-1][$dim[0]]*10)))%10;	#for sorting to reduce wastage of processors

}

if($maxdom == 3) {
	$slabs = 1;
	&slab($slabs);
}

# x - $prediction_mid[$i-1][3] y - $prediction_mid[$i-1][4] start_rank - 7	

#sort the nests based on ratios - for 4,5
my @sorted = &sortnests(1);
@left = (0, 0, 0);
@right = (0, 0, 0);

if ($maxdom == 4) {


	if ($prediction_mid[$sorted[0]][1] + $prediction_mid[$sorted[1]][1] < 0.6) {
		$left[0] = $prediction_mid[$sorted[0]][1] + $prediction_mid[$sorted[1]][1];
		$right[0] = $prediction_mid[$sorted[2]][1];
	}
}
elsif ($maxdom == 5) {


	if ($prediction_mid[$sorted[0]][1] + $prediction_mid[$sorted[3]][1] < 0.6) {

		$left[0] = $prediction_mid[$sorted[0]][1] + $prediction_mid[$sorted[3]][1];
		$right[0] = $prediction_mid[$sorted[1]][1] + $prediction_mid[$sorted[2]][1];
	}
}

if ($maxdom == 4 or $maxdom == 5) {	#for both 5 and 6
	
	if($left[0] < 0.6 and $left[0] > 0) {

		$slabs = 0;

		$left[1] = $left[0]*$maxproc/$shorter;
		$right[1] = $right[0]*$maxproc/$shorter;

		$left[2] = int($left[1]);
		$right[2] = int($right[1]);

		$leftset = $left[2]*$shorter;

	}

	else {

		$slabs = 1;
		&slab($slabs);

	}
}


#stage 2 to get decomposition of processors

if ($maxdom == 4 and $slabs == 0) {

	$start_rank = 0;

	$prediction_mid[$sorted[0]][$dim[0]] = $left[2];						#x
	$prediction_mid[$sorted[0]][$dim[1]] = int(($prediction_mid[$sorted[0]][2]/$leftset) * $shorter);

	$prediction_mid[$sorted[1]][$dim[0]] = $left[2];		
	$prediction_mid[$sorted[1]][$dim[1]] = $shorter - $prediction_mid[$sorted[0]][$dim[1]];

	$prediction_mid[$sorted[2]][$dim[0]] = $right[2];
	$prediction_mid[$sorted[2]][$dim[1]] = $shorter;

	&partition($slabs);


}
elsif ($maxdom == 5 and $slabs == 0) {

	$start_rank = 0;	

	$prediction_mid[$sorted[0]][$dim[0]] = $left[2];							#x
	$prediction_mid[$sorted[0]][$dim[1]] = int(($prediction_mid[$sorted[0]][2]/$leftset) * $shorter);

	$prediction_mid[$sorted[3]][$dim[0]] = $left[2];		
	$prediction_mid[$sorted[3]][$dim[1]] = $shorter - $prediction_mid[$sorted[0]][$dim[1]];

	$prediction_mid[$sorted[2]][$dim[0]] = $right[2];
	$prediction_mid[$sorted[2]][$dim[1]] = int(($prediction_mid[$sorted[2]][2]/$leftset) * $shorter);

	$prediction_mid[$sorted[1]][$dim[0]] = $right[2];
	$prediction_mid[$sorted[1]][$dim[1]] = $shorter - $prediction_mid[$sorted[2]][$dim[1]];


	&partition($slabs);


}

$start_rank = 0;
foreach $i (1 .. $maxdom-1) {

	if($slabs == 0)	{				#vertical slabs suffice
	
#		$prediction_mid[$i-1][$dim[0]] = int($prediction_mid[$i-1][2]/$decomp[1]);		#old value
#		&removefrag($slabs, $dim[0]);

#		$prediction_mid[$i-1][4] = $longer;	#decomp[1];
#		$prediction_mid[$i-1][5] = $start_rank;
#		$start_rank += $prediction_mid[$i-1][3];


	}

}

&check();
&write();
exit;


#subroutines

sub check()
{
#Test to see if all processors have been used
	my $a = 0;
	foreach $i (1 .. $maxdom-1) {
		$a += $prediction_mid[$i-1][3] * $prediction_mid[$i-1][4];
	}

	if($a != $decomp[0]*$decomp[1])
	{
#		print "error\n";
	}
	if($a > $decomp[0]*$decomp[1])
	{
#		print "Major error\n";
	}
}


sub write()
{
	open OUT, ">$outfile" or die "ERROR: could not open file $outfile $!\n";
	printf (OUT "%d\n", $maxdom-1);

	foreach $i (1 .. $maxdom-1) {

		$id = $i + 1;					#nest id: skip 1 for the parent
			$start_rank = $prediction_mid[$i-1][7];
		$nest_x = $prediction_mid[$i-1][3];
		$nest_y = $prediction_mid[$i-1][4];

		printf (OUT "%d %d %d %d\n", $id, $start_rank, $nest_x, $nest_y); 

	}

	close (OUT);
}

sub get_dom()
{
	$maxdom = `grep max_dom $nlinp | awk '{print \$NF}' | awk -F',' '{print \$1}'`;
	chomp($maxdom);
	$e_we = `grep e_we $nlinp | awk -F'=' '{print \$NF}'`;
	chomp($e_we);
	$e_sn = `grep e_sn $nlinp | awk -F'=' '{print \$NF}'`;
	chomp($e_sn);
}

sub mpaspect()
{
	my $P = $_[0];
	my ($M, $N);
	my ($MINI, $MINM, $MINN, $PROCMIN_M, $PROCMIN_N) = (2*$P, 1, $P, 1, 1);

	for $M (1 .. $P)
	{
		if($P % $M == 0)
		{
			$N = $P / $M;
			if((abs($M-$N) < $MINI) and ($M > $PROCMIN_M) and ($N > $PROCMIN_N))
			{
				$MINI = abs($M-$N);
			        $MINM = $M;
				$MINN = $N;
			}
		}
	}

	$decomp[0] = $MINM;
	$decomp[1] = $MINN;

}


sub sortnests()
{
	my $option = $_[0];
	my (@temp, @prediction_324_copy);
	my $min_index;

	foreach $i (1 .. $maxdom-1) {
		$temp[$i-1] = $i-1;		#indices of nest ids
		$prediction_324_copy[$i-1][$option] = $prediction_mid[$i-1][$option];
	}

	foreach $i (1 .. $maxdom-1) {
		#print "\n$prediction_324_copy[$i-1][$option]";
	}

	my $tmp;

	foreach $i (0 .. $maxdom-2) {

	  $min_index = $i;

	  foreach $j ($i+1 .. $maxdom-2) { 	
		if ($prediction_324_copy[$min_index][$option] > $prediction_324_copy[$j][$option]) {
			$min_index = $j;
		}
	  }

	  $tmp = $temp[$i];
	  $temp[$i] = $temp[$min_index];
	  $temp[$min_index] = $tmp;

	  $tmp = $prediction_324_copy[$i];
	  $prediction_324_copy[$i] = $prediction_324_copy[$min_index];
	  $prediction_324_copy[$min_index] = $tmp;

	}

	return @temp;
}


sub removefrag()
{
	my $option = $_[0];
	my $field = $_[1];
	my $tmp = 0;
	
	if($option == 1) {			#option 1 is for maxdom 5
		foreach $i (1 .. $maxdom-1) {
			$tmp += int($prediction_mid[$i-1][$field]);
		}
	}
	else {
		$tmp = $left[2] + $right[2];
	}

	my $remaining = $longer - $tmp;

	if ($remaining > 0)
	{
		if($option == 1) {		#option 1 is for maxdom 5
			my @sorted = &sortnests(6);		#sort based on first decimal place
			for ($j=$#sorted ; $j>=0 ; $j--) {
				$prediction_mid[$sorted[$j]][$field] += 1;
				if(--$remaining == 0) { last; }
			}
		}
		else {				#subpartition the partitions
			my $t1 = ($left[1]*10)%10;
			my $t2 = ($right[1]*10)%10;

			if ($t1 > $t2) {
				$left[2] += 1;
			}
			else {
				$right[2] += 1;
			}
		}
	}

	foreach $i (1 .. $maxdom-1) {
	}

}


sub trim($)
{
	my $string = shift;
	chomp($string);
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub indexArray($@)
{
	my $s=shift;
	$_ eq $s && return @_ while $_=pop;
	-1;
}

sub slab()
{

	my $opt = $_[0];

	foreach $i (1 .. $maxdom-1) {
		$prediction_mid[$i-1][$dim[0]] = int($prediction_mid[$i-1][$dim[0]]);
	}

	&removefrag($opt, $dim[0]);

	my $start_rank = 0;
	foreach $i (1 .. $maxdom-1) {

		$prediction_mid[$i-1][7] = $start_rank;
		if ($xgreatery == 1) {
			$start_rank += $prediction_mid[$i-1][3] * $prediction_mid[$i-1][4];	#x * y		-- horizontal slabs 
		}
		else {
			$start_rank += $prediction_mid[$i-1][3];
		}

	}

	foreach $i (1 .. $maxdom-1) {
		$prediction_mid[$i-1][$dim[0]] = int($prediction_mid[$i-1][$dim[0]]);
	}

}

sub partition()
{

	my $opt = $_[0];

	&removefrag($opt, $dim[0]);


	if($maxdom == 4) {

		$prediction_mid[$sorted[0]][$dim[0]] = $left[2];					#x
		$prediction_mid[$sorted[0]][7] = $start_rank;					#start_rank
		$prediction_mid[$sorted[1]][$dim[0]] = $left[2];		
		$prediction_mid[$sorted[2]][$dim[0]] = $right[2];

		if($xgreatery == 1) {
			$prediction_mid[$sorted[1]][7] = $prediction_mid[$sorted[0]][3]; 	#xdim of sorted[0]
			$prediction_mid[$sorted[2]][7] = $left[2] * $shorter;
		}
		else {
			$prediction_mid[$sorted[1]][7] = $prediction_mid[$sorted[0]][$dim[1]] * $shorter; #start_rank nest 1 y depth * parent x width
			$prediction_mid[$sorted[2]][7] = $left[2];
		}

	}
	else {

		$prediction_mid[$sorted[0]][$dim[0]] = $left[2];							#x
		$prediction_mid[$sorted[3]][$dim[0]] = $left[2];		
		$prediction_mid[$sorted[1]][$dim[0]] = $right[2];
		$prediction_mid[$sorted[2]][$dim[0]] = $right[2];

		if($xgreatery == 1) {

			$prediction_mid[$sorted[0]][7] = $start_rank;					#start_rank
			$prediction_mid[$sorted[3]][7] = $prediction_mid[$sorted[0]][3];		#start_rank nest 1 y depth * parent x width
			$prediction_mid[$sorted[1]][7] = $left[2] * $shorter;				#left width
			$prediction_mid[$sorted[2]][7] = $prediction_mid[$sorted[1]][7] + $prediction_mid[$sorted[1]][3];
		}
		else {
			
			$prediction_mid[$sorted[0]][7] = $start_rank;					#start_rank
			$prediction_mid[$sorted[3]][7] = $prediction_mid[$sorted[0]][4] * $shorter;	#start_rank nest 1 y depth * parent x width
			$prediction_mid[$sorted[1]][7] = $left[2];					#left width
			$prediction_mid[$sorted[2]][7] = $left[2] + $prediction_mid[$sorted[1]][4] * $shorter;

		}

	}

}	

