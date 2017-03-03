#!/usr/bin/perl
# ##############################################################################
# NLStradamus FOR DISTRIBUTION
# Copyright Alex Nguyen Ba 2008
# Version 1.8
#
# LICENCE 
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ##############################################################################

use threads;
use threads::shared;

# USAGE

die "usage: -i fastafile ... or try -h for more information\n" unless (@ARGV);
my $fastafile;
my ($posterior_threshold) = 0.6;
my ($model) = 1; 
my $QUIET=1;
my $ALGORITHM=1;
my $cpu = 1;
my $TAB = 0;
for my $i (0 .. scalar(@ARGV)-1) { 
	if ($ARGV[$i] eq "-i") { $fastafile=$ARGV[$i+1]; }
	elsif ($ARGV[$i] =~m/\-h/) { help(); }
	elsif ($ARGV[$i] eq "-t") { 
		if($ARGV[$i+1] > 0 && $ARGV[$i+1] <= 1){
			$posterior_threshold = $ARGV[$i+1]; 
		}
		else{
			print STDERR "Posterior Threshold must be between 0 and 1\n";
			exit;
		}
	}
	elsif ($ARGV[$i] eq "-a"){
		if($ARGV[$i+1] == 2){
			$ALGORITHM = 2;	
		}
		elsif($ARGV[$i+1] == 1){
			$ALGORITHM = 1;
		}
		elsif($ARGV[$i+1] == 0){
			$ALGORITHM =0;
		}
		else{
			print STDERR "ALGORITHM must be 0,1 or 2\n";
			exit;
		}
	}
	elsif ($ARGV[$i] eq "-tab"){$TAB = 1;}
	elsif ($ARGV[$i] =~m/\-q/){ $QUIET=1; }
	elsif ($ARGV[$i] =~m/\-m/){
		if($ARGV[$i+1] == 2){
			$model = 2;	
		}
		elsif($ARGV[$i+1] == 1){
			$model = 1;
		}
		else{
			print STDERR "Model must be 1 or 2\n";
			exit;
		}
	}
	elsif ($ARGV[$i] eq "-cpu"){
		if($ARGV[$i+1] < 1){
			$cpu = 1;
		}
		else{
			$cpu = $ARGV[$i+1];
		}
	}
}	
if($TAB ==1 ){
	$QUIET = 1;	
}
################################################################################
# Parameters
# 
# Parameters for the HMM are arranged as follows: 
# $e{state}{'emission'} = probability;
# $a[state_1][state_2] = probability;
# $begin[state] = probability;
#
# Emission probabilities ($e) are the probabilities of emitting something that 
# is observable.
#
# Transition probabilities ($a) are the probabilities of going from one state to
# another. The synthax is such : $a[from][to] = $probability.
#
# Beginning probabilities are the probabilities of being into the state at the
# first amino acid.
#
# STATES ARE ZERO-INDEX
################################################################################

# Parameters

## Yeast Genome Frequency
$e{0}{"L"} = 0.0950300339227636;
$e{0}{"S"} = 0.0898291379118776;
$e{0}{"K"} = 0.0734833177642526;
$e{0}{"I"} = 0.0655547462377008;
$e{0}{"E"} = 0.0654254611397116;
$e{0}{"N"} = 0.0616309263671642;
$e{0}{"T"} = 0.0591563341467664;
$e{0}{"D"} = 0.0585263693589517;
$e{0}{"V"} = 0.0556015083490053;
$e{0}{"A"} = 0.0549499388896433;
$e{0}{"G"} = 0.0497648177182998;
$e{0}{"R"} = 0.0443626210376004;
$e{0}{"F"} = 0.0441551476044877;
$e{0}{"P"} = 0.0437563185090993;
$e{0}{"Q"} = 0.0395581536030419;
$e{0}{"Y"} = 0.0337842193992118;
$e{0}{"H"} = 0.0216513102033034;
$e{0}{"M"} = 0.0208056416313105;
$e{0}{"C"} = 0.0125821491915738;
$e{0}{"W"} = 0.0103918470142344;

$e{0}{"B"} = $e{0}{"D"} + $e{0}{"N"};
$e{0}{"X"} = 1;
$e{0}{"Z"} = $e{0}{"E"} + $e{0}{"Q"};

## State Frequencies

$e{1}{'L'} = 0.0548387096774;
$e{1}{'S'} = 0.0677419354839;
$e{1}{'K'} = 0.270967741935;
$e{1}{'I'} = 0.041935483871;
$e{1}{'E'} = 0.041935483871;
$e{1}{'N'} = 0.041935483871;
$e{1}{'T'} = 0.0306451612903;
$e{1}{'D'} = 0.0258064516129;
$e{1}{'V'} = 0.0290322580645;
$e{1}{'A'} = 0.0516129032258;
$e{1}{'G'} = 0.0564516129032;
$e{1}{'R'} = 0.133870967742;
$e{1}{'F'} = 0.0209677419355;
$e{1}{'P'} = 0.0564516129032;
$e{1}{'Q'} = 0.0225806451613;
$e{1}{'Y'} = 0.0112903225806;
$e{1}{'H'} = 0.0225806451613;
$e{1}{'M'} = 0.0112903225806;
$e{1}{'C'} = 0.00161290322581;
$e{1}{'W'} = 0.00645161290323;

$e{1}{"B"} = $e{1}{"D"} + $e{1}{"N"};
$e{1}{"X"} = 1;
$e{1}{"Z"} = $e{1}{"E"} + $e{1}{"Q"};

if($model == 2){
	
	$e{2}{"L"} = 0.0950300339227636;
	$e{2}{"S"} = 0.0898291379118776;
	$e{2}{"K"} = 0.0734833177642526;
	$e{2}{"I"} = 0.0655547462377008;
	$e{2}{"E"} = 0.0654254611397116;
	$e{2}{"N"} = 0.0616309263671642;
	$e{2}{"T"} = 0.0591563341467664;
	$e{2}{"D"} = 0.0585263693589517;
	$e{2}{"V"} = 0.0556015083490053;
	$e{2}{"A"} = 0.0549499388896433;
	$e{2}{"G"} = 0.0497648177182998;
	$e{2}{"R"} = 0.0443626210376004;
	$e{2}{"F"} = 0.0441551476044877;
	$e{2}{"P"} = 0.0437563185090993;
	$e{2}{"Q"} = 0.0395581536030419;
	$e{2}{"Y"} = 0.0337842193992118;
	$e{2}{"H"} = 0.0216513102033034;
	$e{2}{"M"} = 0.0208056416313105;
	$e{2}{"C"} = 0.0125821491915738;
	$e{2}{"W"} = 0.0103918470142344;
	
	$e{2}{"B"} = $e{2}{"D"} + $e{2}{"N"};
	$e{2}{"X"} = 1;
	$e{2}{"Z"} = $e{2}{"E"} + $e{2}{"Q"};

	$e{3}{'L'} = 0.0548387096774;
	$e{3}{'S'} = 0.0677419354839;
	$e{3}{'K'} = 0.270967741935;
	$e{3}{'I'} = 0.041935483871;
	$e{3}{'E'} = 0.041935483871;
	$e{3}{'N'} = 0.041935483871;
	$e{3}{'T'} = 0.0306451612903;
	$e{3}{'D'} = 0.0258064516129;
	$e{3}{'V'} = 0.0290322580645;
	$e{3}{'A'} = 0.0516129032258;
	$e{3}{'G'} = 0.0564516129032;
	$e{3}{'R'} = 0.133870967742;
	$e{3}{'F'} = 0.0209677419355;
	$e{3}{'P'} = 0.0564516129032;
	$e{3}{'Q'} = 0.0225806451613;
	$e{3}{'Y'} = 0.0112903225806;
	$e{3}{'H'} = 0.0225806451613;
	$e{3}{'M'} = 0.0112903225806;
	$e{3}{'C'} = 0.00161290322581;
	$e{3}{'W'} = 0.00645161290323;
	
	$e{3}{"B"} = $e{3}{"D"} + $e{3}{"N"};
	$e{3}{"X"} = 1;
	$e{3}{"Z"} = $e{3}{"E"} + $e{3}{"Q"};
}
if($model == 1){
	## Transition Frequencies
	($a[0][1]) = 0.00263746344819678;
	($a[0][0]) = 1-$a[0][1];
	
	($a[1][0]) = 0.0741935483870968;
	($a[1][1]) = 1 - $a[1][0];
	
	## Initiation probabilities
	($begin[0]) = $a[0][0];
	($begin[1]) = 1 - $begin[0];
}
elsif($model == 2){
	## Transition Frequencies
	($a[0][1]) = 0.00263746344819678;
	($a[0][0]) = 1-$a[0][1];
	($a[0][2]) = 0;
	($a[0][3]) = 0;
	
	($a[1][0]) = 0;
	($a[1][2]) = 0.148387096;
	($a[1][1]) = 1-$a[1][2];
	($a[1][3]) = 0;
	
	($a[2][0]) = 0;
	($a[2][1]) = 0;
	($a[2][2]) = 0.88028169;
	($a[2][3]) = 1-$a[2][2];
	
	($a[3][1]) = 0;
	($a[3][2]) = 0;
	($a[3][0]) = 0.148387096;
	($a[3][3]) = 1-$a[3][0];
	
	## Initiation probabilities
	($begin[0]) = $a[0][0];
	($begin[1]) = 1 - $begin[0];
	($begin[2]) = 0;
	($begin[3]) = 0;
}

################################################################################
# Handle the input fasta file.
################################################################################

open (ORFS, $fastafile) or die "Can't open the fasta file : $!\n";
my ($protein_count) = -1;
my ($fastatext) = "";
while($lineorf = <ORFS>){
	chomp($lineorf);
	$lineorf =~ s/\r//;
	if(index($lineorf,">") == "-1"){
		$protein[$protein_count]{'sequence'} .= uc($lineorf);
		$protein[$protein_count]{'sequence'} =~ s/[^ABCDEFGHIKLMNPQRSTVWXYZ]//g;	
	}
	else{
		@annotation = split(" ",$lineorf);
					
		my($gene_ID) = $annotation[0];
		$gene_ID =~ s/>//g;
		
		++$protein_count;
		$protein[$protein_count]{'ID'} = $gene_ID;
	}
}
close(ORFS);

################################################################################
# Process
################################################################################

my $total_site_count = 0;
share($total_site_count);
my $spread_process = int($protein_count / $cpu);
if($TAB == 1){
	print "#ID\talgorithm\tscore\tstart\tstop\tsequence\n";	
}
for(my ($p) = 0;$p < ($protein_count + 1);$p=$p+$spread_process+1){
	$threads[$p] = threads->create(\&run_hmm_multiple,$p .. min($p+$spread_process,$protein_count));	
}
for(my ($p) = 0;$p < ($protein_count + 1);$p=$p+$spread_process+1){
	print $threads[$p]->join();
}
if($TAB != 1){
	print "===================================================";
	print "\n";
	print "Analyzed ";
	print ($protein_count+1);
	print " proteins.";
	print "\n";
	print "$total_site_count sites were found ";
	if($ALGORITHM != 0){
		print "using the posterior probability threshold.\n";	
	}
	else{
		print "using the viterbi path.\n";	
	}
	print "Input file : $fastafile.\n";
	print "Threshold used : $posterior_threshold.\n";
	print "===================================================";
	print "\n";
}
sub help {
	print STDERR "You are using NLStradamus v1.8 copyright Alex Nguyen Ba 2011\n";
	print STDERR "-i input file\n";
	print STDERR "-t [optional] Posterior Threshold (0...1) default 0.6\n";
	print STDERR "-m [optional] Model (1 for two state, 2 for four state) default 1\n";
	print STDERR "-a [optional] Algorithm (0 for viterbi, 1 for posterior with threshold, 2 for both) default 1\n";
	print STDERR "-q [optional] quiet mode, default prints sequences to screen\n";
	print STDERR "-cpu [optional] multithread capability (1...N) default is 1\n";
	print STDERR "-tab [optional] flag for tab delimited output. default is off\n";
	print STDERR "Please read the README.txt file for an example\n";
	exit;	
}
sub max {
	my($max) = shift(@_);

	foreach $temp (@_) {
		$max = $temp if $temp > $max;
	}
	return($max);
}
sub min {
	my($min) = shift(@_);

	foreach $temp (@_) {
		$min = $temp if $temp < $min;
	}
	return($min);
}
sub run_hmm_multiple{
	my @print_array;
	foreach $p (@_){
		my @gene_txt = run_hmm($p);
		push(@print_array,@gene_txt);	
		
	}
	return join("",@print_array);
}
sub run_hmm{
	my @print_array;
	my $p = @_[0];

	unless ($QUIET) {
		push(@print_array,$protein[$p]{'ID'});
		push(@print_array,"\n");
		push(@print_array,$protein[$p]{'sequence'});
		push(@print_array,"\n");
	}
	my (@v);
	my (@f);
	my (@b);
	
	my (@log_fs);
	my (@log_bs);
	
	my(@f_scaling);
	
	my(@prev_v);
	my(@prev_f);
	
	$count = @a;
	#Initiation of Viterbi Algorithm
		
	for(my ($i) = 0;defined $a[$i];++$i){
		$v[$i][0] = log(1);	
	}
		
	#Initiation of Forward Algorithm
		
	for(my $i = 0;defined $a[$i];++$i){
		$f[$i][0] = 1;
	}
		
	#Process initiation
	my ($log_fs_process) = 0;
	my ($log_bs_process) = 0;
	
	if(length($protein[$p]{'sequence'}) < 1){
		next;	
	}
	
	#Processing of Viterbi and Forward Algorithm : i = amino acids;
	for($i = 1;$i <= length($protein[$p]{'sequence'});++$i){
			
		#Fetch letters
		my ($next_letter) = substr($protein[$p]{'sequence'},$i-1,1);

		#Set Forward Sum Scaling Value
		if($i == 1){
			for(my ($k) = 0;defined $begin[$k];++$k){
				$f_scaling[$i] += ($e{$k}{$next_letter} * $begin[$k]);
			}
		}
		else{
			for(my ($k) = 0;defined $e{$k};++$k){
				for(my ($j) = 0;defined $e{$j};++$j){
					$f_scaling[$i] += $f[$j][$i-1] * $a[$j][$k] * $e{$k}{$next_letter}
				}
			}			
		}
		
		#Go through all the states : j = states. Browse horizontally.
		for(my ($j) = 0;defined $a[$j];++$j){
			
			#Set Initial values for Forward and Viterbi.
			if($i == 1){
				if($begin[$j] != 0){
					$prev_v[$j][$j] = $v[$j][$i-1] + log($begin[$j]); #$v[$k][$i-1] == log(1) == 0;
				}
				else{
					$prev_v[$j][$j] = -10**16;
				}
				
				$prev_f[$j][$j] = $begin[$j];
			}
			else{
				for(my ($k) = 0;defined $a[$k];++$k){
					if($a[$k][$j] != 0){
						$prev_v[$k][$j] = $v[$k][$i-1] + log($a[$k][$j]);
					}
					else{
						$prev_v[$k][$j] = -10**16;	
					}
				}
				for(my ($k) = 0;(defined $a[$k]);++$k){
					$prev_f[$k][$j] = $f[$k][$i-1] * $a[$k][$j];	
				}
			}
			
			#Viterbi Maximisation
			my (@new_array);
			@new_array = ();
			my ($increment) = 0;
			for(my ($k) = 0;defined $a[$k];++$k){
				if(defined $prev_v[$k][$j]){
					$new_array[$increment] = $prev_v[$k][$j];
					++$increment;
				}
			}

			$v[$j][$i] = max(@new_array);
			$v[$j][$i] += log($e{$j}{$next_letter});

			#Forward Sum
			for(my ($k) = 0;(defined $a[$k]);++$k){
				if(defined $prev_f[$k][$j]){
					$f[$j][$i] += $prev_f[$k][$j];
				}
			}
					
			$f[$j][$i] *= $e{$j}{$next_letter}/$f_scaling[$i];
		}
		
		#Scaling value log sum.
		$log_fs_process += log($f_scaling[$i]);
		$log_fs[$i] = $log_fs_process;
		
	}

	# Traceback Initiation of Viterbi
	my @newarray;
	my $maxkey;
	for(my ($j) = 0;defined $a[$j];++$j){
		$newarray[$j] = $v[$j][length($protein[$p]{'sequence'})];
		if(max (@newarray) == $v[$j][length($protein[$p]{'sequence'})]){
			$maxkey = $j;
		}
	}
	$step[length($protein[$p]{'sequence'})] = $maxkey;
	
	# Traceback of Viterbi
	
	for($i = length($protein[$p]{'sequence'});$i >= 1;--$i){
		my (@traces);
		my $max_traces;
		for(my ($k) = 0;defined $a[$k];++$k){		
			if($a[$k][$step[$i]] != 0){
				$traces[$k] = $v[$k][$i-1] + log($a[$k][$step[$i]]);
			}
			else{
				$traces[$k] = -10**16;	
			}
			if(max(@traces) == $traces[$k]){
				$max_traces = $k;
			}
		}
		$step[$i-1] = $max_traces;		
	}
	
	# Initiation of Backward Algorithm [ Scaled ]
	
	for(my ($k) = 0;defined $a[$k];++$k){
		$b[$k][length($protein[$p]{'sequence'})] = 1/$f_scaling[length($protein[$p]{'sequence'})];	
	}
	
	# Implementation of Backward Algorithm
	
	for($i = (length($protein[$p]{'sequence'})-1);$i >= 1;$i = $i-1){
		
		#Scaling values of b
		$log_bs_process += log($f_scaling[$i+1]);
		$log_bs[$i+1] = $log_bs_process;
		
		for(my ($k) = 0;defined $a[$k];++$k){
			for(my ($j) = 0;defined $a[$j];++$j){
				$b[$k][$i] += ($b[$j][$i+1] * $e{$j}{substr($protein[$p]{'sequence'},$i,1)} * $a[$k][$j])/$f_scaling[$i];
			}
		}
	}
	
	#Termination of Scaling Values
	$log_bs_process += log($f_scaling[1]);
	$log_bs[1] = $log_bs_process;
	my $pmatchn=0;
	if ($ALGORITHM!=0) {
		#Termination step of Forward Algorithm [ Scaled Value ]
		my ($f_P);
		for($k = 0;defined $f[$k];++$k){
			$f_P += $f[$k][length($protein[$p]{'sequence'})];	
		}
	
		#Termination step of Backward Algorithm [ Scaled Value ]
		
		my ($b_P);
		for($k = 0;defined $b[$k];++$k){
			$b_P += $e{$k}{substr($protein[$p]{'sequence'},0,1)} * $b[$k][1] * $begin[$k];
		}
		#Posterior in Log Space
		
		my ($posterior_string) == "";
		for($i = 1;$i <= length($protein[$p]{'sequence'});$i = $i + 1){
			#$posterior_state = which posterior state to threshold.
			my $posterior_state;
			if($model == 1){
				$posterior_state = 1;
			}
			else{
				$posterior_state = 0;	
			}
			if($f[$posterior_state][$i] == 0 || $b[$posterior_state][$i] == 0){
				$posterior[$i] = -10**16;
			}
			else{
				$posterior[$i] = log($f[$posterior_state][$i]) + $log_fs[$i] + log($b[$posterior_state][$i]) + $log_bs[$i] - $log_fs[length($protein[$p]{'sequence'})];
			}
			if($model == 1){
				if(exp($posterior[$i]) > $posterior_threshold){	
					push(@print_array,"1") unless ($QUIET);
					$posterior_string .= "1";
				}
				else{
					push(@print_array,"0") unless ($QUIET);
					$posterior_string .= "0";
				}
			}
			else{
				if(1-exp($posterior[$i]) > $posterior_threshold){
					push(@print_array,"1") unless ($QUIET);
					$posterior_string .= "1";
				}
				else{
					push(@print_array,"0") unless ($QUIET);
					$posterior_string .= "0";
				}
			}

		}
		push(@print_array,"\n") unless ($QUIET);
		my %pstarts; my %pstops; 
		while ($posterior_string=~m/([1]+)/g) { 
			$pmatchn++; my $test=$1;
			my $pos=pos($posterior_string); #print "$pos $test??",  length($test),"\n";
			$pstarts{$pmatchn}= $pos-length($test);
			$pstops{$pmatchn}=$pstarts{$pmatchn}+length($test);
		}
		for my $m (1 .. $pmatchn) {
			my $max_posterior = 0;
			for(my $i = $pstarts{$m}+1;$i < $pstops{$m}+1;++$i){
				if($model == 1){
					if(exp($posterior[$i]) > $max_posterior){
						$max_posterior = exp($posterior[$i]);
					}
				}
				else{
					if(1-exp($posterior[$i]) > $max_posterior){
						$max_posterior = 1-exp($posterior[$i]);
					}
				}
			}

			push(@print_array,$protein[$p]{'ID'}."\tposterior\t".sprintf("%.3f",$max_posterior)."\t".($pstarts{$m}+1)."\t".($pstops{$m})."\t".substr($protein[$p]{'sequence'},$pstarts{$m},$pstops{$m}-$pstarts{$m})."\n");
			++$total_site_count;		
		}
	}
	#Print prediction of Viterbi
	my $vmatchn=0;
	if ($ALGORITHM!=1) {
		my $viterbi_string; 
		for(my ($i) = 1;$i<=length($protein[$p]{'sequence'});++$i){
			$viterbi_string .= $step[$i];
			push(@print_array,$step[$i]) unless ($QUIET);
			
		}
			
		print "\n" unless ($QUIET);
		
		my %vstarts; my %vstops;
		while ($viterbi_string=~m/([123]+)/g) { 
			$vmatchn++; my $test=$1;
			my $pos=pos($viterbi_string); #print "$pos $test??",  length($test),"\n";
			$vstarts{$vmatchn}= $pos-length($test);
			$vstops{$vmatchn}=$vstarts{$vmatchn}+length($test);
		}
		for my $m (1 .. $vmatchn) {
			push(@print_array,$protein[$p]{'ID'}."\tviterbi\t\t".($vstarts{$m}+1). "\t".$vstops{$m}."\t".substr($protein[$p]{'sequence'},$vstarts{$m},$vstops{$m}-$vstarts{$m})."\n");
		}
		if($ALGORITHM == 0){
			++$total_site_count;	
		}
	}
	
	if($QUIET && $pmatchn == 0 || $TAB == 1){

	}
	else{
		if($ALGORITHM == 0){
			push(@print_array,"Finished analyzing ".$protein[$p]{'ID'}.". Found $vmatchn sites.\n\n");
		}
		else{
			push(@print_array,"Finished analyzing ".$protein[$p]{'ID'}.". Found $pmatchn sites.\n\n");
		}
	}
	
	
	return @print_array;
}

