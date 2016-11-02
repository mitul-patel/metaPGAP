#!/usr/bin/perl -w

use strict 'subs';
use POSIX;
use Getopt::Std;

# Usage statement
die "Usage 1: $0 [-s|c|b|n] [-f #] inputfile outputfile
Reformats existing COG file.
-or-
Usage 2: $0 -r [-s|c|b|n] [-f #] COGtriangles_params
Runs COGtriangles first with given parameters, then reformats COG file.

Modifies COGtriangles output to restrict proteins to only appear once (on one line) in the outputfile:

	-r	run COGtriangles first with the given parameters, moving the output file
		(-o <filename> option to COGtriangles) to <filename>.orig and replacing <filename> with the modified output.

	-s	strict mode, removes all proteins that appear to belong to multiple COGs, keeping only their first occurance.

	-c	combine mode, combine all COGs that a protein belongs to into the same line, using the special prefix 'multi:'.
		Example - the three lines:
			9626516,Acholeplasma_phage_L2_uid14066,9626516,289,1,289,cluster00001,
			9626516,Acholeplasma_phage_L2_uid14066,9626516,289,1,289,cluster00002,
			9626516,Acholeplasma_phage_L2_uid14066,9626516,289,1,289,cluster02262,
		will get combined into a single entry:
			9626516,Acholeplasma_phage_L2_uid14066,9626516,289,1,289,multi:cluster00001:cluster00002:cluster02262,

	-b	[default] special combine mode, which like the regular one, will not report COGs that are a subset of other,
		larger COGs, but differs in that a COG will only not be reported if it is a complete subset of another single
		COG, but will be reported if it appears joined to several different COGs.

	-n	no change mode, does not remove any proteins. This is the ideal option, if you can handle the additional
		complexity, since proteins appearing in multiple COGs can occur due to various errors: domain shuffling/
		recombination, unresolved paralogy, missing homology edges, etc., and this option preserves such cases,
		which may be useful for diagnostic purposes.

	-f #	filter results to only print clusters of # size or greater. Cluster size refers to # of organisms, not proteins.
		Examples:
			'-f 3' will print only true COGs of 3 or more members.
			'-f 2' will print true COGs plus SymBet pairs.
			'-f 1' [default behavior of no filtering] will print true COGs plus SymBet pairs plus singlets that
				belong to no COGs (note that some proteins in the input set may still not appear here due to
				various filters internal to the COGtriangles program).
			'-f 10' will print only large COGs containing 10 or more organisms.

In any mode, clusters (COGs or SymBet pairs) that are complete subsets of another larger COG will be removed from the output.

" if scalar(@ARGV)<2;

# ---------------------------------
# Get & process common input params
# ---------------------------------

# Run COGtriangles?
$_=$ARGV[0];
if ($_ eq "-r")
{
	$runcogtriangles=1; # run COGtriangles
	shift @ARGV;
	$_=$ARGV[0];
}
else
{
	$runcogtriangles=0; # do not run COGtriangles, default
}

# Which mode am I modifying results in?
if ($_ eq "-s")
{
	$mode="s"; # set 'strict' mode
	shift @ARGV;
}
elsif ($_ eq "-b")
{
	$mode="b"; # set 'special combine' mode
	shift @ARGV;
}
elsif ($_ eq "-c")
{
	$mode="c"; # set 'combine' mode
	shift @ARGV;
}
elsif ($_ eq "-n")
{
	$mode="n"; # set 'no change' mode
	shift @ARGV;
}
else
{
	$mode="b"; # set 'special combine' mode by default
}

# Set filter level?
$_=$ARGV[0];
if ($_ eq "-f")
{
	shift @ARGV;
	$_=shift @ARGV;
	die "Invalid # to filter option: '$_'. Expected integer value.\n" if $_!~/^\d+$/;
	$filter=$_; # set filtering level
}
else
{
	$filter=1; # no filtering, default
}

# -------------------------------------------------------------
# Either run COGtriangles or else set up to just read its input
# -------------------------------------------------------------
$_=$ARGV[0];
if ($runcogtriangles==0)
{ # If not running COGtriangles, get input & output file names
	$filename=shift @ARGV;
	open (IN, "<$filename") or die "Error - cannot read file '$filename' : $!\n";
	die "Error - no output file provided!\n" if scalar(@ARGV)==0;
	$filename=shift @ARGV;
	die "Error - extra parameters provided to $0 - why?\n" if scalar(@ARGV)>0;
}
else
{ # If running COGtriangles, get its params, run it now, and make a backup copy of its output
	my $cmd=$0;
	$cmd=~/(.*)\/\S+$/;
	my $cmdpath=$1;
	$filename="";
	for (my $i=0; $i<scalar(@ARGV); $i++)
	{
	        if ($ARGV[$i]=~/-([iqlon])=(.*)/)
	        {
	                $ARGV[$i]="-$1=\"$2\"";
	                if ($1 eq "o")
	                {
	                        $filename=$2;
	                }
	        }
	}
	$filename="out.csv" if $filename eq "";
	my $args=join(' ',@ARGV);

	print("Running: $cmdpath/COGtriangles $args\n");
	system("$cmdpath/COGtriangles $args") && do { print "Error - quitting.\n"; exit $?; };
	system("/bin/mv -f $filename $filename.orig") && die "Error - could not save COGtriangles output file '$filename' to '$filename.orig' : $!\n";
	open (IN, "<$filename.orig") or die "Error - cannot read file '$filename.orig' : $!\n";
}
open (OUT, ">$filename") or die "Error - cannot output to file '$filename' : $!\n";

# ---------------------------------------------------------
# Now do the actual modification of the COGtriangles output
# ---------------------------------------------------------
# Read first line and identify text prefix that COGs are referred to by ("cluster", "CLS", "COG", etc.)
my $prefix;
my $startingclusternumber;
{
	$_=<IN>;
	chomp;
	my @array=split /,/;
	die "Error in input file at line #1, invalid format. Received \"$_\"\n" if !defined($array[6]);
	$prefix=$array[6];
	$prefix=~s/\d+$//;
	$prefix=~s/.*://;
	$array[6]=~/(\d+)$/;
	$startingclusternumber=$1;
	seek(IN,0,0) or die "Error - could not seek!!\n";
};

# Now read the actual list of proteins that belong to which COGs
while ($line=<IN>)
{
	chomp($line);
	my @array=split /,/, $line;

	my $prot=$array[0];
	my $org=$array[1];
	my $cog=$array[6];

	$prot2org{$prot}=$org;
	if (!defined($cog))
	{
		$singlets{$prot}=$line;
	}
	else
	{
		$cog=~s/$prefix//;
		die "Error - invalid COGtriangles input file. Multiple lines per protein per COG found at line '$line'.\n" if defined($line{$prot}{$cog});
		$line{$prot}{$cog}=$line;
		$cogorgs{$cog}{$org}++;
		push @{$cogprots{$cog}}, $prot;
	}
}
close IN;

# Sort COGs by size: number of organisms, and then number of proteins
@cogs=sort { scalar(keys %{$cogorgs{$b}}) <=> scalar(keys %{$cogorgs{$a}}) 
	|| scalar(@{$cogprots{$b}}) <=> scalar(@{$cogprots{$a}}) 
	|| $a <=> $b } keys %cogorgs;

# Iterating through old COGs, make a set of new COGs according to the new criteria (mode & filtering)
for (my $i=0; $i<scalar(@cogs); $i++)
{
	# count organisms within each COG
	# Mode "n" : count all
	# Mode "b" : count all
	# Mode "c" : count all, and make special note of whether at least one unprocessed protein is present
	# Mode "s" : count only if not seen yet (in a previous, larger COG)
	my %countorgs;
	my $countnewprots=0;
	foreach my $prot (@{$cogprots{$cogs[$i]}})
	{
		next if $mode eq "s" && defined($processedprotein{$prot});
		$countorgs{$prot2org{$prot}}++;
		$countnewprots++ if !defined($processedprotein{$prot});
	}

	# Do not retain COG/cluster if:
	#   size does not pass the given threshold (any mode), or
	#   mode "c" finds no proteins not already present in other, larger COGs (possibly shared between more than one of them), or
	#   mode "s" demotes a COG to a singlet, or
	#   mode "b" finds that this COG is a strict subset of some other single COG (i.e., all proteins also belong to the same other COG)
	next if scalar(keys %countorgs)<$filter;
	next if $mode eq "c" && $countnewprots==0;
	if ($mode eq "s" && $filter<2 && scalar(keys %countorgs)<2)
	{
		foreach my $prot (@{$cogprots{$cogs[$i]}})
		{
			next if defined($processedprotein{$prot});
			my $line=$line{$prot}{$cogs[$i]};
			$line=~s/,$prefix\d+,/,,/;
			die "Error - how is singlet $prot already defined!?\n" if defined($singlets{$prot});
			$singlets{$prot}=$line;
			$processedprotein{$prot}++;
		}
		next;
	}
	if ($mode eq "b" && $countnewprots==0)
	{
		my %cogseen;
		foreach my $prot (@{$cogprots{$cogs[$i]}})
		{
			foreach my $cog (@{$protinnewcog{$prot}})
			{
				$cogseen{$cog}++;
			}
		}
		my @tmp=sort { $cogseen{$b}<=>$cogseen{$a} } keys %cogseen;
		# if all proteins in this COG also belong to the same cog (that is already seen, and hence larger), then do not report this as a distinct COG
		next if scalar(@tmp)>0 && $cogseen{$tmp[0]} == scalar(@{$cogprots{$cogs[$i]}});

		# If COG is accepted in this mode, go ahead and make it (so it can be sorted by size), but still skip printing it off separately (except in multi: lines)
		$skipthiscog{$cogs[$i]}=1;
	}

	# Now that filtering has been applied & this COG/cluster accepted:
	# Mode "n" : accept all proteins
	# Mode "b" : accept all proteins, but when they belong to multiple COGs, only report them with the first (largest) one, although keep a record of all COGs that they belong to
	# Mode "c" : accept all proteins, but when they belong to multiple COGs, only report them with the first (largest) one, although keep a record of all COGs that they belong to
	# Mode "s" : accept only proteins that do not already belong to some other, larger COG
	$cogmapping{$cogs[$i]}=$i;
	foreach my $prot (@{$cogprots{$cogs[$i]}})
	{
		next if $mode eq "s" && defined($processedprotein{$prot});
		if (($mode eq "c" || ($mode eq "b" && $countnewprots>0)) && defined($processedprotein{$prot}))
		{
			push @{$protinnewcog{$prot}}, $cogs[$i];
		}
		else
		{
			$newcogorgs{$cogs[$i]}{$prot2org{$prot}}++;
			push @{$newcogprots{$cogs[$i]}}, $prot;
			push @{$protinnewcog{$prot}}, $cogs[$i];
			$processedprotein{$prot}++;
		}
	}
}

# Calculate amount to pad new COG#s
$numpad=1+floor(log(scalar(keys %newcogorgs))/log(10));

# Sort new COGs again by size: number of organisms, and then number of proteins
@newcogs=sort { scalar(keys %{$newcogorgs{$b}}) <=> scalar(keys %{$newcogorgs{$a}}) 
	|| scalar(@{$newcogprots{$b}}) <=> scalar(@{$newcogprots{$a}})
	|| $a <=> $b } keys %newcogorgs;

# If in mode "c" or "b", re-map all COGs (so that multi's can look ahead and know which COG# to report)
if ($mode eq "c" || $mode eq "b")
{
	for (my $i=0; $i<scalar(@newcogs); $i++)
	{
		my $oldcog=$cogs[$cogmapping{$newcogs[$i]}];
		$newcogmapping{$oldcog}=$i;
	}
}

# Finally, print off the new set of modified COGs
for (my $i=0; $i<scalar(@newcogs); $i++)
{
	next if defined($skipthiscog{$newcogs[$i]}); # for mode "b", skip COGs that are strict subsets of another single, larger COG
	foreach my $prot (@{$newcogprots{$newcogs[$i]}})
	{
		my $line=$line{$prot}{$cogs[$cogmapping{$newcogs[$i]}]};
		my $newcogname=$prefix.sprintf("%0*d", $numpad, $i+$startingclusternumber);
		if (scalar(@{$protinnewcog{$prot}})>1)
		{
			my $tmpnewcogname="multi";
			my $count=0;
			foreach my $j (@{$protinnewcog{$prot}})
			{
				next if !defined($newcogmapping{$j}); # for mode "c", skip COGs that have no members not also seen in larger COGs
				$tmpnewcogname.=":$prefix".sprintf("%0*d", $numpad, $newcogmapping{$j}+$startingclusternumber);
				$count++;
			}
			if ($count>1)
			{
				$newcogname=$tmpnewcogname;
			}
		}
		$line=~s/,$prefix\d+,/,$newcogname,/;
		print OUT "$line\n";
	}
}

# Print singlets (order by organism name)
if ($filter<2)
{
	foreach $prot (sort {$prot2org{$a} cmp $prot2org{$b}} keys %singlets)
	{
		print OUT $singlets{$prot},"\n";
	}
}

close OUT;

exit(0);

