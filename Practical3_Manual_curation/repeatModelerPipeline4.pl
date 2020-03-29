#!/usr/bin/env/perl

# repeatModelerPipeline.pl
# written by Linnéa Smeds                          Mar 2013
# modified and updated by Alexander Suh		   Dec 2015
# =========================================================
#
# =========================================================
# usage:
# perl repeatModelerPipeline.pl <genome_fasta> \
# <BLAST_database> <RMDL_library>

use strict;
use warnings;
use List::Util qw(max min);	#Loading the math modules "max" and "min" 

# Input parameters 
my $ASSEMBLY = $ARGV[0];
my $blastDBpath = $ARGV[1];	#Need full path!
my $FASTA = $ARGV[2];		#Need full path!

# Make folders
system("mkdir blastn");
system("mkdir aligned");

# Other parameters "hard coded" in the script
# If some of these should be changed between the different genomes,
# it's better to have them as input parameters above
my $flank = 2000;
my $maxhitdist = 10000;	#Max length of the sum of close blasthits 
my $minfrac = 0.8;	#Minimum fraction of hits with the major orientation (+ or -)
my $digits = 5;		#set the number of digits in the outfiles
my $hits = 20;

my $blastdir = "blastn";
my $aligndir = "aligned";
my $tempBlastOut = "tempBlast.out";
my $tempMapNames = "tempMapNames.txt";
my $tempMafft = "mafft.out";

print "\nrepeatPipeline.pl started on ". localtime(). "\n";
print "-----------------------------------------------------\n";
 
################################################################################
## RUN THE STEPS

# BLAST 
&step1($FASTA, $blastDBpath, $tempBlastOut, $tempMapNames, $hits);

# FIND HITS FROM ASSEMBLYv AND AD ORIGINAL QUERY
&step2($FASTA, $blastdir, $maxhitdist, $minfrac, $ASSEMBLY);

#RUN THE ALIGNMENT WITH MAFFT
&step3($blastdir, $aligndir, $tempMafft);

# CLEANING UP
# by removing temporary files. 
# NOTE! if running stepwise, better comment away on this one since
# you need the temporary files throughout the whole pipeline 
&cleaning($tempBlastOut, $tempMapNames, $tempMafft);

################################################################################
### SUBROUTINES
# This is a subroutine, it's like a script within a script, that takes an input
# and can give an output.

## BLAST all files in repeat masked directory
sub step1 {

	my $FASTA = shift;
	my $blastDB = shift;
	my $tempBLastOut = shift;
	my $tempMapName = shift;
	my $hits = shift;
	
	print "STEP 1: Blast all sequences in $FASTA against the database $blastDB\n";

	# If there allready is a blast output file, remove it!
	system("rm -f $tempBlastOut");

	#Create a file that keeps track of the file name + fasta header
	open(OUT, ">$tempMapName");

#	print "DEBUG: Now we will run BLAST for all file in $dir\n";

	# Run blast for the full file
	system("blastn -db $blastDB -query $FASTA -outfmt 6 -evalue 10e-10 >>$tempBlastOut");

}

# FIND HITS FROM ASSEMBLY AND ADD ORIGINAL QUERY
sub step2 {
	
	my $FASTA = shift;
	my $blastdir = shift;
	my $maxhitdist = shift;
	my $frac = shift;
	my $ASSEMBLY = shift;
	
	# Blast directory
	my $dir = "$blastdir/";

	#How many bp to print on each line
	my $rowlength = 80;

	# Remove any already existing fasta files in the blast dir
	system("rm -f $dir/*.fa");
	system("rm -f $dir/*.txt*");
	
	print "STEP 2: Go through the blast hits:\n";

	#Make a hash of the assembly sequences
	my %seqHash = ();
	&saveGenomeToHash($ASSEMBLY, \%seqHash);

	# Go through the original infile to get the sequences and names
	open(IN, $FASTA); 
	while(<IN>) {
		if(/>/) {
			print "line is $_\n";
			my $name = $_;
			chomp($name);
			print "name is $name\n";
			$name =~ s/>//;
			my $seq="";

			my $next = <IN>;
			while ($next !~ m/^>/) {
				chomp($next);
				$seq .= $next;
				if(eof(IN)) {
					last;
				}	
				$next = <IN>;
			}
			seek(IN, -length($next), 1);

			#Outfile
			my $file = $name.".fa";
			my $out = "$dir/$file";
			open(OUT, ">$out");
			print OUT ">$name\n";
			my @seqParts = split(/(.{$rowlength})/, $seq);
			for my $seqs (@seqParts) {
				unless($seqs eq "") {
					print OUT $seqs."\n";
				}
			}
			
			#Make a temporary file for each query, and sort the output on (numerical) e value with 150 decimals
			my $tempRows = $name."_temp.out";
			system("awk '(\$1==\"$name\"){print}' $tempBlastOut |cut -f1,2,9,10,11 |perl -ne '\@a=split(/\\t/, \$_);print \$a[0].\"\\t\".\$a[1].\"\\t\".\$a[2].\"\\t\".\$a[3].\"\\t\";printf(\"%.150f\\n\", \$a[4])'| sort -k5n |head -n20 >$tempRows");
			

			#Make a list of unique combinations of column 1 and 2
			my $tempCombo = $name."_combTemp.out";
			system("cut -f1,2 $tempRows |uniq >$tempCombo");

			#
			open(COMB, $tempCombo);
			while(my $row = <COMB>) {	

				my($row1, $row2)=split(/\s+/,$row);
			
				my $tempMini = $name."_miniTemp.out";
				system("awk '(\$1==\"$row1\" && \$2==\"$row2\"){print}' $tempRows |sort -k3n >$tempMini");
			

				#Change the sequence header to the file name, and remove the latin name
				#Now when we open a new file while still having the other file ("IN")
				#open, we cannot use the default parameter "$_" for the line, but instead
				#create a parameter for it, here called "$line"
				open(TMP, $tempMini);
				my ($head, $fullName, $start, $stop);
				my ($minus, $plus) = (0,0);
				my $cnt = 0;
				while(my $line = <TMP>) {
					#First line - just save for comparison
					if($cnt==0) {
						($head, $fullName, $start, $stop) = split(/\s+/, $line);
						print "DEBUG: head is $head\n";
						if($stop>$start) {
							$plus++;
						}
						else {
							$minus++;
							my $temp=$start;
							$start=$stop;
							$stop=$temp;
						}
					}
					#All other lines
					else {
						my @thisline = split(/\s+/, $line);

						#If the line should be merged with the previous
						if(	min($thisline[2],$thisline[3])-max($stop,$start)<=$maxhitdist) {
								$start = min($start,$stop,$thisline[2],$thisline[3]);
								$stop = max($start,$stop,$thisline[2],$thisline[3]);
								if($thisline[3]>$thisline[2]) {
									$plus++;
								}
								else {
									$minus++;
								}
						}
						#if not, print previous result and save this line to compare with
						else {

							#Find the proper direction
							my $direction;
							if($plus>=$minus) {
								$direction="+";
							}
							else {
								$direction="-";
							}
							#Check that a majority of the hits has the same direction
							if(max($plus, $minus)/($plus+$minus) < $frac) {
								print STDERR "Ambiguous orientation of blast hits for ".$fullName.": ".$plus." +, ".$minus." -\n";
							} 

						
							#Extract sequences in subroutine and print directly to the right file
				#			print "DEBUG: out is $out\n";
							&extractFromFasta($fullName, $start-$flank, $stop+$flank, $direction, $out, "", \%seqHash);
						
					
							#Save this line:
							($head, $fullName, $start, $stop) = split(/\s+/, $line);
							if($stop>$start) {
								$plus=1;
							}
							else {
								$minus=1;
								my $temp=$start;
								$start=$stop;
								$stop=$temp;
							}
						}
					}
					$cnt++;
				}
				close(TMP);
				#The last line (or group of lines) has not been printed yet! Do it here:
	#			print "DEBUG: after while, head is $name\n";
			

				#Find the proper direction
				my $direction;
				if($plus>=$minus) {
					$direction="+";
				}
				else {
					$direction="-";
				}
				#Check that a majority of the hits has the same direction
				if(max($plus, $minus)/($plus+$minus) < $frac) {
					print STDERR "Ambiguous orientation of blast hits for ".$fullName.": ".$plus." +, ".$minus." -\n";
				} 

	#			print "DEBUG: Outside while, out is $out\n";
				&extractFromFasta($fullName, $start-$flank, $stop+$flank, $direction, $out, "", \%seqHash);

				#Remove temporary blast file
				#system("rm -f $tempRows");
			}
		}
	}
	close(IN);
}

# RUNNING ALIGNMENTS WITH MAFFT
sub step3 {

	my $blastdir = shift;
	my $aldir = shift;
	my $crapFile = shift;

	# Full path to blast directory
	my $olddir = "$blastdir/";
	my $newdir = "$aldir/";
	
	print "STEP 9: Aligning all the files in $olddir \n";
	print " printing alignment output to $newdir\n"; 

	# Make a handle to the directory (or abort if it doesn't exist)
	opendir(DIR, $olddir) or die "can't opendir $olddir: $!";

	# Go through all files in the directory
	while (defined(my $file = readdir(DIR))) {

		# Only look at fasta files
		if($file =~ m/\.fa/) {

			my $fullPathFile = $olddir."/".$file;
			my $output = $fullPathFile;
			$output =~ s/$olddir/$newdir/; 

			#Run mafft subprogram "einsi"
			system("einsi --thread 3 --adjustdirection $fullPathFile >$output 2>>$crapFile");	
	
		}
	}
	close(DIR);
}


# A little subroutine removing all the files it gets as input
sub cleaning {
	print "Cleaning up: removing temporary files...\n";
	system("rm -rf @_");			
}

#SAVE THE SEQUENCES TO A HASH
sub saveGenomeToHash {

	my $FILE = shift;

	# This part is just for making it possible to input a
	# zipped file (checks if the file ends with ".gz")
	if ($FILE =~ /\.gz$/) {
		open(FILE, "zcat $FILE |"); 
	}
	# If the file isn't compressed (default way of opening file handle)
	else {
		open(FILE, $FILE);
	}
	my ($head, $seq)=("","");
	my $cnt = 0;
	while(<FILE>) {		#The <> means that we read one line at the time

		#Found a header!
		if(/>/) {
			#First row - do nothing
			if($cnt==0) {
			}
			#When we find next sequence, we want to hash the previous one!
			else {
				#This is how to save stuff in a hash:
				${$_[0]}{$head}=$seq;
				$seq = "";
			}
			#Save the header temporarily to $head
			my @tabs = split(/\s+/,$_);	#split the line on spaces 
			$head = $tabs[0];	#takes the first part from the split
			$head =~ s/>//;	#removes the ">"
			#Since the repeat files seem to contain only the first part
			#of the name from the assembly file (for example scaffold10
			#if it in the fasta file was called "scaffold10	13.4") we
			#need to extract only the first part (to find something in
			#the hash you must use the EXACT string, down to the char):
			$cnt++;
		}
		#Found a sequence - save it temporarily to $seq!
		else {
			chomp($_);
			$seq.=$_;
		}			
	}	
	#hash the very last sequence 
	${$_[0]}{$head}=$seq;

	#Close the file handle
	close(FILE);	
}


# EXTRACT SEQUENCE FROM A FASTA FILE
# I copied this code from the "extractFromFasta" script
# and modified it a bit.
sub extractFromFasta {

	my $name = shift;		
	my $start = shift;		
	my $end = shift;
	my $direction = shift;	
	my $outfile = shift;
	my $newHeader=shift;	#If we want another name than the scaffold+start-stop
							#(otherwise set it to "" when calling the sub)
		
	# The following code (looks very strange, I know) is used when a hash is
	# sent to a subroutine as input, creating a local copy of it.
	my(%hashCopy) = %{$_[0]};	

#	 print "DEBUG: OUTFILE is $outfile\n";		#Uncomment this line for display

	my $noOfBases = $end-$start+1;
	my $rowlength = 80;		#How many bases per row


	# Just make sure that the sequence existed in the assembly (=in the hash)
	if(defined $hashCopy{$name}) {

		my $seq = $hashCopy{$name};
		
		# Opens a file handle to the output. One ">" means creating the file
		# from scratch (overwriting any existing file with this name), two (">>")
		# means append the file in the end. We use the latter since we add one
		# seqeuence at the time, opening the same file multiple times.
		open(OUT, ">>$outfile") or die "Couldn't create $outfile - check path!\n";
	
		my $len =  length($seq);

		# Here you can change how the header is printed to the new file
		if($newHeader eq "") {
			print OUT ">".$name."_".$start."-".$end."\n";
		}
		else {
			print OUT ">".$newHeader."\n";
		}

		
	
		# "substr" is a function which is pre-defined in perl
		my $substr = substr($seq, $start-1, $noOfBases);

		if($direction eq "-") {
			$substr = &revComp($substr);
			
		}

		my @seqParts = split(/(.{$rowlength})/, $substr);
		for my $seqs (@seqParts) {
			unless($seqs eq "") {
				print OUT $seqs."\n";
			}
		}
		close(OUT);
	}
	else {	
		print "ERROR: There is no sequence called $name\n";
	}
}

# REVERSE COMPLEMENT A SEQUENCE
sub revComp {
	
	my $DNAstring = shift;

	my $output = "";
	my @a = split(//, $DNAstring);

	for(@a) {
		$_ =~ tr/[A,T,C,G,R,Y,K,M,B,D,H,V,a,t,c,g,r,y,k,m,b,d,h,v]/[T,A,G,C,Y,R,M,K,V,H,D,B,t,a,g,c,y,r,m,k,v,h,d,b]/;
		$output = $_ . $output;
	}

	return $output;
}
	


