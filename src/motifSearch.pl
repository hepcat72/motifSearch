#!/usr/bin/perl

#USAGE: Run with no options to get usage or with --help for basic details

use CommandLineInterface;
use warnings;
use strict;

our $VERSION = '2.1';

setScriptInfo(CREATED => '8/23/2012',
              VERSION => $VERSION,
              AUTHOR  => 'Robert William Leach',
              CONTACT => 'rleach@princeton.edu',
              COMPANY => 'Princeton University',
              LICENSE => 'Copyright 2018',
              HELP    => << 'END_HELP'

Searches for perl regular expressions in sequence files.  This script outputs a tab-delimited file with coordinates.  Example:

    motifSearch.pl -s 'AUAU.{1,100}AUAU.{25,50}(?i:ggg)' -e -i rna.fa --verbose

END_HELP
	      ,
	      DETAILED_HELP => << 'END_DETAIL'

Any fasta-format file can be searched using a perl regular expression without concern for hard returns or white spaces.  If the sequence is DNA or RNA, the reverse complement will also be searched by default.

The perl regular expression must be a valid regular expression.  There may not be an error reported if a regular expression is mal-formed.  In developing your regular expression, try testing it with a test call in the following manner:

    motifSearch.pl -s 'your regular expression' -i 'a sequence that should match'

In other words, instead of submitting a file of sequences to search in, submit a short sequence you expect to match right on the command line.

END_DETAIL
);

setDefaults(HEADER        => 1,
	    ERRLIMIT      => 3,
	    COLLISIONMODE => 'error', #,merge,rename (when outfile conflict)
	    DEFRUNMODE    => 'usage');

my $patfileid =
  addInfileOption(GETOPTKEY   => 's|search-for=s',
		  PRIMARY     => 0,
		  SMRY_DESC   => "Perl regular expression to search for.",
		  DETAIL_DESC => << 'end_detail'

Perl regular expression to search for, either submitted as a string (or space-delimited strings) on the command line or in a file (see --help for file format).  Every sequence in the sequence file will be searched with this regular expression.  Multiple non-overlapping results will be reported per sequence.  Depending on the greediness of your regular expression, a hit could contain other hits.

end_detail
		  ,
		  FORMAT_DESC => << 'end_format'

You can either supply a perl regular expression (i.e. 'pattern') directly to this option on the command line (or a series of [not-escaped-] space-delimited regular expressions) or a text file of regular expressions.  The text file should be formatted such that each line is a separate perl regular expression to search for in all the sequences of the SEARCH-IN files.  No tabs, carriage returns, or newline characters are allowed in a pattern.  Output will indicate the pattern matched by the number of the line on which the pattern was found.  Optionally, a pattern may be preceded by an ID followed by a tab character.  An ID may not have any spaces or tabs in it.  The pattern may be followed by a tab and comments that will be ignored.  Commented lines (lines starting with '#') and lines with nothing on them will be ignored.

end_format
);

my $srcfileid =
  addInfileOption(GETOPTKEY   => 'i|search-in=s',
		  PRIMARY     => 1,
		  SMRY_DESC   => "Fasta sequence string or file.",
		  DETAIL_DESC => << 'end_detail'

Space-separated sequence file(s) or sequence strings to search in.  Note, this option expands glob characters ('*', '?', etc.), so multiple space-delimited glob file patterns can be submitted.  When standard input detected, used as a file name stub (See -o).  See --help for file format and advanced usage.

end_detail
		  ,
		  FORMAT_DESC => << 'end_format'

Fasta sequence file containing unique IDs on the deflines (the first string after '>' up to the first space or the end of the line).

end_format
);

my $exact_only = 0;
addOption(GETOPTKEY => 'e|exact-matches-only!',
	  GETOPTVAL => \$exact_only,
	  DEFAULT   => $exact_only,
	  SMRY_DESC => join('',('Treat ambiguous DNA characters like any ',
				'other character.  Use this option for non-',
				'DNA sequence searches.')));

my $case_sensitive = 0;
addOption(GETOPTKEY => 'c|case-sensitive!',
	  GETOPTVAL => \$case_sensitive,
	  DEFAULT   => $case_sensitive,
	  SMRY_DESC => join('',('Default behavior is case-insensitive.  ',
				'Supply this flag to match case.')));

my $strand = 'both';
addOption(GETOPTKEY => 'd|strand=s',
	  GETOPTVAL => \$strand,
	  DEFAULT   => $strand,
	  ACCEPTS   => ['forward','reverse','both'],
	  SMRY_DESC => join('',('If the sequence supplied to -i is DNA, ',
				'search the indicated strand(s).  If any ',
				'sequence has non-DNA characters, no attempt ',
				'will be made to search the reverse ',
				'complement.')));

my $min_size = 2;
addOption(GETOPTKEY => 'min-size=s',
	  GETOPTVAL => \$min_size,
	  DEFAULT   => $min_size,
	  SMRY_DESC => join('',('Ignore sequences supplied by -s that are ',
				'smaller than this size.')));

my $report_not_found = 0;
addOption(GETOPTKEY => 'report-not-found!',
	  GETOPTVAL => \$report_not_found,
	  DEFAULT   => $report_not_found,
	  SMRY_DESC => join('',('Report empty records when a search-for ',
				'sequence is not found in a search-in ',
				'sequence.')));

my $outfileid =
  addOutfileOption(GETOPTKEY   => 'o|outfile=s',
		   PRIMARY     => 1,
		   DEFAULT     => 'stdout',
		   PAIR_RELAT  => 'ONE',
		   SMRY_DESC   => "Tab delimited output file.",
		   FORMAT_DESC => << 'end_format'

Tab delimited file containing columns: file name/path where the pattern was found, source/parent sequence ID, start coordinate, stop coordinate, and strand (+/-).  Example:

    ../../hg19/hg19.fa      chr4    154156087       154156109       -
    ../../hg19/hg19.fa      chr11   134264352       134264374       -

end_format
);

processCommandLine();

if($strand =~ /^b/i)
  {$strand = 'both'}
elsif($strand =~ /^f/i)
  {$strand = 'forward'}
elsif($strand =~ /^r/i)
  {$strand = 'reverse'}
else
  {
    error("Strand (-d): [$strand] must be one of the following values: ",
	  "[forward,reverse,both].");
    usage(1);
    quit(1);
  }

my $seqs_found = 0;

while(nextFileCombo())
  {
    my $input_file  = getInfile($patfileid);
    my $parent_file = getInfile($srcfileid);
    my $outfile     = getOutfile($outfileid);

    debug("DOING motif [$input_file] and source [$parent_file].");

    my @target_records = ();
    my @source_records = ();

    if(-e $input_file || $input_file eq '-')
      {
	#this mapping is a kludge because this script started out as a copy of
	#getSeqCoords.pl, which expected a fasta file of sequences to search
	#for
	@target_records = map {$_->[0] = ">$_->[0]";$_}
	  grep {ref($_) eq 'ARRAY' && scalar(@$_) >= 2}
	    getPatternsFromFile($input_file);
      }
    else
      {push(@target_records,[">$input_file",$input_file])}

    foreach my $target_rec (@target_records)
      {
	my $defline = $target_rec->[0];
	my $pattern = $target_rec->[1];

	$defline =~ s/^>//;
	
	verboseOverMe("Looking for: [$pattern] in [$parent_file].");
      }

    my $rec_num = 0;

    #If there is a file that exists by this name, grab the fasta records
    if(-e $parent_file)
      {
	openIn(*PARENT,$parent_file) || next;
	unless(isDryRun())
	  {
	    $rec_num++;
	    @source_records = getNextFastaRec(*PARENT);
	    closeIn(*PARENT);
	  }
      }
    #Assume the submitted string is the sequence itself if it looks like DNA
    elsif(scalar(isDNA($parent_file)) || scalar(isRNA($parent_file)))
      {push(@source_records,[">$input_file",$parent_file])}
    else
      {
	error("[$parent_file] does not exist as a file and does not appear ",
	      "to be either a DNA or RNA sequence.  Skipping.");
	next;
      }

    next if(isDryRun());

    my($source_defline,$source_sequence,$source_id);
    my $parent_id_hash = {};
    foreach my $source_rec (@source_records)
      {
	if(scalar(@$source_rec) == 2)
	  {
	    ($source_defline,$source_sequence) = @$source_rec;

	    $source_defline =~ s/^\s*>\s*//;
	    if($source_defline =~ /^(\S+)/)
	      {
		$source_id = $1;
		if(exists($parent_id_hash->{$source_id}))
		  {
		    error("ID [$source_id] on defline [$source_defline] in ",
			  "parent file [$parent_file] already exists.  You ",
			  "appear to have multiple sequences in your ",
			  "parent file with the same ID.  IDs must ",
			  "be unique.  Only the first of each series ",
			  "of sequences with the same ID will be ",
			  "used.");
		  }
		else
		  {$parent_id_hash->{$source_id} = $source_defline}
	      }
	    else
	      {error("Unable to parse defline [$source_defline] in file ",
		     "[$parent_file].")}

	    #Make sure the sequence is DNA if not in exact_only mode
	    if(!$exact_only &&
	       !scalar(isDNA($source_sequence)) &&
	       !scalar(isRNA($source_sequence)))
	      {
		warning("Sequence(s) were encountered in parent file ",
			"[$parent_file] that appear to not be DNA when not ",
			"running in --exact-matches-only mode (e.g.: ",
			"'$source_defline').  Sequence treated as DNA will ",
			"not extract correct coordinates and you will get ",
			"lots of other errors.  Please check your file to ",
			"make sure it is DNA or turn on --exact-matches-",
			"only.");
	      }

	    if(length($source_sequence) < $min_size)
	      {
		debug("Parent sequence: [$source_defline] in parent file: ",
		      "[$parent_file] is below the minimum allowed size ",
		      "[$min_size].  Use --min-size to set the minumin ",
		      "allowed size.");
	      }
	  }
	else
	  {
	    error("Unable to parse parent file [$parent_file] line.\n",
		  "@$source_rec");
	    next;
	  }

	if($strand ne 'reverse')
	  {study($source_sequence)}

	my $rc_source_seq = '';
	if($strand ne 'forward' &&
	   (scalar(isDNA($source_sequence)) ||
	    scalar(isRNA($source_sequence))))
	  {
	    $rc_source_seq =
	      reverseComplement($source_sequence,
				!scalar(isDNA($source_sequence)) &&
				scalar(isRNA($source_sequence)));
	    study($rc_source_seq);
	  }

	my($target_defline,$target_sequence,$target_id);
	my $subseq_check = {};
	#Now go through the search sequences and search for them in the source
	#sequences
	foreach my $target_rec (@target_records)
	  {
	    if(scalar(@$target_rec) == 2)
	      {
		($target_defline,$target_sequence) = @$target_rec;

		if($target_sequence eq '' || $target_sequence =~ /\//)
		  {
		    error("Invalid motif: [$target_sequence] ",
			  (-e $input_file || $input_file eq '-' ?
			   "found in file [$input_file]" :
			   'provided to -s'),".  Skipping.");
		    next;
		  }

		$target_id = '';
		$target_defline =~ s/^\s*>\s*//;
		if($target_defline =~ /^(\S+)/)
		  {
		    $target_id = $1;
		    if(exists($subseq_check->{$target_id}))
		      {
			error("ID [$target_id] on defline [$target_defline] ",
			      "in your motif already exists.  ",
			      "You appear to have multiple sequences in your ",
			      "input file with the same ID.  IDs must ",
			      "be unique.  Only the first of each series ",
			      "of sequences with the same ID will be ",
			      "used.");
			next;
		      }
		    else
		      {$subseq_check->{$target_id} = $target_defline}
		  }
		else
		  {error("Unable to parse defline [$target_defline] in the ",
			 "motif.")}

		my $str = $target_sequence;
		$str =~ s/(?<=.{20}).+/.../;
		$str = "Looking for [$target_id: $str] in [$source_id] in " .
		  "file [$parent_file].";
		verboseOverMe("$str  Sequences found: $seqs_found");

		my $match_seq = $target_sequence;
		if(!$case_sensitive)
		  {$match_seq = "(?i:$match_seq)"}

		if(length($target_sequence) < $min_size)
		  {
		    warning("Skipping pattern: [$target_defline] because ",
			    "it is below the minimum allowed size: ",
			    "[$min_size].  Use --min-size to set the ",
			    "minimum allowed size.");
		    next;
		  }

		my $found_one = 0;
		my $parnum = 0;

		debug("Pattern: [$match_seq] Sequence: [$source_sequence]");

		my($start,$stop,$out);
		if($strand ne 'reverse')
		  {
		    while($source_sequence =~ /$match_seq/g)
		      {
			$start = length($`) + 1;
			$stop  = length($`) + length($&);
			$out   = reportCoords($target_id,
					      $source_id,
					      $start,
					      $stop,
					      '+',
					      $&);
			verboseOverMe("");
			print($out);
			$found_one = 1;
			$seqs_found++;
			verboseOverMe("$str  Sequences found: $seqs_found");
		      }
		  }
		$out = '';
		if($strand ne 'forward')
		  {
		    while($rc_source_seq =~ /$match_seq/g)
		      {
			$start = length($`) + length($&);
			$start = (length($source_sequence) + 1) - $start;
			$stop  = length($`) + 1;
			$stop  = (length($source_sequence) + 1) - $stop;
			$out   = reportCoords($target_id,
					      $source_id,
					      $start,
					      $stop,
					      '-',
					      $&);
			verboseOverMe("");
			print($out);
			$found_one = 1;
			$seqs_found++;
			verboseOverMe("$str  Sequences found: $seqs_found");
		      }
		  }
		if($report_not_found && !defined($out))
		  {
		    $out = reportCoords($target_id,$source_id,'','','','');
		    verboseOverMe("");
		    print($out);
		  }
	      }
	    else
	      {error("Unable to parse line in motif: [$input_file].\n",
		     "$_")}
	  }
      }
  }

verbose("Sequences found: $seqs_found");

##
## End Main
##






























sub getDNAPattern
  {
    my $seq    = $_[0];
    my $seqpat = '';
    my($char);
    while($seq =~ /(.)/g)
      {
	$char = $1;
	if($char =~ /N/i)
	  {$seqpat .= '.'}
	elsif($char =~ /B/i)
	  {$seqpat .= '[^A]'}
	elsif($char =~ /D/i)
	  {$seqpat .= '[^C]'}
	elsif($char =~ /H/i)
	  {$seqpat .= '[^G]'}
	elsif($char =~ /V/i)
	  {$seqpat .= '[^T]'}
	elsif($char =~ /R/i)
	  {$seqpat .= '[^CTY]'}
	elsif($char =~ /Y/i)
	  {$seqpat .= '[^AGR]'}
	elsif($char =~ /K/i)
	  {$seqpat .= '[^ACM]'}
	elsif($char =~ /M/i)
	  {$seqpat .= '[^TGK]'}
	elsif($char =~ /S/i)
	  {$seqpat .= '[^ATW]'}
	elsif($char =~ /W/i)
	  {$seqpat .= '[^GCS]'}
	elsif($char =~ /A/i)
	  {$seqpat .= '[ARMWDHVN]'}
	elsif($char =~ /G/i)
	  {$seqpat .= '[GRKSBDVN]'}
	elsif($char =~ /C/i)
	  {$seqpat .= '[CYMSBHVN]'}
	elsif($char =~ /T/i)
	  {$seqpat .= '[TYKWBDHN]'}
	else
	  {
	    error("Unrecognized character: [$char].");
	    $seqpat .= quotemeta($char);
	  }
      }

    return($seqpat)
  }

sub isDNA
  {
    my $sequence      = $_[0];
    my $no_ambiguous  = $_[1];   #if true, responds false if BDHVRYKMSWN presnt
    my $no_whitespace = $_[2];   #if true, responds false if white space found
    my $allow_mask    = $_[3];   #Allows DNA to have X's
    my $allow_align   = $_[4];   #Allows dashes (-)
    my $allow_ignore  = $_[5];   #Allows DNA to have dots (.)

    if(wantarray)
      {
	my @nonATGCwhite;
	@nonATGCwhite = ($sequence =~ /([^ATGC])/ig);
	@nonATGCwhite = grep {/\s/} @nonATGCwhite if(!$no_whitespace);
	@nonATGCwhite = grep {/[^BDHVRYKMSWN]/i} @nonATGCwhite
	  if(!$no_ambiguous);
	@nonATGCwhite = grep {/[^X]/i} @nonATGCwhite if($allow_mask);
	@nonATGCwhite = grep {/[^\-]/i} @nonATGCwhite if($allow_align);
	@nonATGCwhite = grep {/[^\.]/i} @nonATGCwhite if($allow_ignore);
	return(@nonATGCwhite);
      }

    if($sequence =~ /[^ATGCBDHVRYKMSWNX\-\.\s]/i)
      {return(0)}
    if($no_whitespace && $sequence =~ /\s/)  #not ambiguous
      {return(0)}
    if(!$allow_ignore && $sequence =~ /\./)  #counts as ambiguous
      {return(0)}
    if(!$allow_align  && $sequence =~ /-/)   #not ambiguous
      {return(0)}
    if(!$allow_mask   && $sequence =~ /X/i)  #counts as ambiguous
      {return(0)}
    if($no_ambiguous  && $sequence =~ /[BDHVRYMKSWN]/i)
      {return(0)}

    #At this point, I know there's no non-DNA characters, so all I need to do
    #is find out the ratio of ambiguous characters (including X's and no .'s)
    my $num_ambig = scalar($sequence =~ /([BDHVRYKMSWNX\.])/ig);
    if($num_ambig)
      {
	my $tmp = $sequence;
	$sequence =~ s/[\-\s]//g; #Remove gaps and white space for length calc.
	my $len = length($tmp);
	return(-1) if($len == 0);
	return(($len - $num_ambig) / $len);
      }

    return(1);
  }

sub isRNA
  {
    my $sequence      = $_[0];
    my $no_ambiguous  = $_[1];   #if true, responds false if BDHVRYKMSWN presnt
    my $no_whitespace = $_[2];   #if true, responds false if white space found
    my $allow_mask    = $_[3];   #Allows DNA to have X's
    my $allow_align   = $_[4];   #Allows dashes (-)
    my $allow_ignore  = $_[5];   #Allows DNA to have dots (.)

    if(wantarray)
      {
	my @nonAUGCwhite;
	@nonAUGCwhite = ($sequence =~ /([^AUGC])/ig);
	@nonAUGCwhite = grep {/\s/} @nonAUGCwhite if(!$no_whitespace);
	@nonAUGCwhite = grep {/[^BDHVRYKMSWN]/i} @nonAUGCwhite
	  if(!$no_ambiguous);
	@nonAUGCwhite = grep {/[^X]/i}  @nonAUGCwhite if($allow_mask);
	@nonAUGCwhite = grep {/[^\-]/i} @nonAUGCwhite if($allow_align);
	@nonAUGCwhite = grep {/[^\.]/i} @nonAUGCwhite if($allow_ignore);
	return(@nonAUGCwhite);
      }

    if($sequence =~ /[^AUGCBDHVRYKMSWNX\-\.\s]/i)
      {return(0)}
    if($no_whitespace && $sequence =~ /\s/)  #not ambiguous
      {return(0)}
    if(!$allow_ignore && $sequence =~ /\./)  #counts as ambiguous
      {return(0)}
    if(!$allow_align  && $sequence =~ /-/)   #not ambiguous
      {return(0)}
    if(!$allow_mask   && $sequence =~ /X/i)  #counts as ambiguous
      {return(0)}
    if($no_ambiguous  && $sequence =~ /[BDHVRYMKSWN]/i)
      {return(0)}

    #At this point, I know there's no non-DNA characters, so all I need to do
    #is find out the ratio of ambiguous characters (including X's and no .'s)
    my $num_ambig = scalar($sequence =~ /([BDHVRYKMSWNX\.])/ig);
    if($num_ambig)
      {
	my $tmp = $sequence;
	$sequence =~ s/[\-\s]//g; #Remove gaps and white space for length calc.
	my $len = length($tmp);
	return(-1) if($len == 0);
	return(($len - $num_ambig) / $len);
      }

    return(1);
  }

#Copied from fetch_cog.pl.pl on 8/6/2008 -Rob
sub getNextFastaRec
  {
    my $handle    = $_[0];      #File handle or file name
    my $no_format = $_[1];

    if(exists($main::{FASTABUFFER}) && exists($main::{FASTABUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTABUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTABUFFER}->{$handle}});
		@{$main::{FASTABUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTABUFFER}->{$handle}}));
	  }
	elsif(eof($handle))
	  {return(undef)}
      }

    my $parent_id_check = {};
    my $first_loop      = 0;
    my $line_num        = 0;
    my $line            = '';
    my $defline         = '';
    my $verbose_freq    = 1000;
    my($seq);

    #For each line in the current input file
    while(getLine($handle))
      {
	$line_num++;
	$line = $_;

	verboseOverMe("Reading line [$line_num].")
	  unless($line_num % $verbose_freq);

	next if($line !~ /\S/ || $line =~ /^\s*#/);
	if($line =~ />/)
	  {
	    if($defline)
	      {
		my $solidseq =
		  ($no_format ? $seq :
		   formatSequence($seq));
		chomp($solidseq);
		chomp($defline);

		push(@{$main::{FASTABUFFER}->{$handle}},[$defline,$solidseq]);
	      }
	    $defline = $line;

	    my $tmp_id = $defline;
	    $tmp_id =~ s/^\s*>\s*//;
	    $tmp_id =~ s/\s.*//;
	    if($tmp_id eq '')
	      {warning("No Defline ID on line: [$line_num] of current file.  ",
		       " Universal coordinates will be used if some were ",
		       "supplied either via command line arguments of via ",
		       "coordinate file with no parent sequence ID.")}
	    elsif(exists($parent_id_check->{$tmp_id}))
	      {
		error("Two sequences found with the same ID on the ",
		      "defline: [$tmp_id] in current fasta file.  The same ",
		      "pairs of coordinates will be used for each sequence.");
	      }

	    undef($seq);
	  }
	elsif($line =~ /^([^\t]+?) *\t\s*(.*)/)
	  {
	    $defline = $1;
	    $seq     = $2;

	    my $solidseq =
	      ($no_format ? $seq :
	       formatSequence($seq));
	    chomp($solidseq);
	    chomp($defline);

	    push(@{$main::{FASTABUFFER}->{$handle}},[$defline,$solidseq]);

	    undef($seq);
	  }
	else
	  {$seq .= $line}
      }

    #Handle the last sequence (if there were any sequences)
    if(defined($seq))
      {
	my $solidseq =
	  ($no_format ? $seq :
	   formatSequence($seq));
	chomp($solidseq);
	chomp($defline);

	push(@{$main::{FASTABUFFER}->{$handle}},[$defline,$solidseq]);
      }

    #Return the first sequence (if sequence was parsed)
    if(exists($main::{FASTABUFFER}) && exists($main::{FASTABUFFER}->{$handle}))
      {
	if(scalar(@{$main::{FASTABUFFER}->{$handle}}) > 0)
	  {
	    if(wantarray)
	      {
		my @array = (@{$main::{FASTABUFFER}->{$handle}});
		@{$main::{FASTABUFFER}->{$handle}} = ();
		return(@array);
	      }
	    return(shift(@{$main::{FASTABUFFER}->{$handle}}));
	  }
	else
	  {return(undef)}
      }
    else
      {return(undef)}
  }

#Copied from uniqueSeq.pl on11/3/2009 so as to be independent -Rob
sub formatSequence
  {
    #1. Read in the parameters.
    my $sequence          = $_[0];
    my $chars_per_line    = $_[1];
    my $coords_left_flag  = $_[2];
    my $coords_right_flag = $_[3];
    my $start_coord       = $_[4];
    my $coords_asc_flag   = $_[5];
    my $coord_upr_bound   = $_[6];
    my $uppercase_flag    = $_[7];
    my $print_flag        = $_[8];
    my $nucleotide_flag   = $_[9];

    my($formatted_sequence,
       $sub_string,
       $sub_sequence,
       $coord,
       $max_num_coord_digits,
       $line_size_left,
       $lead_spaces,
       $line);
    my $coord_separator = '  ';
    my $tmp_sequence = $sequence;
    $tmp_sequence =~ s/\s+//g;
    $tmp_sequence =~ s/<[^>]*>//g;
    my $seq_len = length($tmp_sequence);

    #2. Error check the parameters and set default values if unsupplied.
    my $default_chars_per_line    = ''; #Infinity
    my $default_coords_left_flag  = 0;
    my $default_coords_right_flag = 0;
    my $default_start_coord       = (!defined($coords_asc_flag) ||
				     $coords_asc_flag ? 1 : $seq_len);
    my $default_coords_asc_flag   = 1;
    my $default_coord_upr_bound   = undef();  #infinity (going past 1 produces
    my $default_uppercase_flag    = undef();  #          negative numbers)
    my $default_print_flag        = 0;

    if(!defined($chars_per_line) || $chars_per_line !~ /^\d+$/)
      {
        if(defined($chars_per_line) &&
	   $chars_per_line !~ /^\d+$/ && $chars_per_line =~ /./)
	  {print("WARNING:seq-lib.pl:formatSequence: Invalid ",
	         "chars_per_line: [$chars_per_line] - using default: ",
		 "[$default_chars_per_line]<BR>\n")}
        #end if(chars_per_line !~ /^\d+$/)
	$chars_per_line = $default_chars_per_line;
      }
    elsif(!$chars_per_line)
      {$chars_per_line = ''}
    #end if(!defined($chars_per_line) || $chars_per_line !~ /^\d+$/)
    if(!defined($coords_left_flag))
      {$coords_left_flag = $default_coords_left_flag}
    #end if(!defined(coords_left_flag))
    if(!defined($coords_right_flag))
      {$coords_right_flag = $default_coords_right_flag}
    #end if(!defined(coords_right_flag))
    if(!defined($start_coord) || $start_coord !~ /^\-?\d+$/)
      {
        if(defined($start_coord) &&
           ($coords_left_flag || $coords_right_flag))
          {print("WARNING:formatSequence.pl:formatSequence: Invalid ",
                 "start_coord: [$start_coord] - using default: ",
                 "[$default_start_coord]\n")}
        #end if($start_coord !~ /^\d+$/)
        $start_coord = $default_start_coord;
      }
    #end if(!defined($start_coord) || $start_coord !~ /^\d+$/)
    if(!defined($coords_asc_flag))
      {$coords_asc_flag = $default_coords_asc_flag}
    #end if(!defined(coords_right_flag))
    if(defined($coord_upr_bound) && $coord_upr_bound !~ /^\d+$/)
      {undef($coord_upr_bound)}
    if(!defined($print_flag))
      {$print_flag = $default_print_flag}
    #end if(!defined($print_flag))

    if(defined($coord_upr_bound) && $start_coord < 1)
      {$start_coord = $coord_upr_bound + $start_coord}
    elsif($start_coord < 1)
      {$start_coord--}
    elsif(defined($coord_upr_bound) && $start_coord > $coord_upr_bound)
      {$start_coord -= $coord_upr_bound}

    #3. Initialize the variables used for formatting.  (See the DATASTRUCTURES
    #   section.)
    if($coords_asc_flag)
      {
        if(defined($coord_upr_bound) &&
           ($seq_len + $start_coord) > $coord_upr_bound)
          {$max_num_coord_digits = length($coord_upr_bound)}
        else
          {$max_num_coord_digits = length($seq_len + $start_coord - 1)}

        $coord = $start_coord - 1;
      }
    else
      {
        if(defined($coord_upr_bound) && ($start_coord - $seq_len + 1) < 1)
          {$max_num_coord_digits = length($coord_upr_bound)}
        elsif(!defined($coord_upr_bound) &&
              length($start_coord - $seq_len - 1) > length($start_coord))
          {$max_num_coord_digits = length($start_coord - $seq_len - 1)}
        else
          {$max_num_coord_digits = length($start_coord)}

        $coord = $start_coord + 1;
      }
    $line_size_left = $chars_per_line;
    $lead_spaces    = $max_num_coord_digits - length($start_coord);

    #5. Add the first coordinate with spacing if coords_left_flag is true.
    $line = ' ' x $lead_spaces . $start_coord . $coord_separator
      if($coords_left_flag);

    #6. Foreach sub_string in the sequence where sub_string is either a
    #   sub_sequence or an HTML tag.
    foreach $sub_string (split(/(?=<)|(?<=>)/,$sequence))
      {
        #6.1 If the substring is an HTML tag
        if($sub_string =~ /^</)
          #6.1.1 Add it to the current line of the formatted_sequence
          {$line .= $sub_string}
        #end if(sub_string =~ /^</)
        #6.2 Else
        else
          {
            $sub_string =~ s/\s+//g;

	    if($nucleotide_flag)
	      {
		my(@errors);
		(@errors) = ($sub_string =~ /([^ATGCBDHVRYKMSWNX])/ig);
		$sub_string =~ s/([^ATGCBDHVRYKMSWNX])//ig;
		if(scalar(@errors))
		  {print STDERR ("WARNING:formatSequence.pl:formatSequence:",
				 scalar(@errors),
				 " bad nucleotide characters were ",
				 "filtered out of your sequence: [",
				 join('',@errors),
				 "].\n")}
	      }

            #6.2.1 If the sequence is to be uppercased
            if(defined($uppercase_flag) && $uppercase_flag)
              #6.2.1.1 Uppercase the sub-string
              {$sub_string = uc($sub_string)}
            #end if(defined($uppercase_flag) && $uppercase_flag)
            #6.2.2 Else if the sequence is to be lowercased
            elsif(defined($uppercase_flag) && !$uppercase_flag)
              #6.2.2.1 Lowercase the sub-string
              {$sub_string = lc($sub_string)}
            #end elsif(defined($uppercase_flag) && !$uppercase_flag)

            #6.2.3 While we can grab enough sequence to fill the rest of a line
            while($sub_string =~ /(.{1,$line_size_left})/g)
              {
                $sub_sequence = $1;
                #6.2.3.1 Add the grabbed sequence to the current line of the
                #        formatted sequence
                $line .= $sub_sequence;
                #6.2.3.2 Increment the current coord by the amount of sequence
                #        grabbed
                my $prev_coord = $coord;
                if($coords_asc_flag)
                  {
                    $coord += length($sub_sequence);
                    if(defined($coord_upr_bound)      &&
                       $prev_coord <= $coord_upr_bound &&
                       $coord > $coord_upr_bound)
                      {$coord -= $coord_upr_bound}
                  }
                else
                  {
                    $coord -= length($sub_sequence);
                    if(defined($coord_upr_bound) &&
                       $prev_coord >= 1 && $coord < 1)
                      {$coord = $coord_upr_bound + $coord - 1}
                    elsif($prev_coord >= 1 && $coord < 1)
                      {$coord--}
                  }
                #6.2.3.3 If the length of the current sequence grabbed
                #        completes a line
                if($line_size_left eq '' ||
		   length($sub_sequence) == $line_size_left)
                  {
                    $lead_spaces = $max_num_coord_digits - length($coord);
                    #6.2.3.3.1 Conditionally add coordinates based on the
                    #          coords flags
                    $line .= $coord_separator . ' ' x $lead_spaces . $coord
                      if($coords_right_flag);

                    #6.2.3.3.2 Add a hard return to the current line of the
                    #          formatted sequence
                    $line .= "\n";

                    #6.2.3.3.3 Add the current line to the formatted_sequence
                    $formatted_sequence .= $line;
                    #6.2.3.3.4 Print the current line if the print_flag is true
                    print $line if($print_flag);

                    #6.2.3.3.5 Start the next line
                    $lead_spaces = $max_num_coord_digits - length($coord+1);
                    $line = '';
                    $line = ' ' x $lead_spaces
                          . ($coords_asc_flag ? ($coord+1) : ($coord-1))
                          . $coord_separator
                      if($coords_left_flag);

                    #6.2.3.3.6 Reset the line_size_left (length of remaining
                    #          sequence per line) to chars_per_line
                    $line_size_left = $chars_per_line;
                  }
                #end if(length($sub_sequence) == $line_size_left)
                #6.2.3.4 Else
                else
                  #6.2.3.4.1 Decrement line_size_left (length of remaining
                  #          sequence per line) by the amount of sequence
                  #          grabbed
                  {$line_size_left -= length($sub_sequence)}
                #end 6.2.3.4 Else
              }
            #end while($sub_string =~ /(.{1,$line_size_left})/g)
          }
        #end 6.2 Else
      }
    #end foreach $sub_string (split(/(?=<)|(?<=>)/,$sequence))
    #7. Add the last coodinate with enough leadin white-space to be lined up
    #   with the rest coordinates if the coords_right_flag is true
    $lead_spaces = $max_num_coord_digits - length($coord);
    $line .= ' ' x $line_size_left . $coord_separator . ' ' x $lead_spaces
          . $coord
      if($coords_right_flag && $line_size_left != $chars_per_line);
    $line =~ s/^\s*\d+$coord_separator\s*$// if($coords_left_flag);

    #8. Add the ending PRE tag to the last line of the formatted sequence
    $line =~ s/\n+$/\n/s;

    #9. Add the last line to the formatted_sequence
    $formatted_sequence .= $line;
    #10. Print the last line if the print_flag is true
    print "$line\n" if($print_flag);

    if($coord < 1 && ($coords_left_flag || $coords_right_flag))
      {print("WARNING: The sequence straddles the origin.  Coordinates are ",
             "inaccurate.")}

    #11. Return the formatted_sequence
    return $formatted_sequence;
  }

#copied from uniqueSeq on 11/2/2009
sub reverseComplement
  {
    #1. Read in the sequence parameter.
    my $sequence = $_[0];
    my $rna      = defined($_[1]) ? $_[1] : 0;
    my @errors;
    if($rna)
      {if(@errors =
	  ($sequence =~ /([^AUGCBVDHRYKMSWN\.\-\s\r])/isg))
	 {error("Bad character(s) found: ['" . join("','",@errors) .
		"'] in [$sequence]")}}
    else
      {if(@errors =
	  ($sequence =~ /([^ATGCBVDHRYKMSWN\.\-\s\r])/isg))
	 {error("Bad character(s) found: ['" . join("','",@errors) .
		"'] in [$sequence]")}}
    #end if(@errors = ($sequence =~ /([^ATGCBVDHRYKM\s\r])/isg))
    #2. Transcribe the new_sequence.
    if($rna)
      {$sequence =~ tr/AUGCBVDHRYKMaugcbvdhrykm/UACGVBHDYRMKuacgvbhdyrmk/}
    else
      {$sequence =~ tr/ATGCBVDHRYKMatgcbvdhrykm/TACGVBHDYRMKtacgvbhdyrmk/}
    #3. Reverse the new_sequence.
    $sequence = reverse($sequence);
    return $sequence;
  }

sub reportCoords
  {
    my $target_id = $_[0];
    my $source_id = $_[1];
    my $start     = $_[2];
    my $stop      = $_[3];
    my $strand    = $_[4];
    my $sequence  = $_[5];

    return("$target_id\t$source_id\t$start\t$stop\t$strand\n");
  }

sub getPatternFromFile
  {
    my $input_file = $_[0];
    my $pattern    = '';

    openIn(*INPUT,$input_file) || return($pattern);

    my $line_num     = 0;
    my $verbose_freq = 100;

    #For each line in the current input file
    while(getLine(*INPUT))
      {
	$line_num++;
	verboseOverMe("[$input_file] Reading line: [$line_num].")
	  unless($line_num % $verbose_freq);

	chomp;

	next if(/#/ || /^$/);

	if($pattern ne '' && $pattern !~ /(?<!\\)\|$/)
	  {$pattern .= '|'}

	$pattern .= $_;
      }

    closeIn(*INPUT);

    return($pattern);
  }

sub getPatternsFromFile
  {
    my $input_file = $_[0];
    my $patterns   = [];

    openIn(*INPUT,$input_file) || return(wantarray ? @$patterns : $patterns);

    my $line_num     = 0;
    my $verbose_freq = 100;
    my $cnt          = 1;

    #For each line in the current input file
    while(getLine(*INPUT))
      {
	$line_num++;
	verboseOverMe("[$input_file] Reading line: [$line_num].")
	  unless($line_num % $verbose_freq);

	chomp;

	next if(/#/ || /^$/);

	my($id,$pattern);
	my(@data) = split(/\t/,$_,-1);

	if(scalar(@data) == 1)
	  {
	    $id      = $cnt;
	    $pattern = $data[0];
	  }
	elsif(scalar(@data))
	  {
	    $id      = $data[0];
	    $pattern = $data[1];
	  }
	else
	  {next}

	#This is a kludge because I copied code from getSeqCoords which assumes
	#that the sequences we're searching for came from a fasta file and not
	#this tab-delimited file
	$id =~ s/\s+//g;

	$cnt++;

	if(!defined($pattern) || $pattern eq '')
	  {
	    error("Invalid pattern: [$pattern] on line [$line_num] of file ",
		  "[$input_file].  Skipping.");
	    next;
	  }

	push(@$patterns,[$id,$pattern]);
      }

    closeIn(*INPUT);

    return(wantarray ? @$patterns : $patterns);
  }
