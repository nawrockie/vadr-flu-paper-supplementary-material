#!/usr/bin/env perl
# 
# table-errors.pl: gather data and create the latex file for the table mapping FLAN and VADR errors
#                        
# EPN, Tue Aug 15 13:42:25 2023
# 
#
use strict;
use warnings;
use Getopt::Long;

my %flan_words_H = ();
$flan_words_H{"CDS_HAS_FRAMESHIFT,PEPTIDE_FRAMESHIFT"} = "The coding region of /Mature peptide X has a frameshift";
$flan_words_H{"CDS_HAS_STOP_CODON"} = "The coding region of X has stop codon inside exon";
$flan_words_H{"CODING_CAPACITY"} = "This sequence does not have coding capacity"; 
#$flan_words_H{"CONTAMINATION_DOWNSTREAM"} = "contains extra X nts downstream the consensus 3' end sequence of influenza virsuses. Check for possible vector/linker contamination.";
#$flan_words_H{"CONTAMINATION_UPSTREAM"} = "contains extra X nts upstream the consensus 5' end sequence of influenza viruses. Check for possible vector/linker contamination.";
$flan_words_H{"CONTAMINATION_DOWNSTREAM"} = "contains extra X nts downstream the consensus 3' end sequence";
$flan_words_H{"CONTAMINATION_UPSTREAM"} = "contains extra X nts upstream the consensus 5' end sequence";
$flan_words_H{"MUTATION_AT_END"} = "Probable mutation at End of protX";
$flan_words_H{"MUTATION_AT_START"} = "Probable mutation at Start of protX";
$flan_words_H{"NO_BLAST_HITS_FOUND"} = "No blast hits found";
$flan_words_H{"PEPTIDE_OVERLAP,PEPTIDE_SEPARATED"} = "Mature peptides (X) and (Y) have overlap/are separated";
$flan_words_H{"REVCOMPLEM"} = "The input sequence is the reverse complementary strand of the coding sequence";
$flan_words_H{"SPLICE_SITE_NOT_FOUND"} = "Expected splice site consensus sequence not found for protein X";
$flan_words_H{"WRONG_SEGMENT"} = "Wrong exon number X for segment Y protein Z";
    
my $in_data  = ""; # name of input data file

my $usage = "perl table-errors.pl\n\t<parse-and-compare-flan-and-vadr-output.pl script>\n\t<flan error/warning list file>\n\t<vadr alert list file>\n\n\tOPTIONS:\n";
$usage .= "\t-v: verbose mode; output info on each error, whether it is unique or consistent\n\n";
#$usage .= "\t-warn: FLAN list file is warnings not errors and VADR list file is nonfatal alerts\n\n";

my $do_verbose = 0;
&GetOptions( "v"   => \$do_verbose);

if(scalar(@ARGV) != 3) { die $usage; }
my ($parse_script, $flan_list, $vadr_list) = @ARGV;

open(IN, $parse_script) || die "ERROR unable to open $parse_script for reading";

# parse parse-and-compare-flan-and-vadr-output.pl script
my ($type, $flan_error, $vadr_error) = (undef, undef, undef);
my %flan2vadr_H = ();
my %tbl_HH = ();
my %vadr_code2flan_caps_H = ();
my %single_flan_caps2flan_caps_H = ();
while(my $line = <IN>) { 
  chomp $line;
  if($line =~ m/^\s+\#\!\#/) { 
    # print $line . "\n";
    # type 1: mapping
    #!# WARNING vadr:MUTATION_AT_END:mutendcd,mutendns,mutendex == flan:MUTATION_AT_END:Probable mutation at End of prot<d>
    #!# ERROR vadr:PEPTIDE_ADJACENCY_PROBLEM:pepadjcy == flan:PEPTIDE_SEPARATED:Mature peptides (<name>) and <name> are separated
    # type 2: ignored flan 
    #!# IGNORED WARNING flan:NONE:Ambiguity nucleotide(s) found:/) { 
    #!# IGNORED ERROR flan:INSERTION_OF_NT:Insertion of <d+> nt around: <d> - <d>
    ($type, $flan_error, $vadr_error) = (undef, undef, undef);
    if($line =~ m/^\s+\#\!\#\s+(\S+)\s+vadr\:(\S+)\s+\=\=\s+flan\:(.+)$/) { 
      ($type, $vadr_error, $flan_error) = ($1, $2, $3);
    }
    elsif($line =~ m/^\s+\#\!\#\s+(IGNORED)\s+(\S+)\s+flan\:(.+)$/) { 
      ($vadr_error, $type, $flan_error) = ($1, $2, $3);
    }
    elsif($line =~ m/^\s+\#\!\#\s+if \-\-noextra/) { 
      # ignore these lines since we won't use --noextra 
      ($flan_error, $vadr_error) = (undef, undef);
    }
    else { 
      die "ERROR unable to parse map line: $line";
    }
    if((defined $type) && ($type eq "ERROR") && ($vadr_error ne "IGNORED")) { 
      if((defined $flan_error) && (defined $vadr_error)) { 
        # skip WARNINGs
        #print("type: $type\n");
        my ($flan_caps, $flan_explanation) = split(":", $flan_error);
        my ($vadr_caps, $vadr_code) = split(":", $vadr_error);
        if(! defined $vadr_code) { 
          if($vadr_error eq "NONENONE") { 
            $vadr_code = "nonenone";
          }
          else { 
            die "ERROR vadr_code undefined and not NONENONE at line: $line";
          }
        }
        #printf("mapping\n\tflan_error: $flan_error\n\t$flan_caps:$flan_explanation\n\tvadr_error:$vadr_error\n\t$vadr_caps:$vadr_code\n");
        # special case deal with PEPTIDE_HAS_FRAMESHIFT == CDS_HAS_FRAMESHIFT (which both map to the POSSIBLE_FRAMESHIFT vadr errors)
        # and PEPTIDE_SEPARATED == PEPTIDE_OVERLAP (which both map to the PEPTIDE_ADJACENCY vadr errors)
        if(($flan_caps ne "PEPTIDE_FRAMESHIFT") && ($flan_caps ne "PEPTIDE_SEPARATED")) { 
          if($flan_caps eq "CDS_HAS_FRAMESHIFT") {
            $flan_caps .= ",PEPTIDE_FRAMESHIFT";
          }
          if($flan_caps eq "PEPTIDE_OVERLAP") { 
            $flan_caps .= ",PEPTIDE_SEPARATED";
          }
          if(defined $tbl_HH{$flan_caps}) { 
            die "ERROR two entries for flan caps $flan_caps";
          }
          %{$tbl_HH{$flan_caps}} = ();
          $tbl_HH{$flan_caps}{"vadr_caps"} = $vadr_caps;
          $tbl_HH{$flan_caps}{"vadr_code"} = $vadr_code;
          $tbl_HH{$flan_caps}{"flan_exp"}  = $flan_explanation;
          
          if(defined $vadr_code2flan_caps_H{$vadr_code}) { 
            $vadr_code2flan_caps_H{$vadr_code} .= "," . $flan_caps;
            #print "ERROR two flan_caps values for vadr_code: $vadr_code ($flan_caps and " . $vadr_code2flan_caps_H{$vadr_code}. ")\n"; 
          }
          else { 
            $vadr_code2flan_caps_H{$vadr_code} = $flan_caps; 
          }
          #printf("\tset vadr_code2flan_caps_H{$vadr_code} to $flan_caps\n");
          
          $single_flan_caps2flan_caps_H{$flan_caps} = $vadr_code2flan_caps_H{$vadr_code};
          ## unique part should be the FLAN message
          #if(defined $flan2vadr_H{$flan_error}) { 
          #  if($flan2vadr_H{$flan_error} ne $vadr_error) { 
          #    die "ERROR read same FLAN error twice mapped to different VADR alerts: $flan_error";
          #  }
          #}
          #$flan2vadr_H{$flan_error} = $vadr_error;
          #printf("mapped $flan_error to $vadr_error\n");
        } # end of 'if($flan_caps ne "PEPTIDE_HAS_FRAMESHIFT") {'

      }
      else { 
        die "ERROR undefined flan_error or vadr_error for line: $line\n";
      }
      #elsif(defined $flan_error) { 
      #  my ($flan_caps, $flan_explanation) = split(":", $flan_error);
      #  printf("flan-unique $flan_caps:$flan_explanation\n");
      #}
      #elsif(defined $vadr_error) { 
      #  my ($vadr_caps, $vadr_code) = split(":", $vadr_error);
      #    printf("vadr-unique $vadr_caps:$vadr_code\n");
      #}
    }
  }
}
close(IN);

my %flan_perseq_HHH; # key1: $seqname, key2: "sequence" or ftr name, key3: flan_caps error value (e.g. CONTAMINATION_UPSTREAM)
my %flan_percaps_H;  # key1: flan_caps value, value number of instances of flan_caps 
# read in flan list file, and add flan caps based on vadr_code2flan_caps_H
open(IN, $flan_list) || die "ERROR unable to open $flan_list for reading";
while(my $line = <IN>) { 
  chomp $line;
  #NC_006307 sequence extrant3  ERROR: Sequence (gi|1250175388|ref|NC_006307.2|) contains extra 19 nts downstream the consensus 3' end sequence of influenza viruses. Check for possible vector/linker contamination.
  #M15088 CDS.NEP fsthicft,fsthicfi,fstukcft,fstukcfi  ERROR: The coding region of NEP has frameshift
  my @el_A = split(/\s+/, $line);
  my ($seqname, $seq_or_ftr, $code, $desc) = ($el_A[0], $el_A[1], $el_A[2], $el_A[3]);
  #printf("READ FLAN code $code\n");
  for(my $i = 4; $i < scalar(@el_A); $i++) { 
    $desc .= " " . $el_A[$i];
  }
  if(! defined $vadr_code2flan_caps_H{$code}) { 
    die "ERROR parsing FLAN list, unable to lookup $code in vadr_code2flan_caps_H\nline:$line\n";
  }
  my $flan_caps = $vadr_code2flan_caps_H{$code};
  if(($flan_caps !~ m/FRAMESHIFT/) || ($seq_or_ftr !~ m/mat\_peptide/)) { 
    # ignore frameshifts in mat_peptides, they are always accompanied by a frameshift in the CDS (I checked all train/test seqs)
    if(! defined $flan_perseq_HHH{$seqname}) { 
      %{$flan_perseq_HHH{$seqname}} = ();
    }
    if(! defined $flan_perseq_HHH{$seqname}{$seq_or_ftr}) { 
      %{$flan_perseq_HHH{$seqname}{$seq_or_ftr}} = ();
    }
    if(! defined $flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps}) { 
      $flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} = 1;
      if(! defined $flan_percaps_H{$flan_caps}) { 
        $flan_percaps_H{$flan_caps} = 0;
      }
      $flan_percaps_H{$flan_caps}++;
      #printf("updated flan_percaps_H{$flan_caps} to $flan_percaps_H{$flan_caps}\n");
    }
    elsif($flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} == 1) { 
      # do nothing; we don't count multiple errors for the same seq/ftr
    }
    else { 
      die "ERROR unexpected non-1 value for flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps}: " . $flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} . "\n";
    }
  }
}
close(IN);

my %vadr_perseq_HHH; # key1: $seqname, key2: "sequence" or ftr name, key3: flan_caps error value (e.g. CONTAMINATION_UPSTREAM)
my %vadr_percaps_H;  # key1: ***flan_caps*** (not vadr_caps) value, value number of instances of flan_caps mapped vadr error in vadr
my %vadr_nomap_H;    # key1: vadr code with no flan_caps mapping, value number of instances of vadr_code
# read in vadr list file, and add vadr caps based on vadr_code2flan_caps_H
open(IN, $vadr_list) || die "ERROR unable to open $vadr_list for reading";
while(my $line = <IN>) { 
  chomp $line;
  #NC_006307 sequence lowsim5s LOW_SIMILARITY_START significant similarity not detected at 5' end of the sequence [low similarity region of length 21]
  #M22015 CDS.NS1 cdsstopn CDS_HAS_STOP_CODON in-frame stop codon exists 5' of stop position predicted by homology to reference [TGA, shifted M:61]
  my @el_A = split(/\s+/, $line);
  my ($seqname, $seq_or_ftr, $code, $vadr_caps, $desc) = ($el_A[0], $el_A[1], $el_A[2], $el_A[3], $el_A[4]);
  #printf("READ VADR code $code\n");
  for(my $i = 5; $i < scalar(@el_A); $i++) { 
    $desc .= " " . $el_A[$i];
  }
  my $flan_caps = lookup_vadr_code2flan_caps(\%vadr_code2flan_caps_H, $code);
  if(! defined $flan_caps) { 
    if(! defined $vadr_nomap_H{$code}) { 
      $vadr_nomap_H{$code} = 0;
    }
    $vadr_nomap_H{$code}++;
    #printf("updated vadr_nomap_H{$code} to $vadr_nomap_H{$code}\n");
  }
  else { 
    if(! defined $vadr_perseq_HHH{$seqname}) { 
      %{$vadr_perseq_HHH{$seqname}} = ();
    }
    if(! defined $vadr_perseq_HHH{$seqname}{$seq_or_ftr}) { 
      %{$vadr_perseq_HHH{$seqname}{$seq_or_ftr}} = ();
    }
    if(! defined $vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps}) { 
      $vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} = 1;      
      if(! defined $vadr_percaps_H{$flan_caps}) { 
        $vadr_percaps_H{$flan_caps} = 0;
      }
      $vadr_percaps_H{$flan_caps}++;
      #printf("updated vadr_percaps_H{$flan_caps} $vadr_percaps_H{$flan_caps}\n");
    }
    elsif($vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} == 1) { 
      # do nothing; we don't count multiple errors for the same seq/ftr
    }
    else {
      die "ERROR unexpected non-1 value for vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps}: " . $vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} . "\n";
    }
  }
}
close(IN);

# determine how many error instances are in common, and unique to flan
my %ncons_percaps_H = ();
my %nflanunq_percaps_H = ();
my %nvadrunq_percaps_H = ();
foreach my $seqname (sort keys %flan_perseq_HHH) { 
  foreach my $seq_or_ftr (sort keys (%{$flan_perseq_HHH{$seqname}})) { 
    foreach my $flan_caps (sort keys (%{$flan_perseq_HHH{$seqname}{$seq_or_ftr}})) { 
      #printf("FLAN CAPS $flan_caps\n");
      if($flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} ne "1") { 
        die "ERROR unexpected non-1 value for flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps}: " . $flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} . "\n";
      }
      if((defined $vadr_perseq_HHH{$seqname}) && 
         (defined $vadr_perseq_HHH{$seqname}{$seq_or_ftr}) && 
         (defined $vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps})) { 
        if($vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} ne "1") { 
          die "ERROR unexpected non-1 value for vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps}: " . $vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} . "\n";
        }
        if(! defined $ncons_percaps_H{$flan_caps}) { 
          $ncons_percaps_H{$flan_caps} = 0;
        }
        $ncons_percaps_H{$flan_caps}++;
        if($do_verbose){ 
          print("CONSISTENT $seqname $seq_or_ftr $flan_caps\n");
        }
        #printf("\tncons_percaps_H{$flan_caps} $ncons_percaps_H{$flan_caps}\n"); 
      }
      else { 
        if(! defined $nflanunq_percaps_H{$flan_caps}) { 
          $nflanunq_percaps_H{$flan_caps} = 0;
        }
        $nflanunq_percaps_H{$flan_caps}++;
        if($do_verbose){ 
          print("FLAN-unique $seqname $seq_or_ftr $flan_caps\n");
        }
        #printf("\tnflanunq_percaps_H{$flan_caps} $nflanunq_percaps_H{$flan_caps}\n"); 
      }
    }
  }
}

# determine how many error instances are unique to vadr
foreach my $seqname (sort keys %vadr_perseq_HHH) { 
  foreach my $seq_or_ftr (sort keys (%{$vadr_perseq_HHH{$seqname}})) { 
    foreach my $flan_caps (sort keys (%{$vadr_perseq_HHH{$seqname}{$seq_or_ftr}})) { 
      if($vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} ne "1") { 
        die "ERROR unexpected non-1 value for vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps}: " . $vadr_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} . "\n";
      }
      if((defined $flan_perseq_HHH{$seqname}) && 
         (defined $flan_perseq_HHH{$seqname}{$seq_or_ftr}) && 
         (defined $flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps})) { 
        if($flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} ne "1") { 
          die "ERROR unexpected non-1 value for flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps}: " . $flan_perseq_HHH{$seqname}{$seq_or_ftr}{$flan_caps} . "\n";
        }
        ; # do nothing we already counted the instances in common
      }
      else { 
        if(! defined $nvadrunq_percaps_H{$flan_caps}) { 
          $nvadrunq_percaps_H{$flan_caps} = 0;
        }
        if($do_verbose){ 
          print("VADR-unique $seqname $seq_or_ftr $flan_caps\n");
        }
        $nvadrunq_percaps_H{$flan_caps}++;
      }
    }
  }
}

my $caption = "\\textbf{Mapping of FLAN and VADR errors and number of instances in the combined training and testing datasets.} Counts are of number of sequence/feature pairs with at least one of the FLAN or VADR error/alert. Some sequence/feature pairs may have multiple errors for the same feature but these are only counted once. Any mapped FLAN and VADR errors are considered consistent if they occur for the same sequence/feature pair. All FLAN errors that cause a sequence to fail (with exceptions detailed in the text) are listed, but not all VADR fatal alerts are.";

print("\\begin{table}[t]\n");
print("\\caption{$caption}\n");
print("\\begin{tabular}{lllrrrrr}\n");
printf("%-45s & %-55s & %-40s & %13s & %13s & %13s & %13s & %13s \\\\ \\hline\n",
       "FLAN error", "VADR error", "VADR alert code", "\\#FLAN", "\\#VADR" , "\\#consistent", "\\#FLAN-unique", "\\#VADR-unique");
my $tot_flan     = 0;
my $tot_vadr     = 0;
my $tot_ncons    = 0; 
my $tot_nflanunq = 0;
my $tot_nvadrunq = 0;
# FLAN error # VADR alert/error # numcases-in-test-flan # numcases-in-test-vadr # numcases-consistent # nunique-flan # nunique-vadr
foreach my $flan_caps (sort keys (%tbl_HH)) { 
#  printf("flan_caps: $flan_caps\n");
  if(! defined $flan_percaps_H{$flan_caps}) { 
    $flan_percaps_H{$flan_caps} = 0;
  }
  if(! defined $vadr_percaps_H{$flan_caps}) { 
    $vadr_percaps_H{$flan_caps} = 0;
  }
  if(! defined $ncons_percaps_H{$flan_caps}) { 
    $ncons_percaps_H{$flan_caps} = 0;
  }
  if(! defined $nflanunq_percaps_H{$flan_caps}) { 
    $nflanunq_percaps_H{$flan_caps} = 0;
  }
  if(! defined $nvadrunq_percaps_H{$flan_caps}) { 
    $nvadrunq_percaps_H{$flan_caps} = 0;
  }

#  printf("%-45s %-55s %-55s %5d %5d %5d %5d %5d\n", 
#         $flan_caps, $tbl_HH{$flan_caps}{"vadr_caps"}, $tbl_HH{$flan_caps}{"vadr_code"}, $flan_percaps_H{$flan_caps}, $vadr_percaps_H{$flan_caps}, 
#         $ncons_percaps_H{$flan_caps}, $nflanunq_percaps_H{$flan_caps}, $nvadrunq_percaps_H{$flan_caps});

  my $flan_caps2print = $flan_caps;
  my $vadr_caps2print = $tbl_HH{$flan_caps}{"vadr_caps"};
  $flan_caps2print =~ s/\_/\\\_/g;
  $vadr_caps2print =~ s/\_/\\\_/g;
  printf("%-45s & %-55s & %-40s & %13d & %13d & %13d & %13d & %13d \\\\ \n",
         $flan_caps2print, $vadr_caps2print, $tbl_HH{$flan_caps}{"vadr_code"}, $flan_percaps_H{$flan_caps}, $vadr_percaps_H{$flan_caps}, 
         $ncons_percaps_H{$flan_caps}, $nflanunq_percaps_H{$flan_caps}, $nvadrunq_percaps_H{$flan_caps});
  if(! defined $flan_words_H{$flan_caps}) {
    die "flan_words_H{$flan_caps} undef\n";
  }
  printf("%-45s & & & & & & & \\\\ \n",
         $flan_words_H{$flan_caps});
  $tot_flan     += $flan_percaps_H{$flan_caps};
  $tot_vadr     += $vadr_percaps_H{$flan_caps};
  $tot_ncons    += $ncons_percaps_H{$flan_caps};
  $tot_nflanunq += $nflanunq_percaps_H{$flan_caps};
  $tot_nvadrunq += $nvadrunq_percaps_H{$flan_caps};
}

print("\\hline\n");

printf("%-45s & %-55s & %-40s & %13d & %13d & %13d & %13d & %13d \\\\ \n",
    "any", "any", "any", $tot_flan, $tot_vadr, $tot_ncons, $tot_nflanunq, $tot_nvadrunq);

print("\\end{tabular}\n");
print("\\label{tbl:errors}\n");
print("\\end{table}\n");


#################################################################
# Subroutine: lookup_vadr_code2flan_caps
# Incept:     EPN, Wed Jan 24 14:19:55 2024
# 
# Purpose:    Determine if $vadr_code is a subseq of any key in 
#             %{$vadr_code2flan_caps_HR}, if so return that value.
#             
# Arguments:
#   $vadr_code2flan_caps_HR: REF to hash of map from vadr codes to flan caps
#   $vadr_code               vadr code to lookup
#
# Returns:    value of matching key, if any, undef if none match
# 
# Dies:       if $vadr_code is subseq of > 1 key
#
#################################################################
sub lookup_vadr_code2flan_caps {
  my $sub_name = "lookup_vadr_code2flan_caps";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($vadr_code2flan_caps_HR, $vadr_code) = @_;

  my $ret_val = undef;
  foreach my $key (sort keys (%{$vadr_code2flan_caps_HR})) { 
    if($key =~ m/$vadr_code/) { 
      if(defined $ret_val) { 
        die "ERROR in $sub_name, $vadr_code matches to more than one key!";
      }
      $ret_val = $vadr_code2flan_caps_HR->{$key};
    }
  }
  return $ret_val;
}
