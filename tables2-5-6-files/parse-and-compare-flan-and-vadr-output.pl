use strict;
use warnings;
use Getopt::Long qw(:config no_auto_abbrev);

require "sqp_ofile.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

my $usage = "perl parse-and-compare-flan-and-vadr-output.pl [OPTIONS] \n\t<flan feature table file>\n\t<vadr 'pass' feature table file>\n\t<vadr 'fail' feature table file>\n\t<vadr .alt file>\n\t<vadr .sqa file>\n\t<file mapping products to genes>\n\t<output root for naming output files>\n\tOPTIONS:\n";
$usage .= "\t\t--misc2cds: in FLAN output, treat misc_feature's as CDS and deduce product name\n";
$usage .= "\t\t--noextra:  in VADR output, do not expect extrant{5,3} alerts\n";
$usage .= "\t\t--deput:    in FLAN output, remove putative in product name\n\n";
$usage .= "\t\t--weird:    output weird FLAN errors for non-annotated gene/CDS\n";

my $do_misc2cds = 0;   # set to '1' if --misc2cds used
my $do_noextra  = 0;   # set to '1' if --noextra is used
my $do_deput    = 0;   # set to '1' if --deput is used
my $do_weird    = 0;   # set to '1' if --weird is used

&GetOptions( "misc2cds" => \$do_misc2cds, 
             "noextra"  => \$do_noextra,
             "deput"    => \$do_deput,
             "weird"    => \$do_weird);
             

if(scalar(@ARGV) != 7) { die $usage; }

my ($flan_ft_file, $vadr_pass_file, $vadr_fail_file, $vadr_alt_file, $vadr_sqa_file, $product_gene_map_file, $out_root) = (@ARGV);

# process the product:gene map file
my %product2gene_H = ();
my %gene2product_H = ();
process_product_gene_map_file($product_gene_map_file, \%product2gene_H, \%gene2product_H);

# process the flan feature table file
my %flan_seqname_H    = (); # key is sequence name, value is '1' indicating this sequence exists in the FLAN feature table
my %flan_warning_HA   = (); # key is sequence name, value is array of WARNING lines for this sequence 
my %flan_error_HA     = (); # key is sequence name, value is array of ERROR lines for this sequence 
my %flan_ftr_info_HAH = (); # the feature info we read from the flan feature table
# the HHA that includes info on all vadr alerts:
# key 1: sequence name
# key 2: <type>.<name> e.g. "CDS.PB1-F2_protein"
# key 3: array of <code>.<fail>.<desc>
#        <code> is 8 letter vadr alert code
#        <fail> is 'yes' if this is a fatal alert, 'no' if not
#        <desc> alert detail from .alt file
my %flan_alt_info_HHA = (); 
my %flan_seq_info_HH  = ();  # 1st dim: seq name, 2nd dim: "type", "segment" or "serotype"
process_flan_feature_table($flan_ft_file, \%flan_seqname_H, \%flan_seq_info_HH, \%flan_ftr_info_HAH, \%product2gene_H, $do_misc2cds, $do_deput);

# get FLAN info (, warnings and errors
process_flan_warnings_and_errors($flan_ft_file, $out_root, \%flan_seqname_H, \%flan_ftr_info_HAH, \%flan_alt_info_HHA);

# combine and process the vadr pass and fail feature table files
my %vadr_seqname_H    = ();
my %vadr_ftr_info_HAH = ();
process_vadr_feature_tables($vadr_pass_file, $vadr_fail_file, \%vadr_seqname_H, \%vadr_ftr_info_HAH, \%product2gene_H, $do_misc2cds);

# the HHA that includes info on all vadr alerts:
# key 1: sequence name
# key 2: <type>.<name> e.g. "CDS.PB1-F2_protein"
# key 3: array of <code>.<fail>.<desc>
#        <code> is 8 letter vadr alert code
#        <fail> is 'yes' if this is a fatal alert, 'no' if not
#        <desc> alert detail from .alt file
my %vadr_alt_info_HHA = (); 
process_vadr_alt_file($vadr_alt_file, $out_root, \%vadr_alt_info_HHA, \%vadr_ftr_info_HAH, \%product2gene_H);

# read the vadr .sqa file to 
my %vadr_seq_info_HH  = ();  # 1st dim: seq name, 2nd dim: "type", "segment" or "serotype"
process_vadr_sqa_file($vadr_sqa_file, \%vadr_seq_info_HH);

# make sure we have the exact same set of sequences in both HAHs
my $seqname; 
foreach $seqname (sort keys %flan_seqname_H) { 
  if(! defined $vadr_seqname_H{$seqname}) { 
    die "ERROR seq $seqname is in $flan_ft_file but not in $vadr_pass_file or $vadr_fail_file";
  }
}
foreach $seqname (sort keys %vadr_seqname_H) { 
  if(! defined $flan_seqname_H{$seqname}) { 
    die "ERROR seq $seqname is in $vadr_pass_file or $vadr_fail_file but not $flan_ft_file";
  }
}

# Go through each feature we care about (CDS, gene, mat_peptide and sig_peptide)
# and check if they are annotated identically in the two feature table files or not.
#
# CDS  is defined by 'product'
# gene is defined by 'gene'
# mat_peptide is defined by 'product'
# sig_peptide has no value in FLAN output, is a mat_peptide with no product in vadr output

my %handled_HHH = ();
foreach my $seqname (sort keys %flan_ftr_info_HAH) { 
  my $seen_flan_sig_peptide_for_this_seq = 0;
  my $seen_vadr_fatal_alt = 0;
  my $seen_flan_fatal_alt = 0;
  for(my $i = 0; $i < scalar(@{$flan_ftr_info_HAH{$seqname}}); $i++) { 
    my $flan_id = undef; 
    my $vadr_id = undef;
    my $flan_type = undef;
    my $vadr_type = undef;
    my $flan_coords = undef;
    my $vadr_coords = undef;
    if(! defined $flan_ftr_info_HAH{$seqname}[$i]{"coords"}) { 
      die "ERROR in $flan_ft_file sequence $seqname feature $i has no coords";
    }
    if($flan_ftr_info_HAH{$seqname}[$i]{"type"} eq "CDS") { 
      if(! defined $flan_ftr_info_HAH{$seqname}[$i]{"product"}) { 
        die "ERROR in $flan_ft_file sequence $seqname CDS has no product";
      }
      $flan_type = "CDS";
      $flan_id   = $flan_ftr_info_HAH{$seqname}[$i]{"product"};
    }
    elsif($flan_ftr_info_HAH{$seqname}[$i]{"type"} eq "gene") { 
      if(! defined $flan_ftr_info_HAH{$seqname}[$i]{"gene"}) { 
        die "ERROR in $flan_ft_file sequence $seqname gene has no gene";
      }
      $flan_type = "gene";
      $flan_id   = $flan_ftr_info_HAH{$seqname}[$i]{"gene"};
    }
    elsif($flan_ftr_info_HAH{$seqname}[$i]{"type"} eq "mat_peptide") { 
      if(! defined $flan_ftr_info_HAH{$seqname}[$i]{"product"}) { 
        die "ERROR in $flan_ft_file sequence $seqname mat_peptide (i: $i) has no product";
      }
      $flan_type = "mat_peptide";
      $flan_id   = $flan_ftr_info_HAH{$seqname}[$i]{"product"};
    }
    elsif($flan_ftr_info_HAH{$seqname}[$i]{"type"} eq "sig_peptide") { 
      if($seen_flan_sig_peptide_for_this_seq) { 
        die "ERROR saw more than one sig_peptide for same seq $seqname in flan output";
      }
      $flan_type = "sig_peptide";
      $flan_id   = "sig_peptide";
      $seen_flan_sig_peptide_for_this_seq = 1;
    }
    else { 
      die "ERROR found unexpected flan feature for seq $seqname type: " . $flan_ftr_info_HAH{$seqname}[$i]{"type"} . "\n";
    }
    $flan_coords = $flan_ftr_info_HAH{$seqname}[$i]{"coords"};

    # find this feature in the vadr feature table, if it exists
    my $nftr_vadr = scalar(@{$vadr_ftr_info_HAH{$seqname}});
    for(my $j = 0; $j < $nftr_vadr; $j++) { 
      $vadr_type = undef;
      $vadr_id = undef;
      if($vadr_ftr_info_HAH{$seqname}[$j]{"type"} eq $flan_type) { 
        $vadr_type = $flan_type;
        if($vadr_type eq "CDS") { 
          if(! defined $vadr_ftr_info_HAH{$seqname}[$j]{"product"}) { 
            die "ERROR in VADR sequence $seqname CDS has no product";
          }
          $vadr_id = $vadr_ftr_info_HAH{$seqname}[$j]{"product"};
        }
        elsif($vadr_type eq "gene") { 
          if(! defined $vadr_ftr_info_HAH{$seqname}[$j]{"gene"}) { 
            die "ERROR in VADR sequence $seqname gene has no gene";
          }
          $vadr_id = $vadr_ftr_info_HAH{$seqname}[$j]{"gene"};
        }
        elsif($vadr_type eq "mat_peptide") { 
          if(defined $vadr_ftr_info_HAH{$seqname}[$j]{"product"}) { 
            $vadr_id = $vadr_ftr_info_HAH{$seqname}[$j]{"product"};
          }
        }
        elsif($vadr_type eq "sig_peptide") { 
          $vadr_id = "sig_peptide";
        }
        else { 
          die "ERROR found unexpected vadr feature for seq $seqname type: " . $vadr_ftr_info_HAH{$seqname}[$i]{"type"} . "\n";
        }
      }
      # handle case where flan feature is sig_peptide, it could match a vadr mat_peptide
      elsif(($flan_type eq "sig_peptide") && 
            ($vadr_ftr_info_HAH{$seqname}[$j]{"type"} eq "mat_peptide") && 
            (! defined $vadr_ftr_info_HAH{$seqname}[$j]{"product"})) { 
        $vadr_type = "sig_peptide";
        $vadr_id = "sig_peptide";
      }
      # check if ftr $j for vadr matches ftr $i for flan
      if((defined $vadr_id) && 
         (defined $vadr_type) && 
         ($vadr_type eq $flan_type) && 
         ($vadr_id eq $flan_id)) { 
        if(defined $handled_HHH{$seqname}{$flan_type}{$flan_id}) { 
          die "ERROR found two matches for $seqname $flan_type $flan_id";
        }
        $vadr_coords = $vadr_ftr_info_HAH{$seqname}[$j]{"coords"};
        
        my ($flan_alts, $vadr_alts) = 
            get_output_for_warnings_errors_alerts(\%{$flan_ftr_info_HAH{$seqname}[$i]}, 
                                                  \%{$vadr_ftr_info_HAH{$seqname}[$j]}, 
                                                  (defined $flan_alt_info_HHA{$seqname}) ? \%{$flan_alt_info_HHA{$seqname}} : undef, 
                                                  (defined $vadr_alt_info_HHA{$seqname}) ? \%{$vadr_alt_info_HHA{$seqname}} : undef, 
                                                  \%product2gene_H);

        my $coords_compare_string = get_coords_compare_string($flan_coords, $vadr_coords);
        my $alts_compare_string   = get_feature_alert_compare_string($flan_alts, $vadr_alts);
        printf("FEATURE\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $seqname, $flan_type, $flan_id, $flan_coords, $vadr_coords, 
               $coords_compare_string, $flan_alts, $vadr_alts, $alts_compare_string);
        if($flan_alts =~ m/\:\:fatal\:\:/) { 
          $seen_flan_fatal_alt = 1;
        }
        if($vadr_alts =~ m/\:\:fatal\:\:/) { 
          $seen_vadr_fatal_alt = 1;
        }
        $handled_HHH{$seqname}{$flan_type}{$flan_id} = 1;
      }
    }
    if(! defined $handled_HHH{$seqname}{$flan_type}{$flan_id}) { 
      # no vadr ftr matched flan ftr $i
      my ($flan_alts, undef) = 
          get_output_for_warnings_errors_alerts(\%{$flan_ftr_info_HAH{$seqname}[$i]}, 
                                                undef, 
                                                (defined $flan_alt_info_HHA{$seqname}) ? \%{$flan_alt_info_HHA{$seqname}} : undef, 
                                                undef,
                                                \%product2gene_H);
      my $alts_compare_string   = get_feature_alert_compare_string($flan_alts, "-");
      printf("FEATURE\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $seqname, $flan_type, $flan_id, $flan_coords, "-", "flanonly", $flan_alts, "-", $alts_compare_string);
      if($flan_alts =~ m/\:\:fatal\:\:/) { 
        $seen_flan_fatal_alt = 1;
      }
    }
  } # end of for $i loop

  # check if there's any vadr features we didn't have a matching flan feature for, and output those
  my $seen_vadr_sig_peptide_for_this_seq = 0;
  for(my $j = 0; $j < scalar(@{$vadr_ftr_info_HAH{$seqname}}); $j++) { 
    my $vadr_id = undef;
    my $vadr_type = undef;
    my $vadr_coords = undef;
    if(! defined $vadr_ftr_info_HAH{$seqname}[$j]{"coords"}) { 
      die "ERROR in VADR sequence $seqname feature $j has no coords";
    }
    if($vadr_ftr_info_HAH{$seqname}[$j]{"type"} eq "CDS") { 
      if(! defined $vadr_ftr_info_HAH{$seqname}[$j]{"product"}) { 
        die "ERROR in VADR sequence $seqname CDS has no product";
      }
      $vadr_type = "CDS";
      $vadr_id   = $vadr_ftr_info_HAH{$seqname}[$j]{"product"};
    }
    elsif($vadr_ftr_info_HAH{$seqname}[$j]{"type"} eq "gene") { 
      if(! defined $vadr_ftr_info_HAH{$seqname}[$j]{"gene"}) { 
        die "ERROR in VADR sequence $seqname gene has no gene";
      }
      $vadr_type = "gene";
      $vadr_id   = $vadr_ftr_info_HAH{$seqname}[$j]{"gene"};
    }
    elsif(($vadr_ftr_info_HAH{$seqname}[$j]{"type"} eq "mat_peptide") && 
          (defined $vadr_ftr_info_HAH{$seqname}[$j]{"product"})) { 
      $vadr_type = "mat_peptide";
      $vadr_id   = $vadr_ftr_info_HAH{$seqname}[$j]{"product"};
    }
    elsif($vadr_ftr_info_HAH{$seqname}[$j]{"type"} eq "sig_peptide") { 
      if($seen_vadr_sig_peptide_for_this_seq) { 
        die "ERROR saw more than one special mat_peptide for same seq $seqname in vadr output";
      }
      $vadr_id   = "sig_peptide";
      $seen_vadr_sig_peptide_for_this_seq = 1;
    }
    elsif(($vadr_ftr_info_HAH{$seqname}[$j]{"type"} eq "mat_peptide") && 
          (! defined $vadr_ftr_info_HAH{$seqname}[$j]{"product"})) { 
      if($seen_vadr_sig_peptide_for_this_seq) { 
        die "ERROR saw more than one special mat_peptide for same seq $seqname in vadr output";
      }
      $vadr_type = "sig_peptide";
      $vadr_id   = "sig_peptide";
      $seen_vadr_sig_peptide_for_this_seq = 1;
    }
    else { 
      die "ERROR found unexpected vadr feature for seq $seqname type: " . $vadr_ftr_info_HAH{$seqname}[$j]{"type"} . "\n";
    }
    if((defined $vadr_type) && (defined $vadr_id)) { 
      if(! defined $handled_HHH{$seqname}{$vadr_type}{$vadr_id}) { 
        my ($vadr_alts, undef) = 
            get_output_for_warnings_errors_alerts(undef, 
                                                  \%{$vadr_ftr_info_HAH{$seqname}[$j]}, 
                                                  undef, 
                                                  (defined $vadr_alt_info_HHA{$seqname}) ? \%{$vadr_alt_info_HHA{$seqname}} : undef, 
                                                  \%product2gene_H);
        $vadr_coords = $vadr_ftr_info_HAH{$seqname}[$j]{"coords"};
        my $alts_compare_string = get_feature_alert_compare_string("-", $vadr_alts);
        printf("FEATURE\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $seqname, $vadr_type, $vadr_id, "-", $vadr_coords, "vadronly", "-", $vadr_alts, $alts_compare_string);
        if($vadr_alts =~ m/\:\:fatal\:\:/) { 
          $seen_vadr_fatal_alt = 1;
        }
      }
    }
  }
  # output the 'per-sequence' line
  my ($flan_alts, $vadr_alts) = 
      get_output_for_warnings_errors_alerts(undef, 
                                            undef, 
                                            (defined $flan_alt_info_HHA{$seqname}) ? \%{$flan_alt_info_HHA{$seqname}} : undef, 
                                            (defined $vadr_alt_info_HHA{$seqname}) ? \%{$vadr_alt_info_HHA{$seqname}} : undef, 
                                            \%product2gene_H);
  if($vadr_alts =~ m/\:\:fatal\:\:/) { 
    $seen_vadr_fatal_alt = 1;
  }
  if($flan_alts =~ m/\:\:fatal\:\:/) { 
    $seen_flan_fatal_alt = 1;
  }
  my $flan_seq_info_str = $flan_seq_info_HH{$seqname}{"type"} . ":" . $flan_seq_info_HH{$seqname}{"segment"} . ":";
  if(defined $flan_seq_info_HH{$seqname}{"serotype"}) {
    $flan_seq_info_str .= $flan_seq_info_HH{$seqname}{"serotype"};
  }
  else { 
    $flan_seq_info_str .= "-";
  }

  my $vadr_seq_info_str = $vadr_seq_info_HH{$seqname}{"type"} . ":" . $vadr_seq_info_HH{$seqname}{"segment"} . ":";
  if(defined $vadr_seq_info_HH{$seqname}{"serotype"}) {
    $vadr_seq_info_str .= $vadr_seq_info_HH{$seqname}{"serotype"};
  }
  else { 
    $vadr_seq_info_str .= "-";
  }

  my $seq_info_compare_str = ($flan_seq_info_str eq $vadr_seq_info_str) ? "info:identical" : "info:different";

  printf("SEQUENCE\t%s\t%s\t%s\t%s%s%s%s\t%s\t%s\t%s\n", $seqname, "F:" . $flan_seq_info_str, "V:" . $vadr_seq_info_str, "F", ($seen_flan_fatal_alt ? "F" : "P"), "V", ($seen_vadr_fatal_alt ? "F" : "P"), $flan_alts, $vadr_alts, $seq_info_compare_str);
}
    
exit 0;

#################################################################
# Subroutine:  process_flan_feature_table()
# Incept:      EPN, Wed May  3 11:37:59 2023
#
# Purpose:    Parse and save info in a FLAN feature table. Create a 
#             'processed' version of the file that sqf_FeatureTableParse() 
#             can parse and fill %flan_ftr_info_HAH.
#
# Arguments: 
#  $flan_ft_file:          the FLAN feature table file
#  $flan_seqname_HR:       reference to an 'exists' hash of the FLAN sequence names
#  $flan_seq_info_HHR:     ref to two dim hash, first key seq name, second either "type", "segment" or "serotype"
#  $flan_ftr_info_HAHR:    ref to feature info 
#  $product2gene_HR:       ref to hash mapping products to genes
#  $do_misc2cds:           TRUE to treat misc_feature as CDS
#  $do_deput:              TRUE to remove 'putative' from product names
#
# Returns:    void
#
################################################################# 
sub process_flan_feature_table { 
  my $sub_name = "process_flan_feature_table";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($flan_ft_file, $flan_seqname_HR, $flan_seq_info_HHR, $flan_ftr_info_HAHR, $product2gene_HR, $do_misc2cds, $do_deput) = @_;

  my $p_flan_ft_file = $flan_ft_file . ".p";

  open(INFLAN,       $flan_ft_file)   || die "ERROR unable to open $flan_ft_file for reading";
  open(OUTFLAN, ">", $p_flan_ft_file) || die "ERROR unable to open $p_flan_ft_file for writing";

  my $line;
  my $seqname;
  while($line = <INFLAN>) { 
    chomp $line;
    # remove trailing tabs
    $line =~ s/\t+$//;

    # replace spaces after some 'product' values
    if($line =~ /(^.*\t)product (\S+)$/) { 
      $line = $1 . "product\t" . $2;
    } 

    if($line !~ m/./) { 
      ; # skip blank lines;
    }
    elsif($line =~ /^\>/) { 
      if($line =~ m/^\>Feature\s+(\S+)/) { 
        $seqname = $1;
      }
      elsif($line =~ /^\>(\S+).+/) { 
        $seqname = $1;
        $line = ">Feature $seqname"
      }
      else { 
        die "ERROR unable to parse FLAN line starting with >:$line\n";
      }
      $seqname = shorten_seqname_no_version($seqname);
      if(defined $flan_seqname_HR->{$seqname}) { 
        die "ERROR read $seqname twice in $flan_ft_file";
      }
      
      $flan_seqname_HR->{$seqname} = 1;
      %{$flan_seq_info_HHR->{$seqname}} = ();
      print OUTFLAN $line . "\n";
    }
    elsif($line =~ /^\s+INFO\:\s+Serotype:\s+(.+)$/) { 
      #INFO: Serotype: N1
      if(defined $flan_seq_info_HHR->{$seqname}{"serotype"}) { 
        die "ERROR read serotype for $seqname twice in $flan_ft_file";
      }
      $flan_seq_info_HHR->{$seqname}{"serotype"} = $1;
    }
    elsif($line =~ /^\s+INFO\:\s+Virus type:\s+influenza\s+(\S)$/) { 
      #INFO: Virus type: influenza A
      if(defined $flan_seq_info_HHR->{$seqname}{"type"}) { 
        die "ERROR read type for $seqname twice in $flan_ft_file";
      }
      $flan_seq_info_HHR->{$seqname}{"type"} = $1;
    }
    elsif($line =~ /^\s+INFO\:\s+Segment:\s+(\d)/) { 
      #INFO: Segment: 4 (HA)
      if(defined $flan_seq_info_HHR->{$seqname}{"segment"}) { 
        die "ERROR read segment for $seqname twice in $flan_ft_file";
      }
      $flan_seq_info_HHR->{$seqname}{"segment"} = $1;
    }
    elsif($line =~ /^\s+INFO\:/) { 
      ; # do nothing
    }
    elsif($line =~ /^\s+WARNING\:\s+(.+)$/) { 
      ; # do nothing, we'll store this later in process_flan_warnings_and_errors()
    }
    elsif($line =~ /^\s+ERROR\:\s+(.+)$/) { 
      ; # do nothing, we'll store this later in process_flan_warnings_and_errors()
    }
    else { 
      print OUTFLAN $line . "\n";
    }
  }
  close(INFLAN);
  close(OUTFLAN);

  sqf_FeatureTableParse($p_flan_ft_file, $flan_ftr_info_HAHR, undef);
  truncate_coords($flan_ftr_info_HAHR);
  #utl_HAHDump("flan_ftr_info", $flan_ftr_info_HAHR, *STDOUT);

  if($do_misc2cds) { 
    foreach my $seqname (sort keys %flan_ftr_info_HAH) { 
      for(my $i = 0; $i < scalar(@{$flan_ftr_info_HAHR->{$seqname}}); $i++) { 
        if($flan_ftr_info_HAHR->{$seqname}[$i]{"type"} eq "misc_feature") { 
          if(defined $flan_ftr_info_HAHR->{$seqname}[$i]{"note"}) { 
            if($flan_ftr_info_HAHR->{$seqname}[$i]{"note"} =~ /^nonfunctional (.+) due to mutation/) { 
              #printf("note: " . $flan_ftr_info_HAHR->{$seqname}[$i]{"note"} . "\n");
              my $product = $1;
              $product =~ s/^putative //;
              $flan_ftr_info_HAHR->{$seqname}[$i]{"type"}    = "CDS";
              $flan_ftr_info_HAHR->{$seqname}[$i]{"product"} = $product;
              if(! defined $product2gene_HR->{$product}) { 
                die "ERROR in $sub_name, trying to change misc_feature to CDS but product $product not in product2gene map";
              }
              $flan_ftr_info_HAHR->{$seqname}[$i]{"gene"} = $product2gene_HR->{$product};
              delete $flan_ftr_info_HAHR->{$seqname}[$i]{"note"};
            }
            else { 
              die "ERROR unparseable note for misc_feature: " . $flan_ftr_info_HAHR->{$seqname}[$i]{"note"};
            }
          }
          else { 
            die "ERROR no note for misc_feature";
          }
        }
      }
    }
  }

  if($do_deput) { 
    foreach my $seqname (sort keys %{$flan_ftr_info_HAHR}) { 
      for(my $i = 0; $i < scalar(@{$flan_ftr_info_HAHR->{$seqname}}); $i++) { 
        if($flan_ftr_info_HAHR->{$seqname}[$i]{"type"} eq "CDS") { 
          if(defined $flan_ftr_info_HAHR->{$seqname}[$i]{"product"}) { 
            $flan_ftr_info_HAHR->{$seqname}[$i]{"product"} =~ s/^putative\s+//;
          }
        }
      }
    }
  }
     
  #utl_HAHDump("flan_ftr_info", $flan_ftr_info_HAHR, *STDOUT);

}

#################################################################
# Subroutine:  process_flan_warnings_and_errors()
# Incept:      EPN, Thu May  4 10:42:16 2023
#
# Purpose:    Parse and save warning and error info in a FLAN feature table.
#
# Arguments: 
#  $flan_ft_file:          the FLAN feature table file
#  $out_root:              output root for output files
#  $flan_seqname_HR:       hash of sequence names, already filled
#  $flan_ftr_info_HAHR:    the feature info, already filled
#  $flan_alt_info_HHAR:    ref to 2d hash of arrays with information on warnings (non-fatal) and errors (fatal), FILLED HERE
#
# Returns:    void
#
################################################################# 
sub process_flan_warnings_and_errors {
  my $sub_name = "process_flan_warnings_and_errors";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($flan_ft_file, $out_root, $flan_seqname_HR, $flan_ftr_info_HAHR, $flan_alt_info_HHAR) = @_;

  my $flan_err_file  = $out_root . ".flan.err.list";
  my $flan_wrn_file  = $out_root . ".flan.warn.list";

  open(INFLAN, $flan_ft_file) || die "ERROR unable to open $flan_ft_file for reading";
  open(ERRFLAN, ">", $flan_err_file)  || die "ERROR unable to open $flan_err_file for writing";
  open(WRNFLAN, ">", $flan_wrn_file)  || die "ERROR unable to open $flan_wrn_file for writing";

  my $line;
  my $seqname;
  my $saw_frameshift_flag = 0;
  my $frameshift_info = "";
  while($line = <INFLAN>) { 
    chomp $line;

    # replace spaces after some 'product' values
    if($line =~ /(^.*\t)product (\S+)$/) { 
      $line = $1 . "product\t" . $2;
    } 

    if($line !~ m/./) { 
      ; # skip blank lines;
    }
    elsif($line =~ /^\>/) { 
      if($line =~ m/^\>Feature\s+(\S+)/) { 
        $seqname = $1;
      }
      elsif($line =~ /^\>(\S+).+/) { 
        $seqname = $1;
      }
      else { 
        die "ERROR unable to parse FLAN line starting with >:$line\n";
      }
      $seqname = shorten_seqname_no_version($seqname);
      if(! defined $flan_seqname_HR->{$seqname}) { 
        die "ERROR read new sequence $seqname on second pass";
      }
    }
    elsif($line =~ /^\s+INFO\:\s+(.+)$/) { 
      ; # do nothing
    }
    elsif($line =~ /^\s+WARNING\:\s+(.+)$/) { 
      my $warning = $1;
      # classify the warning, we don't store all of them
      my $do_push  = undef;
      my $type_id  = undef;
      my $alt_code = undef;
      #!# IGNORED WARNING flan:NONE:This sequence lacks \d+ nucleotides of coding sequences at the beginning
      if($warning =~ /^This sequence lacks \d+ nucleotides of coding sequences at the beginning$/) { 
        $do_push = 0;
      }
      #!# IGNORED WARNING flan:NONE:This sequence lacks \d+ nucleotides of coding sequences at the end
      elsif($warning =~ /^This sequence lacks \d+ nucleotides of coding sequences at the end$/) { 
        $do_push = 0;
      }
      #!# IGNORED WARNING flan:NONE:Ambiguity nucleotide(s) found:/) { 
      elsif($warning =~ /^Ambiguity nucleotide\(s\) found:/) { 
        $do_push = 0;
      }
      #!# WARNING vadr:MUTATION_AT_START:mutstart == flan:MUTATION_AT_START:Probable mutation at Start of prot(\d+)
      elsif($warning =~ /^Probable mutation at Start of prot(\d+)/) { 
        my $cds_idx = $1;
        my $gene = find_flan_cds_gene_given_idx(\@{$flan_ftr_info_HAHR->{$seqname}}, $cds_idx);
        if($gene ne "") { 
          $type_id = "CDS." . $gene;
          $alt_code = "mutstart";
          $do_push = 1;
        }
        else { # we didn't find this feature so we ignore this (this happens for HE584760.1)
          $do_push = 0;
        }
      }
      #!# WARNING vadr:MUTATION_AT_END:mutendcd,mutendns,mutendex == flan:MUTATION_AT_END:Probable mutation at End of prot<d>
      elsif($warning =~ /^Probable mutation at End of prot(\d+)/) { 
        my $cds_idx = $1;
        my $gene = find_flan_cds_gene_given_idx(\@{$flan_ftr_info_HAHR->{$seqname}}, $cds_idx);
        if($gene ne "") { 
          $type_id = "CDS." . $gene;
          $alt_code = "mutendcd,mutendns,mutendex";
          $do_push = 1;
        }
        else { # we didn't find this feature so we ignore this (this happens for HE584760.1)
          $do_push = 0;
        }
      }
      #!# WARNING vadr:MUTATION_AT_START:mutstart == flan:NONE:The coding region of <CDS> has no start codon
      elsif($warning =~ /^The coding region of (\S+) has no start codon$/) { 
        # only store this if the feature is not 5' truncated
        my $gene = $1;
        my $is_5trunc = undef;
        # only call check_flan_trunc if we have annotation for at least one CDS
        # sometimes we don't even though we should (FLAN doesn't annotate any features when there is a frameshift in flu B seqs for example)
        if((defined $flan_ftr_info_HAHR->{$seqname}) && (scalar(@{$flan_ftr_info_HAHR->{$seqname}}) > 0)) { 
          ($is_5trunc, undef) = check_flan_trunc($flan_ftr_info_HAHR, $seqname, "CDS", $gene);
          if(! $is_5trunc) { 
            $do_push  = 1;
            $alt_code = "mutstart";
            $type_id  = "CDS." . $gene;
          }
          else { 
            $do_push = 0;
          }
        }
        else { 
          $do_push = 0;
        }
      }
      #!# WARNING vadr:MUTATION_AT_END:mutendcd,mutendns,mutendex == flan:NONE:The coding region of <CDS> has no stop codon
      elsif($warning =~ /^The coding region of (\S+) has no stop codon$/) { 
        # only store this if the feature is not 3' truncated
        my $gene = $1;
        my $is_3trunc = undef;
        # only call check_flan_trunc if we have annotation for at least one CDS
        # sometimes we don't even though we should (FLAN doesn't annotate any features when there is a frameshift in flu B seqs for example)
        if((defined $flan_ftr_info_HAHR->{$seqname}) && (scalar(@{$flan_ftr_info_HAHR->{$seqname}}) > 0)) { 
          (undef, $is_3trunc) = check_flan_trunc($flan_ftr_info_HAHR, $seqname, "CDS", $gene);
          if(! $is_3trunc) { 
            $do_push  = 1;
            $alt_code = "mutendcd,mutendns,mutendex";
            $type_id  = "CDS." . $gene;
          }
          else { 
            $do_push = 0;
          }
        }
        else { 
          $do_push = 0;
        }
      }
      #!# WARNING vadr:CDS_HAS_STOP_CODON:cdsstopn == flan:CDS_HAS_STOP_CODON:The coding region of <CDS> has stop codon inside exon
      elsif($warning =~ /^The coding region of (.+) has stop codon inside exon$/) { 
        $type_id  = "CDS." . $1;
        $do_push  = 1;
        $alt_code = "cdsstopn";
      }
      #!# WARNING vadr:POSSIBLE_FRAMESHIFT_HIGH_CONF,POSSIBLE_FRAMESHIFT:fsthicft,fsthicfi,fstukcft,fstukcfi == flan:CDS_HAS_FRAMESHIFT:The coding region of (CDS) has frameshift
      elsif($warning =~ /^The coding region of (.+) has frameshift$/) { 
        $type_id  = "CDS." . $1; 
        $do_push  = 1;
        $alt_code = "fsthicft,fsthicfi,fstukcft,fstukcfi";
        # FOR WARNINGS DON'T add frameshift info we read above for deletions or insertions
        # (unless future instances force us to, only ones I've seen so far are PB1-F2 which 
        # overlaps with PB1...)
        #$error .= ";" . $frameshift_info;
        #$frameshift_info = "";
        #$saw_frameshift_flag = 1;
      }
      #!# WARNING vadr:DELETION_OF_FEATURE:deletins == flan:NONE:Protein <CDS> was not found
      # MAYBE THIS SHOULD BE deletina? 
      elsif($warning =~ /^Protein (.+) was not found$/) { 
        $type_id  = "CDS." . $1;
        $do_push  = 1;
        $alt_code = "deletins";
      }
      #!# WARNING vadr:MUTATION_AT_SPLICE_SITE:mutspst5,mutspspt3 == flan:SPLICE_SITE_NOT_FOUND:Splice site was not found for segment <d> protein <CDS>
      elsif($warning =~ /^Splice site was not found for segment \d+ protein (.+)$/) { 
        $type_id  = "CDS." . $1;
        $do_push  = 1;
        $alt_code = "mutspst5,mutspspt3";
      }
      if(! defined $do_push) { 
        die "ERROR did not recognize WARNING line: $line\n";
      }
      if($do_push) { 
        if(! defined $type_id) { 
          die "ERROR type_id not set for pushable warning on line: $line\n";
        }
        if(! defined $alt_code) { 
          die "ERROR alt_code not set for pushable warning on line: $line\n";
        }
        if(! defined $flan_alt_info_HHAR->{$seqname}) { 
          %{$flan_alt_info_HHAR->{$seqname}} = ();
        }
        if(! defined $flan_alt_info_HHAR->{$seqname}{$type_id}) { 
          @{$flan_alt_info_HHAR->{$seqname}{$type_id}} = ();
        }
        my $code_fail_desc = $alt_code . "::nonfatal::" . $warning;
        push(@{$flan_alt_info_HHAR->{$seqname}{$type_id}}, $code_fail_desc);
      }
    }
    elsif($line =~ /^\s+ERROR\:\s+(.+)$/) { 
      my $error = $1;
      # classify the error, we don't store all of them
      my $do_push  = undef;
      my $is_fatal = undef;
      my $type_id  = undef;
      my $alt_code = undef;
      #!# IGNORED ERROR flan:SEQUENCE_IS_TOO_SHORT:Input sequence is too short
      if($error =~ /^Input sequence is too short$/) { 
        $type_id  = "sequence";
        $do_push  = 0;
        $is_fatal = 0;
        if($do_weird) { 
          printf("WEIRD: $seqname has TOO_SHORT error we are ignoring\n");
        }
      }
      #!# ERROR vadr:CDS_HAS_STOP_CODON:cdsstopn == flan:CDS_HAS_STOP_CODON:The coding region of <CDS> has stop codon inside exon
      elsif($error =~ /^The coding region of (.+) has stop codon inside exon$/) { 
        my $cds_id = $1;
        $type_id  = "CDS." . $cds_id;
        $do_push  = 1;
        $alt_code = "cdsstopn";
        $is_fatal = 1;
        #$is_fatal = $cds_id eq "PA-X" ? 0 : 1; # for some reason FLAN reports this particular ERROR for PA-X but it shouldn't be fatal afaik
      }
      #!# ERROR vadr:NO_FEATURES_ANNOTATED:noftrann,noftrant == flan:CODING_CAPACITY:This sequence does not have coding capacity
      elsif($error =~ /^This sequence does not have coding capacity/) {
        $type_id  = "sequence";
        $do_push  = 1;
        $alt_code = "noftrann,noftrant";
        $is_fatal = 1;
      }
      #!# ERROR vadr:MUTATION_AT_START:mutstart == flan:MUTATION_AT_START:Probable mutation at Start of prot(\d+)
      elsif($error =~ /^Probable mutation at Start of prot(\d+)/) { 
        my $cds_idx = $1;
        my $gene = find_flan_cds_gene_given_idx(\@{$flan_ftr_info_HAHR->{$seqname}}, $cds_idx);
        if($gene ne "") { 
          $type_id = "CDS." . $gene;
          $alt_code = "mutstart";
          $do_push = 1;
          $is_fatal = 1;
        }
        else { # we didn't find this feature so we ignore this (this happens for HE584760.1)
          $do_push = 0;
          $is_fatal = 0;
          if($do_weird) { 
            printf("WEIRD: $seqname has Probable mutation at Start of prot error we are ignoring\n");
          }
        }
      }
      #!# ERROR vadr:MUTATION_AT_END:mutendcd,mutendns,mutendex == flan:MUTATION_AT_END:Probable mutation at End of prot<d>
      elsif($error =~ /^Probable mutation at End of prot(\d+)/) { 
        my $cds_idx = $1;
        my $gene = find_flan_cds_gene_given_idx(\@{$flan_ftr_info_HAHR->{$seqname}}, $cds_idx);
        if($gene ne "") { 
          $type_id = "CDS." . $gene;
          $alt_code = "mutendcd,mutendns,mutendex";
          $do_push = 1;
          $is_fatal = 1;
        }
        else { # we didn't find this feature so we ignore this (this happens for HE584760.1)
          $do_push = 0;
          $is_fatal = 0;
          if($do_weird) { 
            printf("WEIRD: $seqname has Probable mutation at End of prot error we are ignoring\n");
          }
        }
      }
      #!# ERROR vadr:NO_ANNOTATION:noannotn == flan:NO_BLAST_HITS_FOUND:No blast hits found!
      elsif($error =~ /No blast hits found!$/) { 
        $type_id  = "sequence";
        $do_push  = 1;
        $alt_code = "noannotn";
        $is_fatal = 1;
      }
      #!# IGNORED ERROR flan:DELETION_OF_NT:Deletion of <d+> nt around: <d> - <d>
      elsif($error =~ /^Deletion of \d+ nt around: \d+ - \d+$/) { 
        # we actually don't ignore this, we store it to add to the next frameshift
        $frameshift_info .= $error . ",";
        $do_push = 0;
        $is_fatal = 0;
      }
      #!# IGNORED ERROR flan:INSERTION_OF_NT:Insertion of <d+> nt around: <d> - <d>
      elsif($error =~ /^Insertion of \d+ nt around: \d+ - \d+$/) { 
        # we actually don't ignore this, we store it to add to the next frameshift
        $frameshift_info .= $error . ",";
        $do_push = 0;
        $is_fatal = 0;
      }
      #!# ERROR vadr:POSSIBLE_FRAMESHIFT_HIGH_CONF,POSSIBLE_FRAMESHIFT:fsthicft,fsthicfi,fstukcft,fstukcfi == flan:CDS_HAS_FRAMESHIFT:The coding region of (CDS) has frameshift
      elsif($error =~ /^The coding region of (.+) has frameshift$/) { 
        $type_id  = "CDS." . $1; 
        $do_push  = 1;
        $alt_code = "fsthicft,fsthicfi,fstukcft,fstukcfi";
        # add frameshift info we read above for deletions or insertions
        $error .= ";" . $frameshift_info;
        $frameshift_info = "";
        $saw_frameshift_flag = 1;
        $is_fatal = 1;
      }
      #!# ERROR vadr:POSSIBLE_FRAMESHIFT_HIGH_CONF,POSSIBLE_FRAMESHIFT:fsthicft,fsthicfi,fstukcft,fstukcfi == flan:PEPTIDE_FRAMESHIFT:Mature peptide (<name>) has frameshift
      elsif($error =~ /^Mature peptide \((.+)\) has frameshift$/) { 
        $type_id  = "mat_peptide." . $1; 
        $do_push  = 1;
        $alt_code = "fsthicft,fsthicfi,fstukcft,fstukcfi";
        $is_fatal = 1;
      }
      #!# ERROR vadr:PEPTIDE_ADJACENCY_PROBLEM:pepadjcy == flan:PEPTIDE_SEPARATED:Mature peptides (<name>) and <name> are separated
      elsif($error =~ /^Mature peptides \((.+)\) and \(.+\) are separated$/) { 
        $type_id  = "mat_peptide." . $1; 
        $do_push  = 1;
        $alt_code = "pepadjcy";
        $is_fatal = 1;
      }
      #!# ERROR vadr:PEPTIDE_ADJACENCY_PROBLEM:pepadjcy == flan:PEPTIDE_OVERLAP:Mature peptides (<name>) and <name> overlap
      elsif($error =~ /^Mature peptides \((.+)\) and \(.+\) overlap$/) { 
        $type_id  = "mat_peptide." . $1; 
        $do_push  = 1;
        $alt_code = "pepadjcy";
        $is_fatal = 1;
      }
      #!# ERROR vadr:MUTATION_AT_SPLICE_SITE:mutspst5,mutspst3 == flan:SPLICE_SITE_NOT_FOUND:Expected splice site consensus sequence not found for protein <d>.
      elsif($error =~ /^Expected splice site consensus sequence not found for protein (\d+)\.$/) { 
        my $cds_idx = $1;
        my $gene = find_flan_cds_gene_given_idx(\@{$flan_ftr_info_HAHR->{$seqname}}, $cds_idx);
        if($gene ne "") { 
          $type_id = "CDS." . $gene;
          $do_push  = 1;
          $alt_code = "mutspst5,mutspst3";
          $is_fatal = 1;
        }
        else { 
          $do_push = 0;
          $is_fatal = 0;
          if($do_weird) { 
            printf("WEIRD: $seqname has Expected splice site consensus sequence not found for protein we are ignoring\n");
          }
        }
      }
      #!# ERROR vadr:NONEWRONG:nonewrong == flan:WRONG_SEGMENT:Wrong exon number 2 for segment 3 protein (PA)
      elsif($error =~ /^Wrong exon number \d+ for segment 3 protein (.+)$/) { 
        $type_id  = "CDS." . $1;
        $do_push  = 1;
        $alt_code = "nonewrong";
        $is_fatal = 1;
      }
      #!# ERROR vadr:EXTRA_SEQUENCE_START:extrant5 == flan:CONTAMINATION_UPSTREAM:Sequence \((\S+)\) contains extra 1 nts upstream the consensus 5' end sequence of influenza viruses. Check for possible vector\/linker contamination.
      #!# if --noextra: vadr:lowsim5s == flan:Sequence \((\S+)\) contains extra 1 nts upstream the consensus 5' end sequence of influenza viruses. Check for possible vector\/linker contamination.
      elsif($error =~ /^Sequence \((\S+)\) contains extra \d+ nts upstream the consensus 5' end sequence of influenza viruses. Check for possible vector\/linker contamination.$/) { 
        $type_id  = "sequence";
        $do_push  = 1;
        $alt_code = ($do_noextra) ? "lowsim5s" : "extrant5";
        $is_fatal = 1;
      }
      #!# ERROR vadr:EXTRA_SEQUENCE_END:extrant3 == flan:CONTAMINATION_DOWNSTREAM:Sequence \((\S+)\) contains extra 1 nts upstream the consensus 3' end sequence of influenza viruses. Check for possible vector\/linker contamination.
      #!# if --noextra: vadr:lowsim3s == flan:Sequence \((\S+)\) contains extra 1 nts upstream the consensus 3' end sequence of influenza viruses. Check for possible vector\/linker contamination.
      elsif($error =~ /^Sequence \((\S+)\) contains extra \d+ nts downstream the consensus 3' end sequence of influenza viruses. Check for possible vector\/linker contamination.$/) { 
        $type_id  = "sequence";
        $do_push  = 1;
        $alt_code = ($do_noextra) ? "lowsim3s" : "extrant3";
        $is_fatal = 1;
      }
      #!# ERROR vadr:REVCOMPLEM:revcompl == flan:REVCOMPLEM:The input sequence is the reverse complementary strand of the coding sequence.
      elsif($error =~ /^The input sequence is the reverse complementary strand of the coding sequence./) { 
        $type_id  = "sequence";
        $do_push  = 1;
        $alt_code = "revcompl";
        $is_fatal = 1;
      }
      if(! defined $do_push) { 
        die "ERROR did not recognize ERROR line: $line\n";
      }
      if($do_push) { 
        if(! defined $type_id) { 
          die "ERROR type_id not set for pushable warning on line: $line\n";
        }
        if(! defined $alt_code) { 
          die "ERROR alt_code not set for pushable warning on line: $line\n";
        }
        if(! defined $flan_alt_info_HHAR->{$seqname}) { 
          %{$flan_alt_info_HHAR->{$seqname}} = ();
        }
        if(! defined $flan_alt_info_HHAR->{$seqname}{$type_id}) { 
          @{$flan_alt_info_HHAR->{$seqname}{$type_id}} = ();
        }

        # 02.09.24: do not allow fatal errors for PB1-F2, PA-X or NB
        if(($is_fatal) && 
           (($type_id eq "CDS.PA-X") || ($type_id eq "CDS.PB1-F2") || ($type_id eq "CDS.NB"))) { 
          $is_fatal = 0;
          if($do_weird) { 
            printf("WEIRD: $seqname has ERROR for $type_id, we are ignoring it because indexers misc_featurize it\n");
          }
        }

        my $code_fail_desc = sprintf($alt_code . "%s" . $error, ($is_fatal) ? "::fatal::" : "::nonfatal::");
        push(@{$flan_alt_info_HHAR->{$seqname}{$type_id}}, $code_fail_desc);

        if($is_fatal) { 
          print ERRFLAN "$seqname $type_id $alt_code $line\n";
        }
        else { 
          print WRNFLAN "$seqname $type_id $alt_code $line\n";
        }
      }
    }
  }
  # final check of frameshifts
  if(($frameshift_info ne "") && ($saw_frameshift_flag)) { 
    die "ERROR for $seqname, frameshift info non-empty but frameshift flag not raised";
  }
  close(INFLAN);
  close(ERRFLAN);
  close(WRNFLAN);

}

#################################################################
# Subroutine:  process_vadr_feature_tables()
# Incept:      EPN, Wed May  3 12:16:28 2023
#
# Purpose:    Parse and save info in a VADR feature table. Create a 
#             'processed' version of the file that sqf_FeatureTableParse() 
#             can parse and fill %vadr_ftr_info_HAH.
#
# Arguments: 
#  $vadr_pass_file:        the VADR feature table file with passing seqs
#  $vadr_fail_file:        the VADR feature table file with failing seqs
#  $vadr_seqname_HR:       reference to an 'exists' hash of the VADR sequence names
#  $vadr_ftr_info_HAHR:    ref to feature info 
#  $product2gene_HR:       ref to hash mapping products to genes
#  $do_misc2cds:           TRUE to treat misc_feature as CDS
#
# Returns:    void
#
################################################################# 
sub process_vadr_feature_tables { 
  my $sub_name = "process_vadr_feature_tables";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($vadr_pass_file, $vadr_fail_file, $vadr_seqname_HR, $vadr_ftr_info_HAHR, $product2gene_HR, $do_misc2cds) = @_;

  my $p_vadr_ft_file = $vadr_pass_file . ".p";
  $p_vadr_ft_file =~ s/\.pass//;

  open(OUTVADR, ">", $p_vadr_ft_file) || die "ERROR unable to open $p_vadr_ft_file for writing";

  my $line;
  my $seqname;

  for my $file ("pass", "fail") { 
    if($file eq "pass") { 
      open(INVADR, $vadr_pass_file) || die "ERROR unable to open $vadr_pass_file for reading";
    }
    else { 
      open(INVADR, $vadr_fail_file) || die "ERROR unable to open $vadr_fail_file for reading";
    }
    while($line = <INVADR>) { 
      chomp $line;
      if($line =~ /^\>Feature\s+(\S+)/) { 
        $seqname = shorten_seqname_no_version($1);
        $vadr_seqname_HR->{$seqname} = 1;
        print OUTVADR $line . "\n";
      }
      elsif($line =~ /^Additional/) { 
        ;
      }    
      elsif($line =~ /^ERROR/) { 
        ;
      }    
      else {
        print OUTVADR $line . "\n";
      }
    }
    close(INVADR);
  }

  close(OUTVADR);

  sqf_FeatureTableParse($p_vadr_ft_file, $vadr_ftr_info_HAHR, undef);
  truncate_coords($vadr_ftr_info_HAHR);
  #utl_HAHDump("vadr_ftr_info", $vadr_ftr_info_HAHR, *STDOUT);

  if($do_misc2cds) { 
    foreach my $seqname (sort keys %vadr_ftr_info_HAH) { 
      for(my $i = 0; $i < scalar(@{$vadr_ftr_info_HAHR->{$seqname}}); $i++) { 
        if($vadr_ftr_info_HAHR->{$seqname}[$i]{"type"} eq "misc_feature") { 
          if(defined $vadr_ftr_info_HAHR->{$seqname}[$i]{"note"}) { 
            if($vadr_ftr_info_HAHR->{$seqname}[$i]{"note"} =~ /^similar to (.+)$/) { 
              #printf("note: " . $vadr_ftr_info_HAHR->{$seqname}[$i]{"note"} . "\n");
              my $product = $1;
              $vadr_ftr_info_HAHR->{$seqname}[$i]{"type"}    = "CDS";
              $vadr_ftr_info_HAHR->{$seqname}[$i]{"product"} = $product;
              if(! defined $product2gene_HR->{$product}) { 
                die "ERROR in $sub_name, trying to change misc_feature to CDS but product $product not in product2gene map";
              }
              $vadr_ftr_info_HAHR->{$seqname}[$i]{"gene"} = $product2gene_HR->{$product};
              delete $vadr_ftr_info_HAHR->{$seqname}[$i]{"note"};
            }
            else { 
              die "ERROR unparseable note for misc_feature: " . $vadr_ftr_info_HAHR->{$seqname}[$i]{"note"};
            }
          }
          else { 
            die "ERROR no note for misc_feature";
          }
        }
      }
    }
  }
  
  return;
}

#################################################################
# Subroutine:  process_vadr_alt_file()
# Incept:      EPN, Wed May  3 13:03:06 2023
#
# Purpose:    Parse and save info in a VADR .alt file. 
#
# Arguments: 
#  $vadr_alt_file:         the VADR .alt file 
#  $out_root:              output root for output files
#  $vadr_alt_info_HHAR:    ref to 2d hash of arrays
#  $vadr_ftr_info_HAHR:    ref to feature info
#  $product2gene_HR:       ref to hash with map from product to gene values
#
# Returns:    void
#
################################################################# 
sub process_vadr_alt_file { 
  my $sub_name = "process_vadr_alt_file";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($vadr_alt_file, $out_root, $vadr_alt_info_HHAR, $vadr_ftr_info_HAHR, $product2gene_HR) = @_;

  my $vadr_fatal_file    = $out_root . ".vadr.fatal.list";
  my $vadr_nonfatal_file = $out_root . ".vadr.nonfatal.list";

  open(INALT,            $vadr_alt_file)      || die "ERROR unable to open $vadr_alt_file for reading";
  open(OUTFATAL,    ">", $vadr_fatal_file)    || die "ERROR unable to open $vadr_fatal_file for reading";
  open(OUTNONFATAL, ">", $vadr_nonfatal_file) || die "ERROR unable to open $vadr_nonfatal_file for reading";

  ##       seq                                     ftr   ftr             ftr  alert           alert                                 seq   seq          mdl   mdl  alert 
  ##idx    name                          model     type  name            idx  code      fail  description                        coords   len       coords   len  detail
  ##-----  ----------------------------  --------  ----  --------------  ---  --------  ----  -----------------------------  ----------  ----  -----------  ----  ------
  #4.1.1   gi|1928843018|gb|MW220253.1|  CY003646  CDS   polymerase_PB1    3  fsthicft  yes   POSSIBLE_FRAMESHIFT_HIGH_CONF   9..2121:+  2113  187..2298:+  2112  high confidence possible frameshift in CDS (frame not restored before end) [cause:insert,S:9(1),M:186; frame:<3(1); length:<8:(2113); shifted_avgpp:0.974; exp_avgpp:0.934;]
  while(my $line = <INALT>) { 
    chomp $line;
    if($line !~ m/^\#/) { 
      my @el_A = split(/\s+/, $line);
      my ($long_seqname, $ftr_type, $ftr_name, $alt_code, $fail, $desc, $detail) = ($el_A[1], $el_A[3], $el_A[4], $el_A[6], $el_A[7], $el_A[8], $el_A[13]);
      for(my $i = 14; $i < scalar(@el_A); $i++) { 
        $detail .= " " . $el_A[$i];
      }
      my $seqname = shorten_seqname_no_version($long_seqname);
      $ftr_name =~ s/\_/ /g;
      if($ftr_type eq "CDS") { 
        if(! defined $product2gene_HR->{$ftr_name}) { 
          die "ERROR in $sub_name, no gene info for product $ftr_name in the map";
        }
        $ftr_name = $product2gene_HR->{$ftr_name};
      }
      my $type_id = ($ftr_type eq "-") ? "sequence" : ($ftr_type . "." . $ftr_name);
      if(! defined $vadr_alt_info_HHAR->{$seqname}) { 
        %{$vadr_alt_info_HHAR->{$seqname}} = ();
      }
      if(! defined $vadr_alt_info_HHAR->{$seqname}{$type_id}) { 
        @{$vadr_alt_info_HHAR->{$seqname}{$type_id}} = ();
      }
      my $code_fail_detail = sprintf("%s::%s::%s", $alt_code, ($fail eq "yes" ? "fatal" : "nonfatal"), $detail);
      push(@{$vadr_alt_info_HHAR->{$seqname}{$type_id}}, $code_fail_detail);

      if($fail eq "yes") { 
        print OUTFATAL "$seqname $type_id $alt_code $desc $detail\n";
      }
      else { 
        print OUTNONFATAL "$seqname $type_id $alt_code $desc $detail\n";
      }
    }
  }
  close(OUTFATAL);
  close(OUTNONFATAL);

  return;
}

#################################################################
# Subroutine:  process_vadr_sqa_file()
# Incept:      EPN, Fri Aug 18 13:43:43 2023
#
# Purpose:    Parse and save info on group and subgroup in a VADR .sqa file. 
#
# Arguments: 
#  $vadr_sqa_file:         the VADR .sqa file 
#  $vadr_seq_info_HHR:      ref to 2d hash of hashes
#
# Returns:    void
#
################################################################# 
sub process_vadr_sqa_file { 
  my $sub_name = "process_vadr_sqa_file";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($vadr_sqa_file, $vadr_seq_info_HHR) = @_;

  open(INSQA, $vadr_sqa_file) || die "ERROR unable to open $vadr_sqa_file for reading";

  ##seq  seq                              seq             best                  sub                             seq   
  ##idx  name                             len  p/f   ant  model      grp        grp  nfa  nfn  nf5  nf3  nfalt  alerts
  ##---  ------------------------------  ----  ----  ---  ---------  ---------  ---  ---  ---  ---  ---  -----  ------
  #1     gi|211910015|ref|NC_006306.2|    935  PASS  yes  NC_006306  fluC-seg7  -      4    0    0    0      0  -     
  #9     gi|817124656|gb|CY190157.1|     1434  PASS  yes  CY002538   fluA-seg6  N1     2    0    0    0      0  -     
  while(my $line = <INSQA>) { 
    chomp $line;
    if($line !~ m/^\#/) { 
      my @el_A = split(/\s+/, $line);
      my ($long_seqname, $grp, $subgrp) = ($el_A[1], $el_A[6], $el_A[7]);
      my $seqname = shorten_seqname_no_version($long_seqname);
      if(defined $vadr_seq_info_HHR->{$seqname}) { 
        die "ERROR in $vadr_sqa_file, read two lines for $long_seqname";
      }
      if($grp =~ /flu(\S+)\-seg(\d)/) { 
        $vadr_seq_info_HHR->{$seqname}{"type"} = $1;
        $vadr_seq_info_HHR->{$seqname}{"segment"} = $2;
      }
      else { 
        die "ERROR unable to parse line in $vadr_sqa_file:\n$line\n";
      }
      if($subgrp ne "-") { 
        $vadr_seq_info_HHR->{$seqname}{"serotype"} = $subgrp;
      }
    }
  }
  close(INSQA);

  return;
}

#################################################################
# Subroutine:  get_output_for_warnings_errors_alerts()
# Incept:      EPN, Wed May  3 13:55:06 2023
#
# Purpose:    Return a string that gives output for warnings, errors and alerts for a sequence
#
# Arguments: 
#  $flan_ftr_HR:       ref to hash with FLAN feature info
#  $vadr_ftr_HR:       ref to hash with VADR feature info
#  $flan_alt_info_HAR: hash of array of FLAN warnings/errors for current sequence
#  $vadr_alt_info_HAR: hash of array of VADR alerts for current sequence
#  $product2gene_HR:   hash mapping product values to genes
#
# Returns:    void
#
################################################################# 
sub get_output_for_warnings_errors_alerts {
  my $sub_name = "get_output_for_warnings_errors_alerts";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($flan_ftr_HR, $vadr_ftr_HR, $flan_alt_info_HAR, $vadr_alt_info_HAR, $product2gene_HR) = @_;

  my $type = undef;
  my $id   = undef;

  my $flan_output  = "";
  my $flan_in_type = "";
  my $flan_in_id   = "";
  if(defined $flan_alt_info_HAR) { 
    if(defined $flan_ftr_HR) { 
      $flan_in_type = $flan_ftr_HR->{"type"};
      # per-feature alerts
      if($flan_in_type eq "CDS") { 
        if(! defined $product2gene_HR->{$flan_ftr_HR->{"product"}}) { 
          die "ERROR in $sub_name, no value for " . $flan_ftr_HR->{"product"} . " in $sub_name";
        }
        $flan_in_id = $product2gene_HR->{$flan_ftr_HR->{"product"}};
      }
      else { 
        $flan_in_id = $flan_ftr_HR->{"gene"};
      }
    }

    foreach my $key (sort keys (%{$flan_alt_info_HAR})) { 
      if($key =~ m/([^\.]+)\.(.+)$/) { 
        ($type, $id) = ($1, $2);
        if(($flan_in_type eq $type) && 
           ($flan_in_id   eq $id)) { 
          foreach my $alert (@{$flan_alt_info_HAR->{$key}}) { 
            $flan_output .= $key . "::" . $alert . ";;;";
          }
        }
      }
      elsif(($key eq "sequence") && (! defined $flan_ftr_HR)) { 
        foreach my $alert (@{$flan_alt_info_HAR->{$key}}) { 
          $flan_output .= $key . "::" . $alert . ";;;";
        }
      }
    }
  }

  my $vadr_output  = "";
  my $vadr_in_type = "";
  my $vadr_in_id   = "";
  if(defined $vadr_alt_info_HAR) { 
    if(defined $vadr_ftr_HR) { 
      # per feature alerts
      $vadr_in_type = $vadr_ftr_HR->{"type"};
      if($vadr_in_type eq "CDS") { # id is the gene, but we only have the product in $vadr_ftr_HR, so we use the map to get gene
        if(! defined $product2gene_HR->{$vadr_ftr_HR->{"product"}}) { 
          die "ERROR in $sub_name, no value for " . $vadr_ftr_HR->{"product"} . " in $sub_name";
        }
        $vadr_in_id = $product2gene_HR->{$vadr_ftr_HR->{"product"}};
      }
      elsif($vadr_in_type eq "mat_peptide") { # id is the product
        if(defined $vadr_ftr_HR->{"product"}) { 
          $vadr_in_id = $vadr_ftr_HR->{"product"};
        }
        else { 
         $vadr_in_id = "mat_peptide.1"; # the sig_peptide
        }
      }
      elsif($vadr_in_type eq "sig_peptide") { # id is the product
        if(defined $vadr_ftr_HR->{"product"}) { 
          $vadr_in_id = $vadr_ftr_HR->{"product"};
        }
        else { 
         $vadr_in_id = "sig_peptide.1"; # the sig_peptide
        }
      }
      else { 
        $vadr_in_id = $vadr_ftr_HR->{"gene"};
      }
      if(! defined $vadr_in_id) { 
        die "ERROR undefined vadr_in_id for type $type";
      }
    }

    foreach my $key (sort keys (%{$vadr_alt_info_HAR})) { 
      if($key =~ m/([^\.]+)\.(.+)$/) { 
        ($type, $id) = ($1, $2);
        if(($type eq $vadr_in_type) && 
           ($id   eq $vadr_in_id)) { 
          foreach my $alert (@{$vadr_alt_info_HAR->{$key}}) { 
            $vadr_output .= $key . "::" . $alert . ";;;";
          }
        }
      }
      elsif(($key eq "sequence") && (! defined $vadr_ftr_HR)) { 
        foreach my $alert (@{$vadr_alt_info_HAR->{$key}}) { 
          $vadr_output .= $key . "::" . $alert . ";;;";
        }
      }
    }
  }

  if($flan_output eq "") { $flan_output = "-"; }
  if($vadr_output eq "") { $vadr_output = "-"; }

  return($flan_output, $vadr_output);
}
  
#################################################################
# Subroutine:  shorten_seqname_no_version()
# Incept:      EPN, Wed May  3 14:09:03 2023
#
# Purpose:    Return shortened part of a long NCBI sequence name,
#             without a version:
#             Example:
#
#             in:     gi|1928843018|gb|MW220253.1| 
#             return: MW220253
#
# Arguments: 
#  $seqname:           potentially long sequence name
#
# Returns:    short seqname
#
################################################################# 
sub shorten_seqname_no_version {
  my $sub_name = "shorten_seqname_no_version";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($in_name) = (@_);

  if($in_name =~ /^gi\|\d+\|[^\|]+\|([^\|]+)\|$/) { 
    $in_name = $1;
  }
  elsif($in_name =~ /^gi\|\d+\|[^\|]+\|([^\|]+)\|[^\|]+$/) { 
    $in_name = $1;
  }
  if($in_name =~ /(\S+)\.\d+/) { 
    $in_name = $1; 
  }
  return $in_name;
}

#################################################################
# Subroutine:  truncate_coords()
# Incept:      EPN, Wed May  3 14:21:51 2023
#
# Purpose:    Add back in '<' and '>' to coordinates in a 
#             $ftr_info_HAR, based on 'trunc5' and 'trunc3' values.
#
# Arguments: 
#  $ftr_info_HAHR: feature info hash of arrays
#
# Returns:    void
#
################################################################# 
sub truncate_coords {
  my $sub_name = "truncate_coords";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAHR) = (@_);

  foreach my $seqname (sort keys %{$ftr_info_HAHR}) { 
    for(my $i = 0; $i < scalar(@{$ftr_info_HAHR->{$seqname}}); $i++) { 
      if((defined $ftr_info_HAHR->{$seqname}[$i]{"trunc5"}) && 
         $ftr_info_HAHR->{$seqname}[$i]{"trunc5"} == 1) { 
        $ftr_info_HAHR->{$seqname}[$i]{"coords"} = "<" . $ftr_info_HAHR->{$seqname}[$i]{"coords"};
      }
      if((defined $ftr_info_HAHR->{$seqname}[$i]{"trunc3"}) && 
         $ftr_info_HAHR->{$seqname}[$i]{"trunc3"} == 1) { 
        if($ftr_info_HAHR->{$seqname}[$i]{"coords"} =~ /^(.+\.\.)(\d+\:[\+\-])$/) { 
          $ftr_info_HAHR->{$seqname}[$i]{"coords"} = $1 . ">" . $2;
          #printf("ftr_info_HAHR->{$seqname}[$i]{coords} = " . $ftr_info_HAHR->{$seqname}[$i]{"coords"});
        }
        else { 
          die "ERROR unable to parse coords: " . $ftr_info_HAHR->{$seqname}[$i]{"coords"};
        }
      }
    }
  }
  return;
}
  
#################################################################
# Subroutine:  check_flan_trunc()
# Incept:      EPN, Thu May  4 10:49:45 2023
#
# Purpose:    Check if a feature is truncated or not given its
#             <type> and <id> (for <type> = CDS or mat_peptide this 'product', 
#             for gene this is 'gene').
#
# Arguments: 
#  $ftr_info_HAHR: feature info hash of arrays
#  $seqname:       sequence name
#  $type:          type
#  $gene:          gene value
#
# Returns:    Two values:
#             value 1: '1' if this feature is 5' truncated, else '0'
#             value 2: '1' if this feature is 3' truncated, else '0'
#             Dies if feature is not found or two matching features are found
#
################################################################# 
sub check_flan_trunc {
  my $sub_name = "check_flan_trunc";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAHR, $seqname, $type, $gene) = (@_);

  if(($type ne "CDS") && ($type ne "mat_peptide") && ($type ne "gene")) { 
    die "ERROR in $sub_name, unexpected type: $type"; 
  }

  my $found_match = 0;
  my $ret_val1 = undef;
  my $ret_val2 = undef;

  for(my $i = 0; $i < scalar(@{$ftr_info_HAHR->{$seqname}}); $i++) { 
    if($ftr_info_HAHR->{$seqname}[$i]{"type"} eq $type) {
      if(! defined $ftr_info_HAHR->{$seqname}[$i]{"gene"}) { 
        die "ERROR in $sub_name type: $type seqname: $seqname no gene\n";
      }
      if((($type eq "CDS")         && ($ftr_info_HAHR->{$seqname}[$i]{"gene"} eq $gene)) || 
         (($type eq "mat_peptide") && ($ftr_info_HAHR->{$seqname}[$i]{"gene"} eq $gene)) || 
         (($type eq "gene")        && ($ftr_info_HAHR->{$seqname}[$i]{"gene"} eq $gene))) { 
        if($found_match) { 
          die "ERROR in $sub_name, found two matches (seqname: $seqname type: $type gene: $gene)";
        }
        $found_match = 1;
        $ret_val1 = ((defined $ftr_info_HAHR->{$seqname}[$i]{"trunc5"}) && ($ftr_info_HAHR->{$seqname}[$i]{"trunc5"} == 1)) ? 1 : 0;
        $ret_val2 = ((defined $ftr_info_HAHR->{$seqname}[$i]{"trunc3"}) && ($ftr_info_HAHR->{$seqname}[$i]{"trunc3"} == 1)) ? 1 : 0;
      }
    }
  }
  if(! $found_match) { 
    die "ERROR found no match for $seqname $type $gene in $sub_name";
  }
  return ($ret_val1, $ret_val2);
}
  

#################################################################
# Subroutine:  process_product_gene_map_file()
# Incept:      EPN, Thu May  4 13:29:25 2023
#
# Purpose:    Parse the product:gene map file and fill the 
#             product2gene_H hash.
#
# Arguments: 
#  $product_gene_map_file: the input file
#  $product2gene_HR:       ref to hash to fill here
#
# Returns:    void
#             Dies if same product occurs more than once (a gene can occur more than once)
#
################################################################# 
sub process_product_gene_map_file { 
  my $sub_name = "process_product_gene_map_file";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($product_gene_map_file, $product2gene_HR, $gene2product_HR) = (@_);

  %{$product2gene_HR} = ();
  open(IN, $product_gene_map_file) || die "ERROR unable to open $product_gene_map_file for reading";
  while(my $line = <IN>) { 
    #product:BM2 protein gene:BM2
    chomp $line;
    if($line =~ m/^product\:(.+)\s+gene\:(.+)$/) { 
      my ($product, $gene) = ($1, $2);
      if(defined $product2gene_HR->{$product}) { 
        die "ERROR read product $product twice in $sub_name";
      }
      if(defined $gene2product_HR->{$gene}) { 
        # die "ERROR read gene $gene twice in $sub_name";
        # we know this happens, so we put in a flag and if we ever try to use this gene to determine a product, we'll exit
        $gene2product_HR->{$gene} = "FAIL";
      }
      else { 
        $gene2product_HR->{$gene} = $product;
      }
      $product2gene_HR->{$product} = $gene;
      #print("product2gene_HR->{$product} $product2gene_HR->{$product}\n"); 
    }
  }
  close(IN);

  return;
}

#################################################################
# Subroutine:  get_coords_compare_string()
# Incept:      EPN, Wed May 10 11:53:17 2023
#
# Purpose:    Compare flan and vadr coords and return either:
#             "coords-same"  if identical
#             "coords-same*" if identical except for truncation mode (same positions)
#             "coords-diff"  if different positions
#
# Arguments: 
#  $flan_coords:     flan coords string
#  $vadr_coords:     vadr coords string
#
# Returns:    See purpose
#
################################################################# 
sub get_coords_compare_string {
  my $sub_name = "get_coords_compare_string";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($flan_coords, $vadr_coords) = (@_);

  if($flan_coords eq $vadr_coords) { 
    return "coords-same";
  }
  # remove truncation marks
  $flan_coords =~ s/\>//g;
  $flan_coords =~ s/\<//g;
  $vadr_coords =~ s/\>//g;
  $vadr_coords =~ s/\<//g;
  if($flan_coords eq $vadr_coords) { 
    return "coords-same*";
  }
  return "coords-diff";
}

#################################################################
# Subroutine:  find_flan_cds_gene_given_idx()
# Incept:      EPN, Wed May 10 12:48:17 2023
#
# Purpose:    Return the 'gene' value for CDS index <cds_idx> in 
#             @{$flan_ftr_AHR}.
#
# Arguments: 
#  $flan_ftr_info_AHR: array of hash with flan ftr info
#  $cds_idx:           index of cds we're interested in
#
# Returns:    The 'gene' value if we find the CDS, else
#             if we don't find the CDS, return "",
#             this happens for HE584760.1.
#
################################################################# 
sub find_flan_cds_gene_given_idx { 
  my $sub_name = "find_flan_cds_gene_given_idx";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($flan_ftr_info_AHR, $cds_idx) = (@_);

  if(! defined $flan_ftr_info_AHR) { 
    die "ERROR in $sub_name, flan_ftr_info_AHR is undefined";
  }

  my $cur_cds_idx = 0;
  for(my $i = 0; $i < scalar(@{$flan_ftr_info_AHR}); $i++) { 
    if($flan_ftr_info_AHR->[$i]{"type"} eq "CDS") { 
      $cur_cds_idx++;
      if($cur_cds_idx == $cds_idx) { 
        if(! defined $flan_ftr_info_AHR->[$i]{"gene"}) { 
          die "ERROR in $sub_name, found cds $cds_idx but it does not have a gene value";
        }
        return $flan_ftr_info_AHR->[$i]{"gene"};
      }
    }
  }
  return "";
}

#################################################################
# Subroutine:  get_feature_alert_compare_string()
# Incept:      EPN, Thu May 11 11:59:32 2023
#
# Purpose:    Compare alert strings for the same feature from vadr and flan 
#             and classify them as:
#             "same:none:nonfatal":          if identical and empty for both
#             "same:consistent:nonfatal":    if not empty and all non-fatal for both and at least 1 consistent
#             "same:consistent:fatal":       if not empty and all     fatal for both and at least 1 consistent
#             "diff:vadronly:nonfatal":      flan empty, vadr non-empty and non-fatal
#             "diff:vadronly:fatal":         flan empty, vadr non-empty and at least 1 fatal
#             "diff:flanonly:nonfatal":      vadr empty, flan non-empty and non-fatal
#             "diff:flanonly:fatal":         vadr empty, flan non-empty and at least 1 fatal
#             "diff:inconsistent:nonfatal":  both non empty, all non fatal, none are consistent
#             "diff:inconsistent:fatal":     both non empty, both at least 1 fatal, none are consistent
#             "diff:inconsistent:flanfatal": both non empty, vadr all non fatal, flan >= 1 fatal, none are consistent
#             "diff:inconsistent:vadrfatal": both non empty, vadr >= 1 fatal, flan all non fatal, none are consistent
#             "diff:consistent:flanfatal":   both non empty, vadr all non fatal, flan >= 1 fatal, at least 1 consistent
#             "diff:consistent:vadrfatal":   both non empty, vadr >= 1 fatal, flan all non fatal, at least 1 consistent
#
# Arguments: 
#  $flan_alts:   flan alert string
#  $vadr_alts:   vadr alert string
#
# Returns:    See purpose
#
################################################################# 
sub get_feature_alert_compare_string { 
  my $sub_name = "get_feature_alert_compare_string";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($flan_alts, $vadr_alts) = (@_);

  #printf("$sub_name flan_alts: $flan_alts vadr_alts: $vadr_alts\n");

  if(($flan_alts eq "-") && ($vadr_alts eq "-")) { 
    return "same:none:nonfatal";
  }
  else { 
    my @flan_alt_A = ();
    my @vadr_alt_A = ();
    my $flan_has_fatal = ($flan_alts =~ m/\:\:fatal\:\:/) ? 1 : 0; 
    my $vadr_has_fatal = ($vadr_alts =~ m/\:\:fatal\:\:/) ? 1 : 0; 
    if($flan_alts ne "-") { 
      @flan_alt_A = split(";;;", $flan_alts);
    }
    if($vadr_alts ne "-") { 
      @vadr_alt_A = split(";;;", $vadr_alts);
    }
    my $nflan_alt = scalar(@flan_alt_A);
    my $nvadr_alt = scalar(@vadr_alt_A);
    if($nflan_alt == 0) { 
      if($nvadr_alt == 0) { 
        die "ERROR in $sub_name, coding error 1";
      }
      if($vadr_has_fatal) { 
        return "diff:vadronly:fatal";
      }
      else { 
        return "diff:vadronly:nonfatal";
      }
    }
    if($nvadr_alt == 0) { 
      if($nflan_alt == 0) { 
        die "ERROR in $sub_name, coding error 2";
      }
      if($vadr_has_fatal) { 
        return "diff:flanonly:fatal";
      }
      else { 
        return "diff:flanonly:nonfatal";
      }
    }
    # if we get here we have both vadr and flan alerts
    # check if they are all consistent
    # to be consistent: 
    # - at least one flan alert must be consistent with 1 vadr alert
    my $consistent = 0;
    my $flan_alt;
    my $vadr_alt;
    #########################################################################
    foreach my $flan_alt (@flan_alt_A) { 
      my $found_consistent = 0;
      my ($flan_type_id, $flan_code, $flan_fatal, $flan_msg) = (undef, undef, undef);
      if($flan_alt =~ /^([^\:]+)\:\:([^\:]+)\:\:([^\:]+)\:\:(.+)$/) { 
        ($flan_type_id, $flan_code, $flan_fatal, $flan_msg) = ($1, $2, $3, $4);
      }
      else { 
        die "ERROR in $sub_name, unable to parse FLAN alert: $flan_alt\n$flan_alts";
      }
      #printf("consistent: $consistent checking $flan_alt flan_code: $flan_code\n");
      foreach my $vadr_alt (@vadr_alt_A) { 
        my ($vadr_type_id, $vadr_code, $vadr_fatal, $vadr_msg) = (undef, undef, undef);
        if($vadr_alt =~ /^([^\:]+)\:\:([^\:]+)\:\:([^\:]+)\:\:(.+)$/) { 
          ($vadr_type_id, $vadr_code, $vadr_fatal, $vadr_msg) = ($1, $2, $3, $4);
        }
        else { 
          die "ERROR in $sub_name, unable to parse VADR alert: $vadr_alt";
        }
        my @flan_code_A = split(",", $flan_code);
        foreach my $indi_flan_code (@flan_code_A) { 
          my @vadr_code_A = split("," , $vadr_code);
          foreach my $indi_vadr_code (@vadr_code_A) { 
            #printf("\tchecking indi_flan_code vs indi_vadr_code $indi_flan_code vs $indi_vadr_code\n");
            if(($indi_flan_code =~ /$indi_vadr_code/) || ($indi_vadr_code =~ m/$indi_flan_code/)) { 
              $found_consistent = 1;
              #printf("HEYA found consistency comparing $flan_alt ($flan_code) and $vadr_alt ($vadr_code)\n");
            }
          }
        }
      }
      if($found_consistent) { 
        $consistent = 1;
      }
    }
    #########################################################################
    foreach my $vadr_alt (@vadr_alt_A) { 
      my $found_consistent = 0;
      my ($vadr_type_id, $vadr_code, $vadr_fatal, $vadr_msg) = (undef, undef, undef);
      if($vadr_alt =~ /^([^\:]+)\:\:([^\:]+)\:\:([^\:]+)\:\:(.+)$/) { 
        ($vadr_type_id, $vadr_code, $vadr_fatal, $vadr_msg) = ($1, $2, $3, $4);
      }
      else { 
        die "ERROR in $sub_name, unable to parse FLAN alert: $vadr_alt";
      }
      foreach my $flan_alt (@flan_alt_A) { 
        my ($flan_type_id, $flan_code, $flan_fatal, $flan_msg) = (undef, undef, undef);
        if($flan_alt =~ /^([^\:]+)\:\:([^\:]+)\:\:([^\:]+)\:\:(.+)$/) { 
          ($flan_type_id, $flan_code, $flan_fatal, $flan_msg) = ($1, $2, $3, $4);
        }
        else { 
          die "ERROR in $sub_name, unable to parse VADR alert: $flan_alt";
        }
        if(($vadr_code =~ /$flan_code/) || ($flan_code =~ m/$vadr_code/)) { 
          $found_consistent = 1;
          #printf("HEYA found consistency comparing $vadr_alt ($vadr_code) and $flan_alt ($flan_code)\n");
        }
      }
      if($found_consistent) { 
        $consistent = 1;
      }
    }
    ##########################################################################
    if($consistent == 1) { 
      if(($flan_has_fatal) && ($vadr_has_fatal)) { 
        return "same:consistent:fatal";
      }
      elsif((! $flan_has_fatal) && (! $vadr_has_fatal)) { 
        return "same:consistent:nonfatal";
      }
      elsif(($flan_has_fatal) && (! $vadr_has_fatal)) { 
        return "diff:consistent:flanfatal";
      }
      elsif((! $flan_has_fatal) && ($vadr_has_fatal)) { 
        return "diff:consistent:vadrfatal";
      }
      else { 
        die "ERROR in $sub_name, consistent alerts but fatalness differs:\nflan_alts: $flan_alts\nvadr_alts: $vadr_alts\n";
      }
    }
    else { # not consistent 
      if(($flan_has_fatal) && ($vadr_has_fatal)) { 
        return "diff:inconsistent:fatal";
      }
      elsif((! $flan_has_fatal) && (! $vadr_has_fatal)) { 
        return "diff:inconsistent:nonfatal";
      }
      elsif(($flan_has_fatal) && (! $vadr_has_fatal)) { 
        return "diff:inconsistent:flanfatal";
      }
      elsif((! $flan_has_fatal) && ($vadr_has_fatal)) { 
        return "diff:inconsistent:vadrfatal";
      }
      else {
        die "ERROR in $sub_name, inconsistent alerts but fatalness differs:\nflan_alts: $flan_alts\nvadr_alts: $vadr_alts\n";
      }
    }
    die "ERROR in $sub_name, unhandled case:\nflan_alts: $flan_alts\nvadr_alts: $vadr_alts\n";
  }
}
