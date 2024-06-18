#!/usr/bin/env perl
# 
# table-models.pl: gather data and create the latex file with data on FLAN models
#                        
# EPN, Tue Jun 27 14:40:17 2023
# 
#
use strict;
use warnings;
use Getopt::Long;

my $in_data  = ""; # name of input data file

my $usage = "perl table-models.pl <FLAN vadr .minfo file> <blastdb directory (FLAN)> <blastdb directory (FLAN-plus)>\n";

my %fluC_old2new_H = ();
my %fluC_new2old_H = ();
$fluC_old2new_H{"GN364866"}  = "NC_006307";
$fluC_old2new_H{"GM968019"}  = "NC_006308";
$fluC_old2new_H{"GM968018"}  = "NC_006309";
$fluC_old2new_H{"GM968017"}  = "NC_006310";
$fluC_old2new_H{"NC_006311"} = "NC_006311";
$fluC_old2new_H{"GM968015"}  = "NC_006312";
$fluC_old2new_H{"NC_006306"} = "NC_006306";

$fluC_new2old_H{"NC_006307"} = "GN364866";
$fluC_new2old_H{"NC_006308"} = "GM968019";
$fluC_new2old_H{"NC_006309"} = "GM968018";
$fluC_new2old_H{"NC_006310"} = "GM968017";
$fluC_new2old_H{"NC_006311"} = "NC_006311";
$fluC_new2old_H{"NC_006312"} = "GM968015";
$fluC_new2old_H{"NC_006306"} = "NC_006306";

my %new_models_H = ();
$new_models_H{"CY103881"} = 1;
$new_models_H{"CY103882"} = 1;
$new_models_H{"CY103883"} = 1;
$new_models_H{"CY103884"} = 1;
$new_models_H{"CY103885"} = 1;
$new_models_H{"CY103886"} = 1;
$new_models_H{"CY103887"} = 1;
$new_models_H{"CY103888"} = 1;
$new_models_H{"CY125942"} = 1;
$new_models_H{"CY125943"} = 1;
$new_models_H{"CY125944"} = 1;
$new_models_H{"CY125945"} = 1;
$new_models_H{"CY125946"} = 1;
$new_models_H{"CY125947"} = 1;
$new_models_H{"CY125948"} = 1;
$new_models_H{"CY125949"} = 1;
$new_models_H{"ON637239"} = 1;
$new_models_H{"NC_036615"} = 1;
$new_models_H{"NC_036616"} = 1;
$new_models_H{"NC_036617"} = 1;
$new_models_H{"NC_036618"} = 1;
$new_models_H{"NC_036619"} = 1;
$new_models_H{"NC_036620"} = 1;
$new_models_H{"NC_036621"} = 1;
$new_models_H{"NC_036622"} = 1;

# there are 7 model accessions for which the 'version' used was .2, instead of .1
my %version2_H = ();
$version2_H{"CY000449"} = 1;
$version2_H{"NC_006307"} = 1;
$version2_H{"NC_006308"} = 1;
$version2_H{"NC_006309"} = 1;
$version2_H{"NC_006310"} = 1;
$version2_H{"NC_006312"} = 1;
$version2_H{"NC_006306"} = 1;

my @group_ordered1_A = ();
my @group_ordered2_A = ();
push(@group_ordered1_A, "fluA-seg1:1");
push(@group_ordered1_A, "fluA-seg1:2");
push(@group_ordered1_A, "fluA-seg1:3");
push(@group_ordered1_A, "fluA-seg2:1");
push(@group_ordered1_A, "fluA-seg2:2");
push(@group_ordered1_A, "fluA-seg2:3");
push(@group_ordered1_A, "fluA-seg3:1");
push(@group_ordered1_A, "fluA-seg3:2");
push(@group_ordered1_A, "fluA-seg3:3");
push(@group_ordered1_A, "fluA-seg4.H1:1");
push(@group_ordered1_A, "fluA-seg4.H2:1");
push(@group_ordered1_A, "fluA-seg4.H3:1");
push(@group_ordered1_A, "fluA-seg4.H4:1");
push(@group_ordered1_A, "fluA-seg4.H5:1");
push(@group_ordered1_A, "fluA-seg4.H6:1");
push(@group_ordered1_A, "fluA-seg4.H7:1");
push(@group_ordered1_A, "fluA-seg4.H8:1");
push(@group_ordered1_A, "fluA-seg4.H9:1");
push(@group_ordered1_A, "fluA-seg4.H10:1");
push(@group_ordered1_A, "fluA-seg4.H11:1");
push(@group_ordered1_A, "fluA-seg4.H12:1");
push(@group_ordered1_A, "fluA-seg4.H13:1");
push(@group_ordered1_A, "fluA-seg4.H14:1");
push(@group_ordered1_A, "fluA-seg4.H15:1");
push(@group_ordered1_A, "fluA-seg4.H16:1");
push(@group_ordered1_A, "fluA-seg4.H17:1");
push(@group_ordered1_A, "fluA-seg4.H18:1");
push(@group_ordered1_A, "fluA-seg4.H19:1");
push(@group_ordered1_A, "fluA-seg5:1");
push(@group_ordered1_A, "fluA-seg5:2");
push(@group_ordered1_A, "fluA-seg5:3");
push(@group_ordered1_A, "fluA-seg6.N1:1");
push(@group_ordered1_A, "fluA-seg6.N2:1");
push(@group_ordered1_A, "fluA-seg6.N3:1");
push(@group_ordered1_A, "fluA-seg6.N4:1");
push(@group_ordered1_A, "fluA-seg6.N5:1");
push(@group_ordered1_A, "fluA-seg6.N6:1");
push(@group_ordered1_A, "fluA-seg6.N7:1");
push(@group_ordered1_A, "fluA-seg6.N8:1");
push(@group_ordered1_A, "fluA-seg6.N9:1");
push(@group_ordered1_A, "fluA-seg6.N10:1");
push(@group_ordered1_A, "fluA-seg6.N11:1");
push(@group_ordered1_A, "fluA-seg7:1");
push(@group_ordered1_A, "fluA-seg7:2");
push(@group_ordered1_A, "fluA-seg7:3");
push(@group_ordered1_A, "fluA-seg8:1");
push(@group_ordered1_A, "fluA-seg8:2");
push(@group_ordered1_A, "fluA-seg8:3");

push(@group_ordered2_A, "fluB-seg1:1");
push(@group_ordered2_A, "fluB-seg2:1");
push(@group_ordered2_A, "fluB-seg3:1");
push(@group_ordered2_A, "fluB-seg4:1");
push(@group_ordered2_A, "fluB-seg5:1");
push(@group_ordered2_A, "fluB-seg6:1");
push(@group_ordered2_A, "fluB-seg7:1");
push(@group_ordered2_A, "fluB-seg8:1");
push(@group_ordered2_A, "fluC-seg1:1");
push(@group_ordered2_A, "fluC-seg2:1");
push(@group_ordered2_A, "fluC-seg3:1");
push(@group_ordered2_A, "fluC-seg4:1");
push(@group_ordered2_A, "fluC-seg5:1");
push(@group_ordered2_A, "fluC-seg6:1");
push(@group_ordered2_A, "fluC-seg7:1");
push(@group_ordered2_A, "fluD-seg1:1");
push(@group_ordered2_A, "fluD-seg2:1");
push(@group_ordered2_A, "fluD-seg3:1");
push(@group_ordered2_A, "fluD-seg4:1");
push(@group_ordered2_A, "fluD-seg5:1");
push(@group_ordered2_A, "fluD-seg6:1");
push(@group_ordered2_A, "fluD-seg7:1");

if(scalar(@ARGV) != 3) { die $usage; }
my ($minfo_file, $blast_dir, $blast_plus_dir) = @ARGV;
my ($ftr_mdl_name, $type, $coords, $nexons);

my $blastdb_file;
my $blastdb_plus_file;
my ($mdl_name, $blastdb, $group, $length, $subgroup, $product);

# 1D hashes with info per group, key is group
# group is actual a concatenation of group and subgroup (if any) from MODEL lines
# (example of group: 'fluA-seg1' or 'fluA-seg4.H14')
my %model_H;      # model name for this group (can only be one, script checks for this and dies if more than one)            
my %length_H;     # length for model for this group

# 1D hashes with info per model, key is model:
my %ncds_H = ();           # number of CDS for this model
my %mdl_name2group_H = (); # group name for this model

# 1D hashes with info per product, key is product:
my %product2gene_H = ();   # gene that corresponds to this product (can only be one, script checks for this and dies if more than one)

my @group_A = (); # array of all groups

# 2D hashes that contain info on each CDS, 1D key: $mdl_name, 2D key: $product
my %nexons_HH;                # number of exons for this CDS
my %ncoords_HH;               # number of distinct start/stops (coords) in non-plus minfo for this CDS
my %nblastseqs_HH = ();       # number of blast protein seqs for this CDS, in first blastdb dir
my %nblastseqs_plus_HH = ();  # number of blast protein seqs for this CDS, in second blastdb dir
my %nmp_HH = ();              # number of mat_peptides for this CDS
my %mnf_HH = ();              # '1' if this feature has misc_not_failure attribute set, 0 if not

my $ftr_idx = 0;
my $cur_cds_ftr_idx = 0;

open(IN, $minfo_file) || die "ERROR unable to open $minfo_file for reading";
#my $fluD_flag = 0; # set to 1 when we're reading a fluD model
my %seen_group_H = ();
while(my $line = <IN>) { 
  # for model lines, store group (e.g. fluB-seg8), length, and model accession
  if($line =~ m/^MODEL/) { 
    #if($line =~ m/fluD/) { 
    #  $fluD_flag = 1;
    #}
    #else { 
    #  $fluD_flag = 0;
    #}
    #$fluD_flag = 0;
    #if(! $fluD_flag) {
    $ftr_idx = 0;
    chomp $line;
    #MODEL AY504614 blastdb:"AY504614.vadr.protein.fa" group:"fluB-seg8" length:"1097"
    my @el_A = split(/\s+/, $line);
    ($mdl_name, $blastdb, $group, $length) = ($el_A[1], $el_A[2], $el_A[3], $el_A[4]);
    #printf("read MODEL line:\n$line\n");
    if($blastdb =~ /^blastdb\:\"([^\"]+)\"$/) { 
      $blastdb = $1;
    }
    else { 
      die "ERROR unable to parse blastdb out of MODEL line: $line";
    }
    my $blastdb_maybe_swapped = $blastdb;
    if(defined $fluC_new2old_H{$mdl_name}) { 
      $blastdb_maybe_swapped = $fluC_new2old_H{$mdl_name} . ".vadr.protein.fa";
    }        
    #printf("blastdb               $blastdb\n");
    #printf("blastdb_maybe_swapped $blastdb_maybe_swapped\n");
    if($group =~ /^group\:\"([^\"]+)\"$/) { 
      $group = $1;
    }
    else { 
      die "ERROR unable to parse group ($group) out of MODEL line: $line";
    }
    if($length =~ /^length\:\"([^\"]+)\"$/) { 
      $length = $1;
    }
    else { 
      die "ERROR unable to parse length out of MODEL line: $line";
    }
    $subgroup = undef;
    if($line =~ /subgroup\:\"([^\"]+)\"/) { 
      $subgroup = $1;
    }
    if(defined $subgroup) { 
      $group .= "." . $subgroup;
    }

    my $gidx = 1;
    while(defined $seen_group_H{"$group:$gidx"}) { 
      $gidx++;
    }
    $group .= ":" . $gidx;
    $seen_group_H{$group} = 1;

    %{$ncoords_HH{$mdl_name}} = ();

    if(defined $model_H{$group}) { 
      die "ERROR read two MODEL lines for same group: $group";
    }
    $model_H{$group} = $mdl_name;
    $mdl_name2group_H{$mdl_name} = $group;
    $length_H{$group} = $length; 
    # determine number of proteins in blastdb
    $blastdb_file      = $blast_dir      . "/" . $blastdb_maybe_swapped;
    $blastdb_plus_file = $blast_plus_dir . "/" . $blastdb;
    if(! -s $blastdb_file) { 
      if(defined $new_models_H{$mdl_name}) { 
        $blastdb_file = undef;
      }
      else { 
        die "ERROR blastdb_file $blastdb_file does not exist or is empty";
      }
    }
    if(! -s $blastdb_plus_file) { 
      die "ERROR blastdb_plus_file $blastdb_plus_file does not exist or is empty";
    }
  }
########################################
# for feature CDS lines, store
# - product
# - gene
# - misc_not_failure?
# - number of exons
# - number of proteins in the two blast dbs
# - ftr_idx
#elsif(($line =~ m/^FEATURE.+type:\"CDS\"/) && (! $fluD_flag)) { 
  elsif($line =~ m/^FEATURE.+type:\"CDS\"/) { 
    chomp $line;
    #FEATURE NC_006311 type:"CDS" coords:"30..1727:+" parent_idx_str:"GBNULL" gene:"NP" product:"nucleoprotein"
    my @el_A = split(/\s+/, $line);
    ($ftr_mdl_name, $type, $coords) = ($el_A[1], $el_A[2], $el_A[3]);
    if($type =~ /^type\:\"([^\"]+)\"$/) { 
      $type = $1;
    }
    else { 
      die "ERROR unable to parse type out of FEATURE line: $line";
    }
    if($coords =~ /^coords\:\"([^\"]+)\"$/) { 
      $coords = $1;
    }
    else { 
      die "ERROR unable to parse coords out of FEATURE line: $line";
    }
    # determine number of exons
    if($coords =~ /^\d+\.\.\d+\:\+$/) { 
      $nexons = 1;
    }
    elsif($coords =~ /^\d+\.\.\d+\:\+\,\d+\.\.\d+\:\+$/) { 
      $nexons = 2;
    }
    else { 
      die "ERROR unable to determine number of exons in coords string $coords on line: $line";
    }

    if($ftr_mdl_name ne $mdl_name) { 
      die "ERROR unexpected model name for CDS, expected $mdl_name but read $ftr_mdl_name"; 
    }
    if(! defined $blastdb_file) { 
      if(! defined $new_models_H{$mdl_name}) { 
        die "ERROR blastdb_file not defined for (non-new model) line: $line\n";
      }
    }
    if(! defined $blastdb_plus_file) { 
      die "ERROR blastdb_plus_file not defined for line: $line\n";
    }
    my $nblastseqs      = (defined $blastdb_file) ? count_matching_coords_in_fasta($blastdb_file, $coords) : -1;
    my $nblastseqs_plus = count_matching_coords_in_fasta($blastdb_plus_file, $coords);

    $product = undef;
    my $gene    = undef;
    my $mnf     = 0;
    if($line =~ /product\:\"([^\"]+)\"/) { 
      $product = $1;
    }
    else { 
      die "ERROR unable to parse product out of FEATURE line: $line";
    }
    if($line =~ /gene\:\"([^\"]+)\"/) { 
      $gene = $1;
    }
    else { 
      die "ERROR unable to parse gene out of FEATURE line: $line";
    }
    if((defined $product2gene_H{$product}) && 
       ($product2gene_H{$product} ne $gene)) { 
      if(($ftr_mdl_name eq "NC_036621") && ($gene eq "NS2")) { 
        # we know about this CDS with product='nonstructural protein 2' has gene as NS2, instead of NEP, we deal with this manually
      }        
      else { 
        die "ERROR two different gene values read for product $product";
      }
    }
    if(($ftr_mdl_name ne "NC_036621") || ($gene ne "NS2")) { 
      $product2gene_H{$product} = $gene;
    }

    # check if there's a misc_not_failure:
    if($line =~ /misc_not_failure\:\"1\"/) { 
      $mnf = 1;
    }
    
    if(! defined $ncoords_HH{$mdl_name}{$product}) { 
      $ncoords_HH{$mdl_name}{$product} = 1; 
      $nblastseqs_HH{$mdl_name}{$product} = $nblastseqs;
      $nblastseqs_plus_HH{$mdl_name}{$product} = $nblastseqs_plus;
      $nexons_HH{$mdl_name}{$product} = $nexons;
      $mnf_HH{$mdl_name}{$product} = $mnf;
      
      if(! defined $ncds_H{$mdl_name}) { 
        $ncds_H{$mdl_name} = 1;
      }
      else { 
        $ncds_H{$mdl_name}++;
      }
    }
    else { 
      $nexons_HH{$mdl_name}{$product} = $nexons;
      $ncoords_HH{$mdl_name}{$product}++;
      $nblastseqs_HH{$mdl_name}{$product} += $nblastseqs;
      $nblastseqs_plus_HH{$mdl_name}{$product} += $nblastseqs_plus;
    }
    $cur_cds_ftr_idx = $ftr_idx;
    $ftr_idx++;
  } # end of if FEATURE.*CDS
#elsif(($line =~ m/^FEATURE.*mat_peptide/) && (! $fluD_flag)) { 
  elsif($line =~ m/^FEATURE.*mat_peptide/) { 
    my $parent_ftr_idx = undef;
    if($line =~ m/^FEATURE.*parent_idx_str\:\"(\d+)\"/) { 
      $parent_ftr_idx = $1;
      if($parent_ftr_idx ne $cur_cds_ftr_idx) { 
        die "ERROR unable to match up mat_peptide with CDS";
      }
      if(! defined $nmp_HH{$mdl_name}{$product}) { 
        $nmp_HH{$mdl_name}{$product} = 1;
      }
      else {
        $nmp_HH{$mdl_name}{$product}++;
      }
    }
    else {
      die "ERROR unable to find parent_ftr_idx in mat_peptide line $line";
    }
  }
  elsif($line =~ m/^FEATURE/) { 
    $ftr_idx++;
  }
}

my $caption1 = "\\textbf{List of VADR and FLAN influenza A model reference sequences and attributes of associated proteins.} The ``\\#coords\'' column indicates the number of distinct pairs of start and stop genome nucleotide coordinates for all proteins in the set. The ``\\# proteins'' column indicates the number of proteins in the set, and the ``\\#extra'' column indicates the number of additional proteins added to the VADR protein set not present in the FLAN set based on analysis of training set results. For CDS that have italicized ``CDS product'' and ``gene'' names, FLAN errors are converted to warnings and VADR converts them to ``misc\\_feature'' features if they have certain usually fatal alerts, instead of failing the sequence. Bold model accessions indicate models without an analog in FLAN.";

my $caption2 = "\\textbf{List of VADR and FLAN influenza B, C, and D model reference sequences and attributes of associated proteins.} The ``\\#coords'' column indicates the number of distinct pairs of start and stop genome nucleotide coordinates for all proteins in the set. The ``\\# proteins'' column indicates the number of proteins in the set. For CDS that have italicized ``CDS product'' and ``gene'' names, FLAN errors are converted to warnings and VADR converts them to ``misc\\_feature'' features if they have certain usually fatal alerts, instead of failing the sequence. Bold model accessions indicate models without an analog in FLAN. Italicized model accessions indicate models for which VADR uses a different accession than FLAN.";

for(my $t = 0; $t < 2; $t++) { 
  my @group_ordered_A = ($t == 0) ? (@group_ordered1_A) : (@group_ordered2_A);
  my $caption         = ($t == 0) ? $caption1 : $caption2;

  print("\\begin{table}[t]\n");
  print("\\caption{$caption}\n");
  my $prv_type = "";
  if($t == 0) { 
    print("\\begin{tabular}{llllrlllrrr}\n");
    printf("%-4s & %-3s & %-7s & %-10s & %7s & %-25s & %-6s & %7s & %9s & %9s & %9s \\\\ \n", 
           "",     "",   "sub-", "model", "model", "", "", "", "", "", "");
    printf("%-4s & %-3s & %-7s & %-10s & %7s & %-25s & %-6s & %7s & %9s & %9s & %9s \\\\ \\hline\n", 
           "type", "seg", "type", "accession", "length", "CDS product", "gene", "intron?", "\\#proteins", "\\#coords", "\\#extra");
  }
  else {  # no sub-type column, no #extra column
    print("\\begin{tabular}{lllrlllrr}\n");
    printf("%-4s & %-3s & %-10s & %7s & %-25s & %-6s & %7s & %9s & %9s \\\\ \n", 
           "",     "",    "model", "model", "", "", "", "", "");
    printf("%-4s & %-3s & %-10s & %7s & %-25s & %-6s & %7s & %9s & %9s \\\\ \\hline\n", 
           "type", "seg", "accession", "length", "CDS product", "gene", "intron?", "\\#proteins", "\\#coords");
  }
  my $prv_group = "";
  foreach $group (@group_ordered_A) { 
    my $group2parse = $group;
    $group2parse =~ s/\:\d+$//;
    my ($type, $segment, $subtype) = (undef, undef, undef);
    if($group2parse =~ /^flu(\S)\-seg(\d)(.*)$/) { 
      ($type, $segment, $subtype) = ($1, $2, $3);
      $subtype =~ s/^\.//;
    }
    else { 
      die "ERROR unable to parse group $group";
    }
    $mdl_name   = $model_H{$group};
    my $mdl_name2print = $mdl_name;
    if(defined $version2_H{$mdl_name}) { 
      $mdl_name2print .= ".2"; 
    }
    else { 
      $mdl_name2print .= ".1"; 
    }
    $mdl_name2print =~ s/\_/\\\_/g;
    my $mdl_is_new     = (defined $new_models_H{$mdl_name})   ? 1 : 0;
    my $mdl_is_swapped = ((defined $fluC_new2old_H{$mdl_name}) && ($fluC_new2old_H{$mdl_name} ne $mdl_name)) ? 1 : 0;

    if($mdl_is_new && $mdl_is_swapped) { 
      die "ERROR mdl $mdl_name is new and swapped";
    }
    if($mdl_is_new) { 
      $mdl_name2print = "\\textbf{" . $mdl_name2print . "}"
    }
    if($mdl_is_swapped) { 
      $mdl_name2print = "\\textit{" . $mdl_name2print . "}"
    }

    my $mdl_len = $length_H{$group};
    foreach my $product (sort keys %{$ncoords_HH{$mdl_name}}) { 
      if(($prv_type ne "") && ($prv_type ne $type)) { 
        print "\n";
      }
      my $product2print = ($mnf_HH{$mdl_name}{$product} == 1) ? "\\textit{$product}" : "$product";
      my $gene2print    = ($mnf_HH{$mdl_name}{$product} == 1) ? "\\textit{$product2gene_H{$product}}" : "$product2gene_H{$product}";
      if(($mdl_name eq "NC_036621") && ($product2print eq "nonstructural protein 2")) { 
        $gene2print = "NS2";
      }
      # special hard-coded changes for CY002284, it has only 9 coords in FLAN, but 10 in flan-plus (one of the 3 added proteins is unique coords)
      my $ncoords2print = (($mdl_name eq "CY002284") && ($ncoords_HH{$mdl_name}{$product} == 10)) ? "9" :  $ncoords_HH{$mdl_name}{$product};
      my $nextra2print = ($nblastseqs_HH{$mdl_name}{$product} == $nblastseqs_plus_HH{$mdl_name}{$product}) ? 
          "-" : ($nblastseqs_plus_HH{$mdl_name}{$product} - $nblastseqs_HH{$mdl_name}{$product});
      my $nblastseqs_plus2print = $nblastseqs_plus_HH{$mdl_name}{$product};
      if($nblastseqs_HH{$mdl_name}{$product} == -1) { # special case for flu D seqs
        $nextra2print = "-";
      }
      if($mdl_name eq "NC_006307") { 
        # we can't count the number of proteins for the old FLAN model libraries because coords don't match so we just hardcode this
        $nextra2print = "-";
      }
      if($t == 0) { 
        printf("%-4s & %-3s & %-7s & %-10s & %7s & %-25s & %-6s & %7s & %9d & %9d & %9s \\\\\n", 
               ($group eq $prv_group) ? "" : $type, 
               ($group eq $prv_group) ? "" : $segment, 
               ($group eq $prv_group) ? "" : ($subtype eq "") ? "-" : $subtype,  
               ($group eq $prv_group) ? "" : $mdl_name2print, 
               ($group eq $prv_group) ? "" : $mdl_len, 
               $product2print, $gene2print,
               ($nexons_HH{$mdl_name}{$product} == 1) ? "no" : "yes", 
               $nblastseqs_plus2print,
               $ncoords2print, 
               $nextra2print);
      }
      else { 
        printf("%-4s & %-3s & %-10s & %7s & %-25s & %-6s & %7s & %9d & %9d \\\\\n", 
               ($group eq $prv_group) ? "" : $type, 
               ($group eq $prv_group) ? "" : $segment, 
               ($group eq $prv_group) ? "" : $mdl_name2print, 
               ($group eq $prv_group) ? "" : $mdl_len, 
               $product2print, $gene2print,
               ($nexons_HH{$mdl_name}{$product} == 1) ? "no" : "yes", 
               $nblastseqs_plus2print,
               $ncoords2print);
      }
      $prv_group = $group;
      $prv_type  = $type;
    }
  }
  print("\\end{tabular}\n");
  my $name = ($t == 0) ? "models1" : "models2";
  print("\\label{tbl:$name}\n");
  print("\\end{table}\n");
}
#######################################3




#################################################################
# Subroutine: count_matching_coords_in_fasta
# Incept:     EPN, Tue Jun 27 16:08:14 2023
#
# Purpose:    Given a fasta file path and a 'coords' string, 
#             return the number of sequence names in the fasta file
#             that match the coords string (have the coords string as 
#             a substring of the name)
#
# Arguments:
#  $fasta_file:      fasta sequence file
#  $coords_string:   coords string to look for in sequence names
#
# Returns:  number of sequences in $fasta_file with $coords_string
#           as a subset of their name
#
#################################################################
sub count_matching_coords_in_fasta { 
  my $sub_name = "count_matching_coords_in_fasta";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($fasta_file, $coords_string) = (@_);

  #printf("in $sub_name, fasta_file: $fasta_file, coords_string: $coords_string\n");

  open(FA, $fasta_file) || die "ERROR unable to open $fasta_file for reading";

  my $nseq_read = 0;
  my $ret_val = 0;

  while(my $line = <FA>) { 
    chomp $line;
    if($line =~ /^\>(\S+)/) { 
      my $seqname = $1;
      $nseq_read++;
      if($seqname =~ m/\Q$coords_string/) { 
        $ret_val++;
      }
    }
  }
  close(FA);
  
  #print("HEYA in $sub_name, read $nseq_read seqs, returning $ret_val\n");
  return $ret_val;
}
