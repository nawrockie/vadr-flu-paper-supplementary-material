#!/usr/bin/env perl
# 
# table-passfail-test.pl: create the pass/fail table for the testing datasets given the table-passfail-test.data
#                        
# EPN, Tue Nov 21 13:38:18 2023
# 
#
use strict;
use warnings;
use Getopt::Long;

my $in_data  = ""; # name of input data file

my $usage = "perl table-passfail-test.pl <table-passfail-test.data>\n";

if(scalar(@ARGV) != 1) { die $usage; }
my ($data_file) = (@ARGV);

my %tbl_HH = ();
my %tot_H = ();
my %tot_nseq_H = ();
my %tot_FPVP_H = ();
my %tot_FPVF_H = ();
my %tot_FFVP_H = ();
my %tot_FFVF_H = ();
my %emb_nseq_H = (); # key is number of total seqs, value is number that are EMBL/ENA
my %dbj_nseq_H = (); # key is number of total seqs, value is number that are DDBJ

open(IN, $data_file) || die "ERROR unable to open $data_file for reading";
while(my $line = <IN>) { 
  chomp $line;
  my @el_A = split(/\s+/, $line);
  if(scalar(@el_A) == 3) { 
    # 3 tokens, these tell us how many dbj vs ena seqs there are
    #ts.fluA.fr10000 emb 3965
    my ($type_nseq, $db, $nseq_db) = (@el_A);
    my $nseq;
    if($type_nseq =~ /(\d+)$/) {
      $nseq = $1;
    }
    else { 
      die "ERROR, could not parse $type_nseq token line $line"; 
    }
    if($db eq "emb") { 
      if(defined $emb_nseq_H{$nseq}) { 
        die "ERROR, read two values for ENA nseq $nseq"; 
      }
      $emb_nseq_H{$nseq} = $nseq_db;
    }
    elsif($db eq "dbj") { 
      if(defined $dbj_nseq_H{$nseq}) { 
        die "ERROR, read two values for DBJ nseq $nseq"; 
      }
      $dbj_nseq_H{$nseq} = $nseq_db;
    }
  }
  elsif(scalar(@el_A) == 4) { 
    #ts.fluA.fr10000   FPVP   fp  9790
    chomp $line;
    my @el_A = split(/\s+/, $line);
    my ($type_nseq, $cat, $minfo, $ct) = (@el_A);
    if($minfo eq "fin") { # we are only concerned with final
      my @el2_A = split(/\./, $type_nseq);
      my ($date, $type, $nseq, $nseq1, $nseq2);
      if(scalar(@el2_A) == 3) { 
        ($date, $type, $nseq) = (@el2_A);
      }
      elsif(scalar(@el2_A) == 4) { 
        ($date, $type, $nseq1, $nseq2) = (@el2_A);
        $nseq = $nseq1 . "." . $nseq2;
      }
      else { 
        die "ERROR unable to parse type_nseq: $type_nseq";
      }

      my $key = "";
      my $nseq2print = "";
      if($type eq "fluA") { 
        $key .= "flu A";
      }
      elsif($type eq "fluB") { 
        $key .= "flu B";
      }
      elsif($type eq "fluC") { 
        $key .= "flu C";
      }
      if($nseq =~ /^gb.fr(\d+)/) { 
        $key .= " GenBank";
        $nseq2print .= " " . $1;
      }
      elsif($nseq =~ /^gb.(\d+)/) { 
        $key .= " GenBank";
        $nseq2print .= " " . $1;
      }
      elsif($nseq =~ /^fr(\d+)/) { 
        $key .= " ENA+DDBJ";
        $nseq2print .= " " . $1;
      }
      elsif($nseq =~ /^nongb.(\d+)/) { 
        $key .= " ENA+DDBJ";
        $nseq2print .= " " . $1;
      }
      else { 
        die "ERROR unable to parse line (2) $line";
      }
      if($date eq "ts") { 
        $key .= ":::" . "existing";
      }
      elsif($date eq "test1") { 
        $key .= ":::" . "existing";
      }
      elsif($date eq "ts23") { 
        $key .= ":::" . "new";
      }
      elsif($date eq "test2") { 
        $key .= ":::" . "new";
      }
      else { 
        die "ERROR could not parse date: $date";
      }

      $tbl_HH{$key}{"nseq"} = $nseq2print;
      $tbl_HH{$key}{$cat}   = $ct;
      $tot_H{$key} += $ct;
      #printf("tot_H{$key} is $tot_H{$key}\n");
    }  
  }
  elsif(scalar(@el_A) != 4) { 
    die "ERROR could not parse line: $line";
  }
}
my @key_order_A = ("flu A GenBank:::existing","flu A GenBank:::new",
                   "flu A ENA+DDBJ:::existing","flu A ENA+DDBJ:::new",
                   "flu B GenBank:::existing","flu B GenBank:::new",
                   "flu B ENA+DDBJ:::existing","flu B ENA+DDBJ:::new",
                   "flu C GenBank:::existing","flu C GenBank:::new",
                   "flu C ENA+DDBJ:::existing","flu C ENA+DDBJ:::new");
                   
my $caption = "\\textbf{Comparison of pass/fail outcomes for FLAN and VADR on the influenza sequence testing datasets.} For ENA+DDBJ datasets, number of ENA and DDBJ sequences in the set are indicated in parantheses in \"num seqs\" column, ENA listed first and DDBJ listed second.";

print("\\begin{table}[t]\n");
print("\\caption{$caption}\n");
print("\\begin{tabular}{lrlrrrrr}\n");
printf("%-30s & %-15s & %-10s & %9s & %9s & %9s & %9s & %9s \\\\\n",
       "",        "",      "", "fraction",  "",           "",         "\\#FLAN-pass", "\\#FLAN-fail"); 
printf("%-30s & %-15s & %-10s & %9s & %9s & %9s & %9s & %9s \\\\ \\hline\n",
       "dataset", "release date", "num seqs", "pass both", "\\#pass both", "\\#fail both", "VADR-fail", "VADR-pass"); 
my $tot_nseq = 0;
my $tot_FPVP = 0;
my $tot_FFVF = 0;
my $tot_FPVF = 0;
my $tot_FFVP = 0;
foreach my $key (@key_order_A) { 
#  printf("key: $key\n");
  my ($dataset, $date) = split(":::", $key);
  my $dataset2print = ($date eq "existing") ? $dataset : "";
  my $date2print = ($date eq "existing") ? "before 3/18/2023" : "after 3/17/2023"; 
  my $nseq2print = $tbl_HH{$key}{"nseq"};
  my $nseq2use = $nseq2print;
  $nseq2use =~ s/^\s+//;
  my $emb_seq2print = "";
  my $dbj_seq2print = "";
  if(($key !~ m/GenBank/) && ($nseq2use ne "")) { 
    if(! defined $emb_nseq_H{$nseq2use}) { 
      die "ERROR emb_nseq_H{$nseq2use} undefined"; 
    }
    if(! defined $dbj_nseq_H{$nseq2use}) { 
      die "ERROR dbj_nseq_H{$nseq2use} undefined"; 
    }
    $emb_seq2print = $emb_nseq_H{$nseq2use};
    $dbj_seq2print = $dbj_nseq_H{$nseq2use};
  }
  printf("%-30s & %-15s & %10d%12s & %9.3f & %9d & %9d & %9d & %9d \\\\ ",
         $dataset2print, $date2print, $nseq2print, 
         ($emb_seq2print eq "") ? "" : sprintf(" (%d+%d)", $emb_seq2print, $dbj_seq2print), 
         ($tbl_HH{$key}{"FPVP"} / $tot_H{$key}), 
         $tbl_HH{$key}{"FPVP"},
         $tbl_HH{$key}{"FFVF"},
         $tbl_HH{$key}{"FPVF"},
         $tbl_HH{$key}{"FFVP"});
  #sprintf("%9d (%.3f)", $tbl_HH{$key}{"FPVP"}, $tbl_HH{$key}{"FPVP"} / $tot_H{$key}), 
  #sprintf("%9d (%.3f)", $tbl_HH{$key}{"FFVF"}, $tbl_HH{$key}{"FFVF"} / $tot_H{$key}), 
  #sprintf("%9d (%.3f)", $tbl_HH{$key}{"FPVF"}, $tbl_HH{$key}{"FPVF"} / $tot_H{$key}), 
  #sprintf("%9d (%.3f)", $tbl_HH{$key}{"FFVP"}, $tbl_HH{$key}{"FFVP"} / $tot_H{$key})); 
  $tot_nseq_H{$date} += $tbl_HH{$key}{"nseq"};
  $tot_FPVP_H{$date} += $tbl_HH{$key}{"FPVP"};
  $tot_FFVF_H{$date} += $tbl_HH{$key}{"FFVF"};
  $tot_FPVF_H{$date} += $tbl_HH{$key}{"FPVF"};
  $tot_FFVP_H{$date} += $tbl_HH{$key}{"FFVP"};
  printf("\n");
  if($key =~ m/new/) { 
    printf(" & & & & & & \\\\\n");
  }
}
printf("\\hline\n");

my @date_A = ("existing", "new");
foreach my $date (@date_A) { 
  my $date2print = ($date eq "existing") ? "before 3/18/2023" : "after 3/17/2023"; 
  my $dataset2print = ($date eq "existing") ? "total" : "";
  printf("%-30s & %-15s & %22d & %9.3f & %9d & %9d & %9d & %9d \\\\\n",
         $dataset2print, $date2print, $tot_nseq_H{$date},
         ($tot_FPVP_H{$date} / $tot_nseq_H{$date}), 
         $tot_FPVP_H{$date},
         $tot_FFVF_H{$date},
         $tot_FPVF_H{$date},
         $tot_FFVP_H{$date});
         #sprintf("%9d (%.3f)", $tot_FPVP_H{$date}, $tot_FPVP_H{$date} / $tot_nseq_H{$date}), 
         #sprintf("%9d (%.3f)", $tot_FFVF_H{$date}, $tot_FFVF_H{$date} / $tot_nseq_H{$date}), 
         #sprintf("%9d (%.3f)", $tot_FPVF_H{$date}, $tot_FPVF_H{$date} / $tot_nseq_H{$date}), 
         #sprintf("%9d (%.3f)", $tot_FFVP_H{$date}, $tot_FFVP_H{$date} / $tot_nseq_H{$date}));
}
print("\\end{tabular}\n");
print("\\label{tbl:passfail-test}\n");
print("\\end{table}\n");
