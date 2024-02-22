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

open(IN, $data_file) || die "ERROR unable to open $data_file for reading";
while(my $line = <IN>) { 
#ts.fluA.fr10000   FPVP   fp  9790
  chomp $line;
  my @el_A = split(/\s+/, $line);
  if(scalar(@el_A) != 4) { 
    die "ERROR unable to parse line $line";
  }
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
      $key .= " EMBL+DDBJ";
      $nseq2print .= " " . $1;
    }
    elsif($nseq =~ /^nongb.(\d+)/) { 
      $key .= " EMBL+DDBJ";
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
my @key_order_A = ("flu A GenBank:::existing","flu A GenBank:::new",
                   "flu A EMBL+DDBJ:::existing","flu A EMBL+DDBJ:::new",
                   "flu B GenBank:::existing","flu B GenBank:::new",
                   "flu B EMBL+DDBJ:::existing","flu B EMBL+DDBJ:::new",
                   "flu C GenBank:::existing","flu C GenBank:::new",
                   "flu C EMBL+DDBJ:::existing","flu C EMBL+DDBJ:::new");
                   
my $caption = "\\textbf{Comparison of pass/fail outcomes for FLAN and VADR on the influenza sequence testing datasets.}";

print("\\begin{table}[t]\n");
print("\\caption{$caption}\n");
print("\\begin{tabular}{lrrrrrrr}\n");
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
  printf("key: $key\n");
  my ($dataset, $date) = split(":::", $key);
  my $dataset2print = ($date eq "existing") ? $dataset : "";
  my $date2print = ($date eq "existing") ? "before 3/18/2023" : "after 3/17/2023"; 
  printf("%-30s & %-15s & %10d & %9.3f & %9d & %9d & %9d & %9d \\\\ ",
         $dataset2print, $date2print, $tbl_HH{$key}{"nseq"},
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
  printf("%-30s & %-15s & %10d & %9.3f & %9d & %9d & %9d & %9d \\\\\n",
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
