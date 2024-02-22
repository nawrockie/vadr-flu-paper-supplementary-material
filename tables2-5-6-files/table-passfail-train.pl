#!/usr/bin/env perl
# 
# table-passfail.pl: create the pass/fail table for the training datasets given the table-passfail.data
#                        
# EPN, Tue Nov 21 13:38:18 2023
# 
#
use strict;
use warnings;
use Getopt::Long;

my $in_data  = ""; # name of input data file

my $usage = "perl table-passfail-train.pl <table-passfail-train.data>\n";

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
#fluA.fr10000   FPVP   fp  9784
  chomp $line;
  my @el_A = split(/\s+/, $line);
  if(scalar(@el_A) != 4) { 
    die "ERROR unable to parse line $line";
  }
  my ($type_nseq, $cat, $minfo, $ct) = (@el_A);
  my @el2_A = split(/\./, $type_nseq);
  my ($type, $nseq, $nseq1, $nseq2, $train_or_test);
  if(scalar(@el2_A) == 2) { 
    ($type, $nseq) = (@el2_A);
  }
  elsif(scalar(@el2_A) == 3) { 
    ($type, $nseq1, $nseq2) = (@el2_A);
    $nseq = $nseq1 . "." . $nseq2;
  }
  elsif(scalar(@el2_A) == 4) { 
    ($train_or_test, $type, $nseq1, $nseq2) = (@el2_A);
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
    die "ERROR unable to parse line $line";
  }
  if($minfo eq "fin") { 
#    $key .= ":::" . "FLAN-plus";
    $key .= ":::" . "final";
  }
  elsif($minfo eq "fo") { 
    $key .= ":::" . "FLAN";
  }
  elsif($minfo eq "fp") { 
    $key .= ":::" . "old-FLAN-plus";
  }
  elsif($minfo eq "nto") { 
    $key .= ":::" . "FLAN-ntonly";
  }
  elsif($minfo eq "nto") { 
    $key .= ":::" . "FLAN-ntonly";
  }
  else { 
    die "ERROR could not parse minfo: $minfo";
  }

  $tbl_HH{$key}{"nseq"} = $nseq2print;
  $tbl_HH{$key}{$cat}   = $ct;
  $tot_H{$key} += $ct;
}  

my @key_order_A = ("flu A GenBank:::final", "flu A GenBank:::FLAN", "flu A GenBank:::FLAN-ntonly", 
                   "flu A EMBL+DDBJ:::final", "flu A EMBL+DDBJ:::FLAN", "flu A EMBL+DDBJ:::FLAN-ntonly", 
                   "flu B GenBank:::final", "flu B GenBank:::FLAN", "flu B GenBank:::FLAN-ntonly", 
                   "flu B EMBL+DDBJ:::final", "flu B EMBL+DDBJ:::FLAN", "flu B EMBL+DDBJ:::FLAN-ntonly", 
                   "flu C GenBank:::final", "flu C GenBank:::FLAN", "flu C GenBank:::FLAN-ntonly", 
                   "flu C EMBL+DDBJ:::final", "flu C EMBL+DDBJ:::FLAN", "flu C EMBL+DDBJ:::FLAN-ntonly"); 

my $caption = "\\textbf{Comparison of pass/fail outcomes for FLAN and VADR on the influenza sequence training datasets.}";

print("\\begin{table}[t]\n");
print("\\caption{$caption}\n");
print("\\begin{tabular}{lrlrrrrr}\n");
printf("%-30s & %10s & %-15s & %9s & %9s & %9s & %9s & %9s \\\\\n",
       "",        "",         "VADR",  "fraction",           "",     "",         "\\#FLAN-pass", "\\#FLAN-fail"); 
printf("%-30s & %10s & %-15s & %9s & %9s & %9s & %9s & %9s \\\\ \\hline\n",
       "dataset", "num seqs", "model", "pass both", "\\#pass both", "\\#fail both", "VADR-fail", "VADR-pass"); 
my $tot_nseq = 0;
my $tot_FPVP = 0;
my $tot_FFVF = 0;
my $tot_FPVF = 0;
my $tot_FFVP = 0;
foreach my $key (@key_order_A) { 
  my ($dataset, $model) = split(":::", $key);
  my $dataset2print = ($model eq "final") ? $dataset : "";
  my $nseq2print    = ($model eq "final") ? $tbl_HH{$key}{"nseq"} : "";
  printf("%-30s & %10s & %-15s & %9.3f & %9d & %9d & %9d & %9d \\\\ ",
         $dataset2print, $nseq2print, $model, 
         ($tbl_HH{$key}{"FPVP"} / $tot_H{$key}),
         $tbl_HH{$key}{"FPVP"}, 
         $tbl_HH{$key}{"FFVF"}, 
         $tbl_HH{$key}{"FPVF"}, 
         $tbl_HH{$key}{"FFVP"});
#         sprintf("%9d (%.3f)", $tbl_HH{$key}{"FPVP"}, $tbl_HH{$key}{"FPVP"} / $tot_H{$key}), 
#         sprintf("%9d (%.3f)", $tbl_HH{$key}{"FFVF"}, $tbl_HH{$key}{"FFVF"} / $tot_H{$key}), 
#         sprintf("%9d (%.3f)", $tbl_HH{$key}{"FPVF"}, $tbl_HH{$key}{"FPVF"} / $tot_H{$key}), 
#         sprintf("%9d (%.3f)", $tbl_HH{$key}{"FFVP"}, $tbl_HH{$key}{"FFVP"} / $tot_H{$key})); 
  $tot_nseq_H{$model} += $tbl_HH{$key}{"nseq"};
  $tot_FPVP_H{$model} += $tbl_HH{$key}{"FPVP"};
  $tot_FFVF_H{$model} += $tbl_HH{$key}{"FFVF"};
  $tot_FPVF_H{$model} += $tbl_HH{$key}{"FPVF"};
  $tot_FFVP_H{$model} += $tbl_HH{$key}{"FFVP"};
  printf("\n");
  if($model eq "FLAN-ntonly") {
    printf(" & & & & & & \\\\\n");
  }
}
printf("\\hline\n");

my @model_A = ("final", "FLAN", "FLAN-ntonly");
foreach my $model (@model_A) { 
  my $dataset2print = ($model eq "final") ? "total" : "";
  my $nseq2print    = ($model eq "final") ? $tot_nseq_H{$model} : "";
  printf("%-30s & %10s & %-15s & %9.3f & %9d & %9d & %9d & %9d \\\\\n",
         $dataset2print, $nseq2print, $model, 
         ($tot_FPVP_H{$model} / $tot_nseq_H{$model}), 
         $tot_FPVP_H{$model}, 
         $tot_FFVF_H{$model},
         $tot_FPVF_H{$model},
         $tot_FFVP_H{$model});
#         sprintf("%9d (%.3f)", $tot_FPVP_H{$model}, $tot_FPVP_H{$model} / $tot_nseq_H{$model}), 
#         sprintf("%9d (%.3f)", $tot_FFVF_H{$model}, $tot_FFVF_H{$model} / $tot_nseq_H{$model}), 
#         sprintf("%9d (%.3f)", $tot_FPVF_H{$model}, $tot_FPVF_H{$model} / $tot_nseq_H{$model}), 
#         sprintf("%9d (%.3f)", $tot_FFVP_H{$model}, $tot_FFVP_H{$model} / $tot_nseq_H{$model}));
}
print("\\end{tabular}\n");
print("\\label{tbl:passfail-train}\n");
print("\\end{table}\n");
