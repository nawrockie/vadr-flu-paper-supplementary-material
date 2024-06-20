#!/usr/bin/env perl
# 
# table-nseqs-data2tex.pl: create the latex table-nseqs table from input data file
#                        
# EPN, Tue Jun 27 12:34:48 2023
# 
#
use strict;
use warnings;
use Getopt::Long;

my $in_data  = ""; # name of input data file
#my $date = "20230626";
#my $date = "20231130";
my $date = "20240220";

my $usage = "perl table-nseqs-data2tex.pl <data file>";
#$usage .= "\tOPTIONS:\n";
#$usage .= "\t\t-T <n>       : minimum bit score to include is <n>\n";
#$usage .= "\t\t-E <x>       : maximum E-value to include is <x>\n";

#&GetOptions( "T=s"       => \$minscore,
#             "E=s"       => \$maxevalue);

if(scalar(@ARGV) != 1) { die $usage; }
($in_data) = @ARGV;

open(IN, $in_data) || die "ERROR unable to open $in_data for reading";
my @month_order_A = ();

my %nseq_HHH = ();

while(my $line = <IN>) { 
#fluC 2022 gb 39
  if($line !~ m/^\#/) { 
    chomp $line;
    my @el_A = split(/\s+/, $line); 
    if(scalar(@el_A) != 4) { 
      die "ERROR unable to parse $line into 4 tokens";
    }
    my ($type, $year, $db, $nseq) = (@el_A);
    if(! defined $nseq_HHH{$db}) { 
      %{$nseq_HHH{$db}} = ();
    }
    if(! defined $nseq_HHH{$db}{$year}) { 
      %{$nseq_HHH{$db}{$year}} = ();
    }
    $nseq_HHH{$db}{$year}{$type} = $nseq;
#    printf("added nseq_HHH{$db}{$year}{$type} as $nseq_HHH{$db}{$year}{$type}\n");
  }
}
close(IN);

my $caption = "\\textbf{Number of influenza virus sequences deposited in GenBank database since 2018.}";
$caption   .= " Sequence counts were obtained using the NCBI Virus Data Hub filtering by release date. ``all'' row includes total counts of all sequences (including before 2018). ``GenBank'' row indicates the fraction of all influenza sequences in INSDC databases that were submitted to GenBank (not to EMBL or DDBJ). NCBI taxonomy ids: 11320 (influenza A); 11520 (influenza B); 11552 (influenza C); 1511084 (influenza D).";

print("\\begin{table}[t]\n");
print("\\caption{$caption}\n");
print("\\begin{tabular}{rrrrr}\n");
#printf("%-18s&%-18s&%-18s&%-18s&%-18s\\\\ \\hline\n", "year", "influenza A", "influenza B", "influenza C", "influenza D");
printf("%-18s&%-18s&%-18s&%-18s&%-18s\\\\ \\hline\n", "year", "influenza A", "B", "C", "D");

my @year_order_A = ("2018", "2019", "2020", "2021", "2022", "2023", "all.$date");
my @type_order_A = ("fluA", "fluB", "fluC", "fluD");

foreach my $year (@year_order_A) { 
  my $year2print = $year;
  $year2print =~ s/\.$date//;
  printf("%18s", $year2print);
  foreach my $type (@type_order_A) { 
    if(! defined $nseq_HHH{"gb"}{$year}{$type}) { 
      die "ERROR nseq_HHH{gb}{$year}{$type} not defined";
    }
    printf("&%18s", sprintf("%6d", $nseq_HHH{"gb"}{$year}{$type}));
  }
  print(" \\\\");
  if($year2print eq "2023") { 
    print(" \\hline");
  }
  print("\n");
}

# print final fraction line
my $year = "all.$date";
my $year2print = $year;
$year2print =~ s/\.$date//;
printf("%18s", "GenBank");
foreach my $type (@type_order_A) { 
  my $fract = $nseq_HHH{"gb"}{$year}{$type} / $nseq_HHH{"all"}{$year}{$type};
  printf("&%18s", sprintf("%5.3f", $fract));
}
print(" \\\\\n");


print("\\end{tabular}\n");
print("\\label{tbl:nseqs}\n");
print("\\end{table}\n");
