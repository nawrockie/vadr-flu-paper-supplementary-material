#!/usr/bin/env perl
# 
# fig-nseq-data2R.pl: create an R script that will generate the number of sequences figure (nseq).
#                        
# EPN, Thu Jun 13 14:23:37 2024
# 
#
# colors obtained from color brewer 2.0
# https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
my $mygreen  = "#1b9e77";
my $myorange = "#d95f02";
my $mypurple = "#7570b3";

my %nseq_HH = ();
while($line = <>) { 
# A 2005 18619
  chomp $line;
  if($line !~ m/^\#/) { 
    my @el_A = split(/\s+/, $line);
    if(scalar(@el_A) != 3) { 
      die "ERROR not 3 tokens on line: $line";
    }
    my ($type, $year, $nseq) = (@el_A);
    if(! defined $nseq_HH{$type}) { 
      $nseq_HH{$type} = ();
    }
    $nseq_HH{$type}{$year} = $nseq;
  }
}

# output R data
my $year; 
my $R_line = "years <- c(";
for($year = 1999; $year <= 2023; $year++) { 
  $R_line .= sprintf("%d", ($year+1));
  if($year < 2023) { 
    $R_line .= ", ";
  }
}
$R_line .= ")\n";
print $R_line;

foreach $type ("A", "B", "C", "D") { 
  $R_line = "$type <- c(";
  for($year = 1999; $year <= 2023; $year++) { 
    $R_line .= "$nseq_HH{$type}{$year}";
    if($year < 2023) { 
      $R_line .= ", ";
    }
  }
  $R_line .= ")\n";
  print $R_line;
}

# create the graph
print("pdf('nseq.pdf', width=8,height=5)\n");
print("plot(years, A, type=\"o\", col=\"black\", xlim=c(2000,2026), ylim=c(0, 1000000), pch=1, xaxt='n', yaxt='n')\n");
print("axis(side=1, at=c(2000, 2005, 2010, 2015, 2020, 2024), labels=c(\"2000\",\"2005\",\"2010\",\"2015\",\"2020\",\"2024\"))\n");
print("axis(side=2, at=c(10000, 100000, 250000, 500000, 750000, 1000000), labels=c(\"10000\",\"100000\",\"250000\",\"500000\",\"750000\",\"1000000\"))\n");
print("axis(side=4, at=c(10000, 100000, 250000, 500000, 750000, 1000000), labels=c(\"10000\",\"100000\",\"250000\",\"500000\",\"750000\",\"1000000\"))\n");
print("lines(years, B, type=\"o\", col=\"$mygreen\", pch=5)\n");
print("lines(years, C, type=\"o\", col=\"$myorange\", pch=0)\n");
print("lines(years, D, type=\"o\", col=\"$mypurple\", pch=3)\n");
print("legend(2000, 900000, legend=c(\"influenza A\", \"influenza B\", \"influenza C\", \"influenza D\"), col=c(\"black\", \"$mygreen\", \"$myorange\", \"$mypurple\"), lty=c(1,1,1,1), pch=c(1,5,0,3), cex=0.8)\n");
print("dev.off()\n");

# to execute in R and create nseq.pdf:
# $ R --vanilla < <output from this script>
