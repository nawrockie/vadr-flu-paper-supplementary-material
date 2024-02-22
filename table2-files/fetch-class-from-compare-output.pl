my $usage = "perl fetch-class-from-compare-output.pl <parse-and-compare-flan-and-vadr-output.pl output file> <pass/fail class (e.g. FPVF)>";
if(scalar(@ARGV) != 2) { 
  die $usage; 
}
my ($compare_file, $class) = (@ARGV);

if(($class ne "FPVP") && ($class ne "FPVF") && ($class ne "FFVP") && ($class ne "FFVF")) { 
  die "ERROR class must be one of FPVP FPVF FFVP FFVF";
}

open(COMPARE, $compare_file) || die "ERROR unable to open $compare_file for reading";
my @out_A = ();
while($line = <COMPARE>) { 
  chomp $line;
  if($line =~ /SEQUENCE/) { 
    if($line =~ m/$class/) { 
      foreach my $out_line (@out_A) { 
        print $out_line . "\n";
      }
      print $line . "\n";
    }
    @out_A = ();
  }
  else {
    push(@out_A, $line);
  }
}
close(COMPARE);
