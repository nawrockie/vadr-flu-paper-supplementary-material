my $identical = 0;
my $different = 0;
my $tot = 0;
while($line = <>) { 
  chomp $line;
  if($line =~ /^\S+\s+\S+\s+\d+$/) { 
    my @el_A = split(/\s+/, $line);
    if($el_A[1] eq "identical") { 
      $identical += $el_A[2];
    }
    elsif($el_A[1] eq "different") { 
      $different += $el_A[2];
    }
    $tot += $el_A[2];
  }
  elsif($line =~ m/\w/) { 
    die "ERROR unable to parse line 1: $line";
  }
}
printf("identical: %10d (%.3f)\n", $identical, ($identical / $tot));
printf("different: %10d (%.3f)\n", $different, ($different / $tot));
printf("total:     %10d (%.3f)\n", $tot,       ($tot / $tot));

