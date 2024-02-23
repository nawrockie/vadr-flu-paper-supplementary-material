my $coords_same = 0;
my $coords_diff = 0;
my $vadronly = 0;
my $flanonly = 0;
my $tot = 0;
while($line = <>) { 
  chomp $line;
  if($line =~ /^\S+\s+\S+\s+\d+$/) { 
    my @el_A = split(/\s+/, $line);
    if($el_A[1] eq "coords-same") { 
      $coords_same += $el_A[2];
    }
    elsif($el_A[1] eq "coords-diff") { 
      $coords_diff += $el_A[2];
    }
    elsif($el_A[1] eq "vadronly") { 
      $vadronly += $el_A[2];
    }
    elsif($el_A[1] eq "flanonly") { 
      $flanonly += $el_A[2];
    }
    else { 
      die "ERROR unable to parse line 0: $line";
    }
    $tot += $el_A[2];
  }
  elsif($line =~ m/\w/) { 
    die "ERROR unable to parse line 1: $line";
  }
}
printf("coords-same: %10d (%.3f)\n", $coords_same, ($coords_same / $tot));
printf("coords-diff: %10d (%.3f)\n", $coords_diff, ($coords_diff / $tot));
printf("vadronly:    %10d (%.3f)\n", $vadronly,    ($vadronly    / $tot));
printf("flanonly:    %10d (%.3f)\n", $flanonly,    ($flanonly    / $tot));
printf("total:       %10d (%.3f)\n", $tot,         ($tot / $tot));

