package CodonFunc;

# Semaphore controlling single thread access to 
# SpliceAI Tensor flow resources
my $sps_lock :shared = 0;

my $debug = 0;

my %bicodon_relative_adap;

my $myCodonTable   = Bio::Tools::CodonTable->new();
$myCodonTable->id(1);  


###############################
# Reads bicodon usage

sub load_bicodon_usage{
    my $bicodon_usage_file = shift;

    my %bicodon_usage;
    my @bicodon_pairs;
    my @bicodon_values;

    my $line_count = 1;
    open IN, "$bicodon_usage_file" or die "Unable to open bicodon usage $bicodon_usage_file!\n";
    while(<IN>){
    my $line = uc($_);
    chomp($line);
    if( $line_count == 1 ){
        @bicodon_pairs = split "\t", $line;
    }elsif( $line_count == 2 ){
        @bicodon_values = split "\t", $line;
    }else{
        die "Expected only two lines on file $bicodon_usage_file";
    }
    $line_count++;
    }
    close(IN);

    # Convert bicodon to pair of aminoacids
    # And finds the bicodon with the highest count/frequency per each pair of AA
    my @pairs_aa;
    my %pair_aa_max_value;
    for my $index ( 0 .. $#bicodon_pairs ){

        my $bicodon_pair = $bicodon_pairs[ $index ];

        my $pair_aa = $myCodonTable->translate( $bicodon_pairs[ $index ] );
        push( @pairs_aa, $pair_aa );

        if( not defined $pair_aa_max_value{$pair_aa} ){
        	$pair_aa_max_value{$pair_aa} = $bicodon_values[ $index ];
        }else{
        	$pair_aa_max_value{$pair_aa} = $bicodon_values[ $index ] if( $pair_aa_max_value{$pair_aa} < $bicodon_values[ $index ] );
        }

        #print "Bicodon: " . $bicodon_pair . " Pair AA:" . $pair_aa . "\n";
        #getc();
    }

    # Calculate relative adaptiveness for each bicodon
    # https://en.wikipedia.org/wiki/Codon_Adaptation_Index

    for my $index ( 0 .. $#bicodon_pairs ){
        $bicodon_relative_adap{ $bicodon_pairs[ $index ] } = $bicodon_values[ $index ] / $pair_aa_max_value{ $pairs_aa[ $index ] };
        #print $bicodon_pairs[ $index ] . " " . $pairs_aa[ $index ] . " " .  
        #  $bicodon_values[ $index ] . " " . $pair_aa_max_value{ $pairs_aa[ $index ] } . " " 
        #  . $bicodon_relative_adap{ $bicodon_pairs[ $index ] } . "\n";
    }

}


#########################
# Logging 


sub reset_log{
  my $prefix_log = shift;
  `rm -f $prefix_log.log`
}

sub log_msg{
  my $prefix_log = shift;
  my $msg = shift;
  open OUT, ">$prefix_log.log" or die $!;
  print OUT  $msg . "\n";
  close(OUT);
}


#########################
# Functions to optmize

sub get_bai{
  my ($gene_seq) = @_;

  my $sum_logs = 0;
  my $count_bicodons = 0;
  for( my $bicodon_start = 0; $bicodon_start <= length( $gene_seq ) - 7; $bicodon_start += 3 ){
    my $curr_bicodon = uc(substr( $gene_seq, $bicodon_start, 6 ));


    # Using sum of logs to calculate geometric mean. Otherwise we will have underflow when
    # calculating BAI for long nt sequences
    $sum_logs += log( $bicodon_relative_adap{ $curr_bicodon } );

    if(0){
    print "Pos: $bicodon_start bicodon:  $curr_bicodon " . 
          "Adap.index: " . $bicodon_relative_adap{ $curr_bicodon } . 
          " Sum logs: $sum_logs " . 
          "\n";
    getc();
    }

    $count_bicodons++;
  }

  my $bai = exp( $sum_logs / $count_bicodons ) ;

  #print "BAI: $bai sum_logs: $sum_logs count: $count_bicodons\n";
  #getc();

  return $bai;
}

sub get_cai{
  my ($gene_seq, $codon_usage_file, $tempdir) = @_;

  my ($fh, $filename) = File::Temp::tempfile(DIR => $tempdir );
  print "CAI: " . $filename . "\n" if $debug;

  print $fh ">nt_seq\n";
  print $fh $gene_seq;
  close($fh);

  my $out = `cai $filename -cfile $codon_usage_file -outfile $filename.out.cai >& /dev/null`;
  my $cai =  `awk '{print \$NF}' $filename.out.cai`;

  `rm $filename $filename.out.cai` if $remove_temp;

  chomp($cai);

  return $cai;
}

sub get_num_cpg{
  my ($gene_seq, $tempdir) = @_;

  #print "Gene>>$gene_seq<<\n";
  #print "Temp>>$tempdir<<\n";
  #getc();
  my $temp = lc( $gene_seq );
  $gene_seq =  $temp;

  my @c = $gene_seq =~ /cg/g;
  my $count = @c;
  return $count;
}

sub get_pas{
  my ($gene_seq, $tempdir) = @_;

}

sub get_cpg{
  my ($gene_seq, $tempdir) = @_;

  #print "Gene>>$gene_seq<<\n";
  #print "Temp>>$tempdir<<\n";
  #getc();
  my $temp = lc( $gene_seq );
  $gene_seq =  $temp;

  my @c = $gene_seq =~ /cg/g;
  my $count = @c * 2;
  my $seq_len = length( $gene_seq );

  my $cpg = ($seq_len - $count) / $seq_len;

  return $cpg;
}

sub get_pas{
  my ($gene_seq, $tempdir) = @_;

  my ($fh, $filename) = File::Temp::tempfile(DIR => $tempdir );
  log_msg($filename, "Starting get_pas...");
  print "PAS: " . $filename . "\n" if $debug;
  print $fh ">nt_seq\n";
  print $fh $gene_seq;
  close($fh);

  log_msg($filename, "Starting DeepPasta...");
  my $out = `./deep_pasta.sh $filename $filename 2>&1`;
  log_msg($filename, $out );
  if( $debug ){
    print STDERR $out;
  }
  log_msg($filename, "Done DeepPasta");

  my $pas =  `awk '{print \$NF}' $filename.pas.txt`;

  # Removing all generated files
  `rm $filename $filename.pas.txt $filename.fai` if $remove_temp;
  `rm $filename.bed $filename.stag.bed $filename.stag.fasta` if $remove_temp;
  `rm $filename.ss.txt $filename.ss.comb.txt $filename.ss.filt.txt` if $remove_temp;
  `rm $filename.ss.per_nt.txt $filename.report.txt` if $remove_temp;

  chomp($pas);

  my $norm_pas = ($pas - 0.95 > 0 ) ? ($pas - 0.95) * 20 : 0;
  reset_log($filename);

  return 1 - $norm_pas;
}

#sub get_sps{
#  my ($gene_seq, $tempdir) = @_;

#  # direct
#  my $direct_sps = get_sps_core_function( $gene_seq, $tempdir );

  # revcomp
#  my $gene_seq_revcomp = reverse( $gene_seq );
#  $gene_seq_revcomp =~ tr/ATGCatgc/TACGtacg/; 
#  my $revcomp_sps = get_sps_core_function( $gene_seq_revcomp, $tempdir );

#  return $direct_sps * $revcomp_sps;
#}

sub get_sps{
  my ($gene_seq, $tempdir) = @_;

  my $out;
  # Only one thread will pass through the next commands
 
  my ($fh, $filename) = File::Temp::tempfile(DIR => $tempdir );
  print "SPS: " . $filename . "\n" if $debug;
  print $fh ">nt_seq\n";
  print $fh $gene_seq;
  close($fh);

  # Only one thread will pass through the next commands
  {
    lock($sps_lock);
    if( $debug ){
      $out = `./spliceAI.py -i $filename -o $filename.splice.txt`;
    }else{
      $out = `./spliceAI.py -i $filename -o $filename.splice.txt >& /dev/null`;
    }
  }

  $out =  `tail -n +2 $filename.splice.txt | awk '\$2 > 0.8 {print "acceptor",\$1,\$2} \$3 > 0.8 {print "donnor",\$1,\$3} ' OFS="\t" > $filename.splice.filtered.txt`;
  
  # Read splice probabilities
  my %spl_prob;
  my @donnor_coords;
  my @acceptor_coords;
  my @probs;


  $out = `wc -l  $filename.splice.filtered.txt | awk '{print \$1}'`;
  
  chomp($out);

  if( $out == 0 ){
    return 1;
  }

  open( $fh, '<', "$filename.splice.filtered.txt" )
    or die "Could not open file $filename.splice.filtered.txt $!";
  
  while (my $line = <$fh>) {
    chomp $line;
    my ($type, $coord, $prob) = split(/\t/, $line);

    $spl_prob{$type}{$coord} = $prob; 
    push( @donnor_coords, $coord ) if $type eq 'donnor';
    push( @acceptor_coords, $coord ) if $type eq 'acceptor';
    push( @probs, $prob );
  }  
  close($fh);
  `rm $filename $filename.splice.txt $filename.splice.filtered.txt` if $remove_temp;

  my $donnor_acceptor_found = 0;

  # Search for donnor followed by acceptor
  # Score will be set to 0 if found
  foreach my $curr_donor ( @donnor_coords ){
    foreach my $curr_acceptor ( @acceptor_coords ){
      if( $curr_donor < $curr_acceptor ){
        $donnor_acceptor_found = 1;
        last;
      }
      last if $donnor_acceptor_found;
    }
  }

  # Probability of a splice signal
  my $avg_prob;
  if( $donnor_acceptor_found == 1 ){
    $avg_prob = 0.99;
  }else{
    if( scalar( @probs ) == 0 ){
      $avg_prob = 0;
    }else{
      $avg_prob = List::Util::sum(@probs) / scalar( @probs );
    }
  }

  return (1 - $avg_prob);
}


sub get_res{
  my ($gene_seq, $re_comma_sep, $tempdir) = @_;

  my ($fh, $filename) = File::Temp::tempfile(DIR => $tempdir );
  log_msg($filename, "Starting get_res...");
  print "RES: " . $filename . "\n" if $debug;
  print $fh ">nt_seq\n";
  print $fh $gene_seq;
  close($fh);


  log_msg($filename, "Starting Restrict...");
  my $out = `restrict -sequence $filename -enzymes $re_comma_sep -sitelen 4 -fragment -outfile $filename.re.txt 2>&1`;
  log_msg($filename, $out );
  if( $debug ){
    print STDERR $out;
  }
  log_msg($filename, "Done Restrict");

  $out = `grep -v "#" $filename.re.txt | awk '\$0 != "" {print \$0}' | tail -n +2 | wc -l`;
  `rm $filename $filename.re.txt` if $remove_temp;

  chomp($out);

  if( $out == 0 ){
    reset_log($filename);
    return 1;
  }else{
    reset_log($filename);
    return 0.01;
  }

}

return 1;
