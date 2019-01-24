task countFrags {
  String finalBAM
  String Root
  String Outdir
  String name
  command <<<
    ${Root}/third_party/sambamba view ${finalBAM} | awk '{if ($9>0){print $9}}' | \
    perl -we 'my %hash=();
    while(<>) {
      chomp;
      my $frag = $_;
      if (exists $hash{$frag}) {
        $hash{$frag}++;
      }
      else {
        $hash{$frag}=1;
      }
    }
    open OUT,">${Outdir}/temp/${name}_frags.txt" or die;
    foreach my $key (sort keys %hash){
      print OUT "$key\t$hash{$key}\n";
    }
    close OUT;'
  >>>
  output {
    String fragsTab="${Outdir}/temp/${name}_frags.txt"
  }
}
