task usableCB {
  String Root
  String Outdir
  String TSS
  String name
  String finalBAM
  command <<<
    ${Root}/third_party/sambamba view ${finalBAM} | sed -r "s|([A-Z0-9_]*)\t.*CR:Z:([A-Z]*)\t.*|\1\t\2|" > ${Outdir}/temp/${name}_ID.txt
  
    ${Root}/third_party/bedtools intersect -a ${TSS} -b ${finalBAM} -wb | cut -f7 |sed -r "s|(.*)/.|\1|" > ${Outdir}/temp/${name}_TSS_count.txt
  
    perl -we 'open ID,"${Outdir}/temp/${name}_ID.txt" or die; open TSS, "${Outdir}/temp/${name}_TSS_count.txt" or die; my %hash=(); my %cb = (); while(<ID>) { chomp; my ($id, $cb) = split; $hash{$id} = $cb; if ( exists $cb{$cb} ){$cb{$cb}[0]++;} else { my @val=(1,0); $cb{$cb}=\@val;}} close ID; while(<TSS>){ chomp; my $id=$_; if(!exists $hash{$id}) { die "Unporperly read name.";} else { my $cb = $hash{$id}; $cb{$cb}[1]++;}} close TSS; open OUT,">","${Outdir}/temp/${name}_readcounts.txt" or die; foreach my $key (keys %cb) { print OUT "$key\t$cb{$key}[0]\t$cb{$key}[1]\n";}'
  >>>
  output {
    String readsCount="${Outdir}/temp/${name}_readcounts.txt"
  }
}
workflow usableCount {
  String Root
  String Outdir
  String name
  String TSS
  String finalBAM
  call usableCB {
    input:
    Root=Root,
    Outdir=Outdir,
    name=name,
    TSS=TSS,
    finalBAM=finalBAM,
  }
  output {
    String readsCount=usableCB.readsCount
  }
}
