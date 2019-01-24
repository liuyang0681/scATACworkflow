task count {
  String Outdir
  String CBtable
  String Root
  command <<<
    perl -we 'open BT,"${CBtable}"; open OT,">","${Outdir}/temp/raw_cell_count.txt"; my %hash=(); while(<BT>) {chomp; my ($key, $qvalue) = split/\t/; if (exists $hash{$key}){$hash{$key}++;}else{$hash{$key}=1;}} foreach my $key (sort keys %hash){ if($hash{$key}>0){print OT "$key\t$hash{$key}\n";}} close BT; close OT;'
  >>>
  output {
    String rawCount="${Outdir}/temp/raw_cell_count.txt"
  }  
}

workflow countCB {
  String Outdir
  String CBtable
  String Root
  call count {
    input:
    Outdir=Outdir,
    CBtable=CBtable,
    Root=Root,
  }
  output {
    String rawCount=count.rawCount
  }
}
