workflow scATAC_standard {
  String Root
  String reference
  String Sample
  String Fastq1
  String Fastq2
  String Outdir
  ## till now, only accept sci/drop/skip mode
  ## sci for single cell index method
  ## drop for droplet method
  ## skip for not parse the cell barcode, works because the cell barcodes are mannully added
  String mode 
  call makedir {
    input: Outdir=Outdir
  }
  call fastq2bam {
    input:
    Outdir=Outdir,
    Fastq1=Fastq1,
    Fastq2=Fastq2,
    reference=reference,
    Root=Root,
    start=makedir.timer,
    mode=mode,
    Sample=Sample,
  }
  call countCB {
    input:
    Outdir=Outdir,
    Root=Root,
    CBtable=fastq2bam.CBtable,
  }
  call sortBAM {
    input:
    Outdir=Outdir,
    Root=Root,
    rawBAM=fastq2bam.rawBAM,
    Sample=Sample,
  }
  call callPeak {
    input:
    Outdir=Outdir,
    Root=Root,
    Sample=Sample,
    finalBAM=sortBAM.finalBAM,
  }
  call peakCount {
    input:
    Outdir=Outdir,
    Root=Root,
    Sample=Sample,
    peak=callPeak.peak,
    finalBAM=sortBAM.finalBAM,
  }

  call usableCount {
    input:
    Outdir=Outdir,
    Root=Root,
    finalBAM=sortBAM.finalBAM,
  }
  call barcodeQC {
    input:
    Root=Root,
    Outdir=Outdir,
    usableCount=usableCount.usableCount,
    peakCount=peakCount.peakCount,
    rawCount=countCB.rawCount,
  }
  call mitoVCF {
    input:
    Outdir=Outdir,
    Root=Root,
    mitoBAM=fastq2bam.mitoBAM,
  }
  call countMTX {
    input:
    Outdir=Outdir,
    Root=Root,
    BClist=barcodeQC.BClist,
    peak=callPeak.peak,
    }
  call runCicero {
    input:
    Outdir=Outdir,
    Root=Root,
    matrix=countMTX.matrix,
    peak=callPeak.peak,
    Sample=Sample,    
  }
  call runChromVAR {
    input:
    Outdir=Outdir,
    Root=Root,
    matrix=countMTX.matrix,
    Sample=Sample,    
  }  
}
task makedir {
  String Outdir
  command {
    mkdir -p ${Outdir}
    mkdir -p ${Outdir}/report
    mkdir -p ${Outdir}/alignment
    mkdir -p ${Outdir}/temp
    mkdir -p ${Outdir}/fastq
    mkdir -p ${Outdir}/report/Cicero
    mkdir -p ${Outdir}/report/chromVAR    
    echo "[`date +%F` `date +%T`] workflow start" >> ${Outdir}/workflowtime.log
  }
  output {
    String timer="`perl -we 'print time()'`"    
  }
}
task fastq2bam {
  String Outdir
  String Fastq1
  String Fastq2
  String mode
  String Sample
  String reference
  String Root
  String start
  command {
    ${Root}/bin/sc_atac_parse_cellbarcode -t 10 ${Fastq1} ${Fastq2} -platform ${mode} -report ${Outdir}/report/cb_report.md -cb ${Outdir}/temp/barcode_table.txt| \
    ${Root}/bin/parse_Tn5me_adaptors -t 10 -d -p -report ${Outdir}/report/trim_me_report.md |\
    ${Root}/third_party/bwa mem -t 20 -p ${reference} /dev/stdin | \
    ${Root}/bin/parse_alignment_readname -filter ${Outdir}/temp/filter.bam -o ${Outdir}/temp/aln.bam -report ${Outdir}/report/aln.md -t 5 -maln ${Outdir}/temp/mito.bam
  }
  output {
    String CBtable="${Outdir}/temp/barcode_table.txt"
    #String timer="`perl -we 'print time()'`"
    String mitoBAM="${Outdir}/temp/mito.bam"
    String rawBAM="${Outdir}/temp/aln.bam"
  }
}
task countCB {
  String Outdir
  String CBtable
  String Root
  command <<<
    perl -we 'open BT,"${Outdir}/temp/barcode_table.txt"; open OT,">","${Outdir}/temp/raw_cell_count.txt"; my %hash=(); while(<BT>) {my ($key, $qvalue) = split/\t/; if (exists $hash{$key}){$hash{$key}++;}else{$hash{$key}=1;}} foreach my $key (sort keys %hash){ if($hash{$key}>100){print OT "$key\t$hash{$key}\n";}} close BT; close OT;'
  >>>
  output {
    #String timer="`perl -we 'print time()'`"
    String rawCount="${Outdir}/temp/raw_cell_count.txt"
  }  
}
task sortBAM {
  String Outdir
  String Root
  #String start
  String rawBAM
  String Sample
  command <<<
    ${Root}/third_party/sambamba sort -t 24 -o ${Outdir}/temp/aln_sorted.bam ${rawBAM}
    
    java -jar ${Root}/third_party/picard.jar MarkDuplicates I=${Outdir}/temp/aln_sorted.bam O=${Outdir}/alignment/aln_rmdup.bam M=${Outdir}/report/duplicates.matrix Barcode_tag=CR REMOVE_DUPLICATES=true

    ${Root}/third_party/sambamba index ${Outdir}/alignment/aln_rmdup.bam

    #rm -f ${Outdir}/alignment/aln.bam ${Outdir}/alignment/aln_sorted.bam ${Outdir}/alignment/aln_sorted.bam.bai
    
    echo "[`date +%F` `date +%T`] human reference alignment finish" >> ${Outdir}/workflowtime.log
  >>>
  output {
    #String timer="`perl -we 'print time()'`"
    String finalBAM="${Outdir}/alignment/aln_rmdup.bam"
  }
}
task callPeak {
  String Outdir
  String Root
  String Sample
  String finalBAM
  command <<<
    ${Root}/third_party/macs2 callpeak -t ${finalBAM} -f BAM -n ${Sample} -q 0.01 --nomodel --shift -75 --extsize 150 --outdir ${Outdir}/temp/

    cut -f1,2,3 ${Outdir}/temp/${Sample}_peaks.narrowPeak | grep -v "_" > ${Outdir}/report/${Sample}_peak.bed
  >>>
  output {
    String peak="${Outdir}/report/${Sample}_peak.bed"
  }
}
task peakCount {
  String Outdir
  String Root
  String Sample
  String peak
  String finalBAM
  command <<<
    ${Root}/third_party/sambamba view ${finalBAM} | sed -r "s|.*CR:Z:([A-Z]*).*|\1\t\0|" | awk '{if ($10>0) {printf("%s\t%d\t%d\t%s\n",$4,$5,$5+$10-1,$1)}}' > ${Outdir}/temp/aln.bed

    ${Root}/third_party/bedtools intersect -a ${peak} -b ${Outdir}/temp/aln.bed -wa -wb | cut -f7 | perl -we 'open OT,">","${Outdir}/temp/peak_cell_count.txt"; my %hash=(); while(<>) {my $key = $_; chomp $key; if (exists $hash{$key}){$hash{$key}++;}else{$hash{$key}=1;}} foreach my $key (sort keys %hash){ print OT "$key\t$hash{$key}\n";} close OT;'
  >>>
  output {
    String peakCount="${Outdir}/temp/peak_cell_count.txt"
  }
}
task usableCount {
  String Outdir
  String Root
  String finalBAM
  command <<<
    ${Root}/third_party/sambamba view ${finalBAM} | sed -r "s|.*CR:Z:([A-Z]*).*|\1|" | perl -we 'open OT,">","${Outdir}/temp/usable_cell_count.txt"; my %hash=(); while(<>) {my $key = $_; chomp $key; if (exists $hash{$key}){$hash{$key}++;}else{$hash{$key}=1;}} foreach my $key (sort keys %hash){ if($hash{$key}){print OT "$key\t$hash{$key}\n";}} close OT;'
  >>>
  output {
    String usableCount="${Outdir}/temp/usable_cell_count.txt"
  }
}
task barcodeQC {
  String Root
  String Outdir
  String usableCount
  String peakCount
  String rawCount
  command <<<
    perl -we 'my %hash = {}; open BC,">", "${Outdir}/temp/bc_list.txt"; open OUT,">","${Outdir}/temp/result_barcode_counts.txt"; open RAW,"${rawCount}" or die "Cannot open ${rawCount}"; open USABLE,"${usableCount}" or die "Cannot open ${usableCount}"; open PEAK,"${peakCount}" or die "${peakCount}"; while(<RAW>) { chomp; my ($key, $val) = split; if (exists $hash{$key}) {$hash{$key}[0]+=$val;} elsif ($key) { my @val = ($val,0,0); $hash{$key} = \@val;} } while(<USABLE>) { chomp; my ($key,$val) = split; if ($key && exists $hash{$key}) { $hash{$key}[1]+=$val;}} while(<PEAK>){ chomp; my ($key,$val) = split; if ($key && exists $hash{$key}){$hash{$key}[2]+=$val;}} foreach my $key (keys %hash) { print OUT "$key\t$hash{$key}[0]\t$hash{$key}[1]\t$hash{$key}[2]\n"; if ($hash{$key}[1]>200 && $hash{$kvey}[2]>100){print BC "$key\n"}} close OUT; close BC; close RAW; close USABLE; close RAW; close PEAK;'
  >>>
  output {
    String BClist="${Outdir}/temp/bc_list.txt"  
  }
}
task countMTX {
  String Outdir
  String BClist
  String Root
  String peak
  command <<<    
    python ${Root}/scripts/sc_atac_count.py ${Outdir}/alignment/aln_rmdup.bam ${BClist} ${peak} ${Outdir}/report/count.mtx True
  >>>
  output {
    String matrix="${Outdir}/report/count.mtx"
  }
}
task mitoVCF {
  String Outdir
  String Root
  String mitoBAM
  command{
    ${Root}/third_party/sambamba sort -t 24 -o ${Outdir}/temp/mito_sorted.bam ${mitoBAM}
    java -jar ${Root}/third_party/picard.jar MarkDuplicates I=${Outdir}/temp/mito_sorted.bam O=${Outdir}/alignment/mito_rmdup.bam M=${Outdir}/report/duplicates_mito.matrix Barcode_tag=CR REMOVE_DUPLICATES=true
  }
  output {
    String timer="`perl -we 'print time()'`"
  }
}
task runCicero {
  String Outdir
  String Root
  String peak
  String matrix
  String Sample
  command {
    ${Root}/third_party/Rscript ${Root}/scripts/run_cicero.R -P ${peak} -C ${matrix} -O ${Outdir}/report/Cicero -S ${Sample}
  }
}
task runChromVAR {
  String Outdir
  String Root
  String matrix
  String Sample
  command {
    ${Root}/third_party/Rscript ${Root}/scripts/run_chromVAR.R ${matrix} ${Outdir}/report/chromVAR/${Sample}
  }
}
