workflow scATAC_benchmark_hs_mm {
  String Root
  String hs_reference
  String mm_reference
  String Sample
  String Fastq1
  String Fastq2
  String Outdir
  call makedir {
    input:
    Dir=Outdir
  }
  call split_barcode {
    input:
    Root=Root,
    Fastq1=Fastq1,
    Fastq2 = Fastq2,
    Outdir=Outdir,
    flag=makedir.success
  }
  call countCB {
    input:
    Outdir=Outdir,
    CBtable=split,
    Root=Root,
  }
  call align_hs {
    input:
    Root=Root,
    reference=hs_reference,
    Outdir=Outdir,
    flag=split_barcode.finished
  }
  call align_mm {
    input:
    Root=Root,
    reference=mm_reference,
    Outdir=Outdir,
    flag=split_barcode.finished
  }
  call report {
    input:
    Root=Root,
    Outdir=Outdir,
    Sample=Sample,
    flag1=align_hs.finished,
    flag2=align_mm.finished,
    start=makedir.start,
  }
}
task makedir {
  String Dir
  command {
    echo "[`date +%F` `date +%T`] workflow start" > ${Dir}/workflowtime.log
    mkdir -p ${Dir}
    mkdir -p ${Dir}/report
    mkdir -p ${Dir}/alignment
    mkdir -p ${Dir}/temp
    mkdir -p ${Dir}/fastq
  }
  output {
    String start=`perl -we "print time()"`
  }
}
task split_barcode {
  String flag
  String Root
  String Outdir
  String Fastq1
  String Fastq2
  command {
    #$${Root}/bin/sc_atac_parse_cellbarcode -5 ${Fastq1} ${Fastq2} -platform stream -cb ${}
    ${Root}/bin/parse_barcode_ME_adaptor -d -t 5 -b ${Outdir}/report/barcode_table.txt -report ${Outdir}/report/split_report.txt ${Fastq1} ${Fastq2} 1> ${Outdir}/fastq/clean.fq 2>> ${Outdir}/workflowtime.log
  }
  output {
    String CBtable="${Outdir}/report/barcode_table.txt"
  }
}
task countCB {
  String Outdir
  String CBtable
  String Root
  command <<<
    perl -we 'open BT,"${CBtable}"; open OT,">","${Outdir}/temp/raw_cell_count.txt"; my %hash=(); while(<BT>) {my ($key, $qvalue) = split/\t/; if (exists $hash{$key}){$hash{$key}++;}else{$hash{$key}=1;}} foreach my $key (sort keys %hash){ if($hash{$key}>100){print OT "$key\t$hash{$key}\n";}} close BT; close OT;'
  >>>
  output {
    #String timer="`perl -we 'print time()'`"
    String rawCount="${Outdir}/temp/raw_cell_count.txt"
  }  
}

task align_hs {
  String Root
  String reference
  String flag
  String Outdir
  command {
    ${Root}/bin/bwa mem -t 20 -p ${reference} ${Outdir}/fastq/clean.fq |\
    ${Root}/bin/parse_alignment_readname -t 4 -o ${Outdir}/alignment/hs.bam -filter ${Outdir}/alignment/hs_filter.bam -report ${Outdir}/report/hs_alignment_report.txt /dev/stdin
    ${Root}/bin/sambamba sort -t 24 -o ${Outdir}/alignment/hs_sorted.bam ${Outdir}/alignment/hs.bam 
    java -jar ${Root}/bin/picard.jar MarkDuplicates I=${Outdir}/alignment/hs_sorted.bam O=${Outdir}/alignment/hs_rmdup.bam M=${Outdir}/report/hs_dup.matrix Barcode_tag=CR REMOVE_DUPLICATES=true
    ${Root}/bin/samtools view ${Outdir}/alignment/hs_rmdup.bam | sed -r "s|([A-Z0-9_]*)\t.*CR:Z:([A-Z]*)\t.*|\1\t\2|" > ${Outdir}/report/ID_hs.txt
    ${Root}/bin/bedtools intersect -a ${Root}/file/merge_sort_2K_hg19_gencode_tss_unique.bed -b ${Outdir}/alignment/hs_rmdup.bam -wb |cut -f7|sed -r "s|(.*)/.|\1|" > ${Outdir}/report/TSS_count_hs.txt
    rm -f ${Outdir}/alignment/hs_sorted.bam ${Outdir}/alignment/hs_sorted.bam.bai ${Outdir}/alignment/hs.bam
    echo "[`date +%F` `date +%T`] human reference alignment finish" >> ${Outdir}/workflowtime.log
  }
  output {
    String success="success"
  }
}
task align_mm {
  String Root
  String reference
  String Outdir
  String flag
  command {
    ${Root}/bin/bwa mem -t 20 -p ${reference} ${Outdir}/fastq/clean.fq|\
    ${Root}/bin/parse_alignment_readname -t 4 -o ${Outdir}/alignment/mm.bam -filter ${Outdir}/alignment/mm_filter.bam -report ${Outdir}/report/mm_alignment_report.txt /dev/stdin 2>>${Outdir}/workflowtime.log
    ${Root}/bin/sambamba sort -t 24 -o ${Outdir}/alignment/mm_sorted.bam ${Outdir}/alignment/mm.bam 
    java -jar ${Root}/bin/picard.jar MarkDuplicates I=${Outdir}/alignment/mm_sorted.bam O=${Outdir}/alignment/mm_rmdup.bam M=${Outdir}/report/mm_dup.matrix BARCODE_TAG=CR REMOVE_DUPLICATES=true
    ${Root}/bin/samtools view ${Outdir}/alignment/mm_rmdup.bam | sed -r "s|([A-Z0-9_]*)\t.*CR:Z:([A-Z]*)\t.*|\1\t\2|" > ${Outdir}/report/ID_mm.txt	
    ${Root}/bin/bedtools intersect -a ${Root}/file/merge_sort_2K_mm9_gencode_tss_unique.bed -b ${Outdir}/alignment/mm_rmdup.bam -wb |cut -f7|sed -r "s|(.*)/.|\1|" > ${Outdir}/report/TSS_count_mm.txt
    rm -f ${Outdir}/alignment/mm_sorted.bam ${Outdir}/alignment/mm_sorted.bam.bai ${Outdir}/alignment/mm.bam
    echo "[`date +%F` `date +%T`] mouse reference alignment finish" >> ${Outdir}/workflowtime.log
  }
  output {
    String success="success"
  }
}
task report {
  String Outdir
  String Root
  String Sample
  String rawCount
  String flag1
  String flag2
  command <<<
    Rscript ${Root}/script/QC_V1.R -C ${rawCount} -H_r ${Outdir}/report/ID_hs.txt -H_t ${Outdir}/report/TSS_count_hs.txt -M_r ${Outdir}/report/ID_mm.txt -M_t ${Outdir}/report/TSS_count_mm.txt -E 1000 -R 0.8 -U 500 -T 0.1 -O ${Outdir}/report -B ${Sample} && echo "[`date +%F` `date +%T`] workflow end" >> ${Outdir}/workflowtime.log
    if [ $? ]
    then
    echo "[`date +%F` `date +%T`] failure" >> ${Outdir}/workflowtime.log
  >>>
  output {
    String success="success"
  }
}
