workflow hs_mm_384well {
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
    start=makedir.start,
  }
  call countCB {
    input:
    Root=Root,
    CBtable=split_barcode.CBtable,
    Outdir=Outdir,
  }
  call align_hs {
    input:
    Root=Root,
    reference=hs_reference,
    Outdir=Outdir,
    FQ1=split_barcode.FQ1,
    FQ2=split_barcode.FQ2,
  }
  call align_mm {
    input:
    Root=Root,
    reference=mm_reference,
    Outdir=Outdir,
    FQ1=split_barcode.FQ1,
    FQ2=split_barcode.FQ2,
  }
  call report {
    input:
    Root=Root,
    Outdir=Outdir,
    rawCount=countCB.rawCount,
    Sample=Sample,
    hs_ID=align_hs.ID,
    mm_ID=align_mm.ID,
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
    String start="success"
  }
}
task split_barcode {
  String Outdir
  String Fastq1
  String Fastq2
  String Root
  String start
  command {
    ${Root}/bin/sc_atac_parse_cellbarcode -t 10 ${Fastq1} ${Fastq2} -platform stream -report ${Outdir}/report/cb_report.md -cb ${Outdir}/temp/barcode_table.txt| \
    ${Root}/bin/parse_Tn5me_adaptors -t 10 -d -p -report ${Outdir}/report/trim_me_report.md -1 ${Outdir}/fastq/read_1.fq -2 ${Outdir}/fastq/read_2.fq
  }
  output {
    String CBtable="${Outdir}/temp/barcode_table.txt"
    String FQ1="${Outdir}/fastq/read_1.fq"
    String FQ2="${Outdir}/fastq/read_2.fq"
  }
}
task countCB {
  String Outdir
  String CBtable
  String Root
  command <<<
    perl -we 'open BT,"${CBtable}"; open OT,">","${Outdir}/temp/raw_cell_count.txt"; my %hash=(); while(<BT>) {chomp; my ($key, $qvalue) = split/\t/; if (exists $hash{$key}){$hash{$key}++;}else{$hash{$key}=1;}} foreach my $key (sort keys %hash){ if($hash{$key}>100){print OT "$key\t$hash{$key}\n";}} close BT; close OT;'
  >>>
  output {
    String rawCount="${Outdir}/temp/raw_cell_count.txt"
  }  
}
task align_hs {
  String Root
  String reference
  String Outdir
  String FQ1
  String FQ2
  command {
    ${Root}/third_party/bwa mem -t 20 -p ${reference} ${FQ1} ${FQ2} | \
    ${Root}/bin/parse_alignment_readname -filter ${Outdir}/temp/hs_filter.bam -o ${Outdir}/temp/hs_aln.bam -report ${Outdir}/report/hs_aln.md -t 5 -maln ${Outdir}/temp/hs_mito.bam
    ${Root}/third_party/sambamba sort -t 24 -o ${Outdir}/temp/hs_sorted.bam ${Outdir}/temp/hs_aln.bam 
    java -jar ${Root}/third_party/picard.jar MarkDuplicates I=${Outdir}/temp/hs_sorted.bam O=${Outdir}/alignment/hs_rmdup.bam M=${Outdir}/report/hs_dup.matrix Barcode_tag=CR REMOVE_DUPLICATES=true
    ${Root}/third_party/sambamba view ${Outdir}/alignment/hs_rmdup.bam | sed -r "s|([A-Z0-9_]*)\t.*CR:Z:([A-Z]*)\t.*|\1\t\2|" > ${Outdir}/report/ID_hs.txt
    ${Root}/third_party/bedtools intersect -a ${Root}/file/merge_sort_2K_hg19_gencode_tss_unique.bed -b ${Outdir}/alignment/hs_rmdup.bam -wb |cut -f7|sed -r "s|(.*)/.|\1|" > ${Outdir}/report/TSS_count_hs.txt
    echo "[`date +%F` `date +%T`] human reference alignment finish" >> ${Outdir}/workflowtime.log
  }
  output {
    String ID = "${Outdir}/report/ID_hs.txt"
  }
}
task align_mm {
  String Root
  String reference
  String Outdir
  String FQ1
  String FQ2
  command {
    ${Root}/third_party/bwa mem -t 20 -p ${reference} ${FQ1} ${FQ2} | \
    ${Root}/bin/parse_alignment_readname -filter ${Outdir}/temp/mm_filter.bam -o ${Outdir}/temp/mm_aln.bam -report ${Outdir}/report/mm_aln.md -t 5 -maln ${Outdir}/temp/mm_mito.bam
    ${Root}/third_party/sambamba sort -t 24 -o ${Outdir}/temp/mm_sorted.bam ${Outdir}/temp/mm_aln.bam 
    java -jar ${Root}/third_party/picard.jar MarkDuplicates I=${Outdir}/temp/mm_sorted.bam O=${Outdir}/alignment/mm_rmdup.bam M=${Outdir}/report/mm_dup.matrix Barcode_tag=CR REMOVE_DUPLICATES=true
    ${Root}/third_party/sambamba view ${Outdir}/alignment/mm_rmdup.bam | sed -r "s|([A-Z0-9_]*)\t.*CR:Z:([A-Z]*)\t.*|\1\t\2|" > ${Outdir}/report/ID_mm.txt
    ${Root}/third_party/bedtools intersect -a ${Root}/file/merge_sort_2K_mm9_gencode_tss_unique.bed -b ${Outdir}/alignment/mm_rmdup.bam -wb |cut -f7|sed -r "s|(.*)/.|\1|" > ${Outdir}/report/TSS_count_mm.txt
    echo "[`date +%F` `date +%T`] human reference alignment finish" >> ${Outdir}/workflowtime.log
  }
  output {
    String ID = "${Outdir}/report/ID_mm.txt"
  }
}
task report {
  String Outdir
  String Root
  String Sample
  String hs_ID
  String mm_ID
  String rawCount
  command <<<
    ${Root}/third_party/Rscript ${Root}/scripts/QC_V1.R -C ${rawCount} -H_r ${hs_ID} -H_t ${Outdir}/report/TSS_count_hs.txt -M_r ${mm_ID} -M_t ${Outdir}/report/TSS_count_mm.txt -E 1000 -R 0.8 -U 500 -T 0.1 -O ${Outdir}/report -B ${Sample}
    echo "[`date +%F` `date +%T`] workflow end" >> ${Outdir}/workflowtime.log
  >>>
  output {
    String success="success"
  }
}
