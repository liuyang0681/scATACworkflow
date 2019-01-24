import "fastq2bam.wdl" as f2b
import "sortBam.wdl" as srt
import "countCB.wdl" as countCB
import "usableCount.wdl" as usCB
import "scATAC_report.wdl" as report
import "countMTX.wdl" as countMTX
import "atac_bam2peak.wdl" as bam2peak
import "countFrags.wdl" as countFrag
import "accessEnhancer.wdl" as accEnh
import "indexbam.wdl" as indexbam

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
    String Outdir="${Dir}"
  }
}

workflow core {
  String Root
  String TestID
  String FQ1
  String FQ2
  String Outdir
  String mode
  String gsize
  String reference
  String tssbed
  String tssflank
  String? lane
  call makedir {
    input:
    Dir=Outdir
  }
  call f2b.fastq2bam as aln{
    input:
    Root=Root,
    Fastq1=FQ1,
    Fastq2=FQ2,
    mode=mode,
    Sample=TestID,
    reference=reference,
    Outdir=makedir.Outdir,
    lane=lane,
  }
  call srt.sortBAM as sort {
    input:
    Root=Root,
    bamfile=aln.rawBAM,
    Outdir=Outdir,
    name=TestID,
  }
  call countCB.count as count {
    input:
    Root=Root,
    CBtable=aln.CBtable,
    Outdir=Outdir
  }
  call usCB.usableCB as usCB {
    input:
    Root=Root,
    Outdir=Outdir,
    TSS=tssflank,
    finalBAM=sort.final,
    name=TestID,
  }
  call barcodeQC {
    input:
    Root=Root,
    Outdir=Outdir,
    readsCount=usCB.readsCount,
    rawCount=count.rawCount,
  }
  call bam2peak.callpeak as callpeak {
    input:
    finalBAM=sort.final,
    gsize=gsize,
    Outdir=Outdir,
    Root=Root,
    Sample=TestID,
  }
  call indexbam.indexbam as idxBam {
    input:
    Root=Root,
    bam=sort.final,
  }
  call countFrag.countFrags as cntfg {
    input:
    Root=Root,
    Outdir=Outdir,
    name=TestID,
    finalBAM=sort.final,
  }
  call accEnh.accessEnhancer as accEnh {
    input:
    Root=Root,
    Outdir=Outdir,
    name=TestID,
    finalBAM=idxBam.final,
    TssBed=tssbed,
  }
  call countMTX.countMTX as cntMTX {
    input:
    Outdir=Outdir,
    Root=Root,
    BClist=barcodeQC.BClist,
    peak=callpeak.peak,
  }
  call report.report {
    input:
    Outdir=Outdir,
    Root=Root,
    matrix=cntMTX.matrix,
    Sample=TestID,
    accEnhPlot=accEnh.accEnhPlot,
    frag=cntfg.fragsTab,
    CBreport=aln.CBreport,
    TMreport=aln.TMreport,
    ALNreport=aln.ALNreport,
  }
}
task barcodeQC {
  String Root
  String Outdir
  String readsCount
  String rawCount
  command <<<
    perl -we 'my %hash = ();
    open BC,">", "${Outdir}/temp/bc_list.txt";
    open OUT,">","${Outdir}/report/result_barcode_counts.txt";
    open RAW,"${rawCount}" or die;
    open USABLE,"${readsCount}" or die;
    while(<RAW>) {
      chomp;
      my ($key, $val) = split;
      if (exists $hash{$key}) {
        $hash{$key}[0]+=$val;
      }
      elsif ($key) {
        my @val = ($val,0,0);
        $hash{$key} = \@val;
      }
    }
    while(<USABLE>) {
      chomp;
      my ($key,$val1,$val2) = split;
      if ($key && exists $hash{$key}) {
        $hash{$key}[1]+=$val1;
        $hash{$key}[2]+=$val2;
      }
    }
    print OUT "barcode_ID\tbarcode_reads\tusable_reads\treads_in_promoters\n";
    foreach my $key (keys %hash) {
      print OUT "$key\t$hash{$key}[0]\t$hash{$key}[1]\t$hash{$key}[2]\n";
      if ($hash{$key}[1]>200 && $hash{$key}[2]>100){
        print BC "$key\n";
      }
    }
    close OUT;
    close BC;
    close RAW;
    close USABLE;
    close RAW;'    
  >>>
  output {
    String BClist="${Outdir}/temp/bc_list.txt"
    String result="${Outdir}/report/result_barcode_counts.txt"
  }
}
