import "splitCB.wdl" as splitCB
import "alignment.wdl" as align
import "usableCount.wdl" as usCB
import "hs_mm_mixed_report.wdl" as rpt
import "countCB.wdl" as countCB
import "sortBam.wdl" as srt
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

workflow hs_mm_mixed_core {
  String Root
  String TestID
  String Fastq1
  String Fastq2
  String Outdir
  String mode
  call makedir {
    input:
    Dir=Outdir
  }
  call splitCB.split as split{
    input:
    Root=Root,
    Fastq1=Fastq1,
    Fastq2=Fastq2,
    mode=mode,
    Outdir=makedir.Outdir,
  }
  call countCB.count as count {
    input:
    Root=Root,
    CBtable=split.CBtable,
    Outdir=Outdir,
  }
  call align.align as align_hs {
    input:
    Root=Root,
    FQ1=split.FQ1,
    FQ2=split.FQ2,
    Outdir=Outdir,
    name="human",
    reference="${Root}/reference/hg19/hg19.fa",
  }
  call align.align as align_mm {
    input:
    Root=Root,
    FQ1=split.FQ1,
    FQ2=split.FQ2,
    Outdir=Outdir,
    name="mouse",
    reference="${Root}/reference/mm9/mm9.fa",    
  }
  call srt.sortBAM as sort_hs {
    input:
    Outdir=Outdir,
    Root=Root,
    name="human",
    bamfile=align_hs.bamfile,
  }
  call srt.sortBAM as sort_mm {
    input:
    Outdir=Outdir,
    Root=Root,
    name="mouse",
    bamfile=align_mm.bamfile,
  }
  call indexbam.indexbam as idx_hs {
    input:
    Root=Root,
    bam=sort_hs.final,
  }
  call indexbam.indexbam as idx_mm {
    input:
    Root=Root,
    bam=sort_mm.final,
  }
  call countFrag.countFrags as frag_hs {
    input:
    Root=Root,
    Outdir=Outdir,
    name="human",
    finalBAM=sort_hs.final,
  }
  call countFrag.countFrags as frag_mm {
    input:
    Root=Root,
    Outdir=Outdir,
    name="mouse",
    finalBAM=sort_mm.final,
  }
  call usCB.usableCB as usCB_hs {
    input:
    Root=Root,
    Outdir=Outdir,
    name="human",
    TSS="${Root}/file/merge_sort_2K_hg19_gencode_tss_unique.bed",
    finalBAM=sort_hs.final,
  }
  call usCB.usableCB as usCB_mm {
    input:
    Root=Root,
    Outdir=Outdir,
    name="mouse",
    TSS="${Root}/file/merge_sort_2K_mm9_gencode_tss_unique.bed",
    finalBAM=sort_mm.final,
  }
  call accEnh.accessEnhancer as accEnh_hs {
    input:
    Root=Root,
    Outdir=Outdir,
    name="human",
    finalBAM=idx_hs.final,
    TssBed="${Root}/file/hg19_gencode_tss_unique.bed",
  }
  call accEnh.accessEnhancer as accEnh_mm {
    input:
    Root=Root,
    Outdir=Outdir,
    name="mouse",
    finalBAM=idx_mm.final,
    TssBed="${Root}/file/mm9_gencode_tss_unique.bed",
  }
  call rpt.report  {
    input:
    Root=Root,
    Outdir=Outdir,
    RawCount=count.rawCount,
    TestID=TestID,
    mm_UsableCount=usCB_mm.readsCount,
    hs_UsableCount=usCB_hs.readsCount,
    CBreport=split.CBreport,
    TMreport=split.TMreport,
    HsAln=align_hs.aln_report,
    MmAln=align_mm.aln_report,
    TestID=TestID,
    accEnh_hs_plot=accEnh_hs.accEnhPlot,
    accEnh_mm_plot=accEnh_mm.accEnhPlot,
    fragHs=frag_hs.fragsTab,
    fragMm=frag_mm.fragsTab,
  }

}
