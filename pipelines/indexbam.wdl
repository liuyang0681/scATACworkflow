task indexbam {
  String bam
  String Root
  command <<<
    ${Root}/third_party/sambamba index ${bam}
  >>>
  output {
    String final="${bam}"
  }
}
