/*
provide all the GVCFs together, the script can split the samples separately


nextflow run get_filter_vcf.nf --reads "001.plate.raw.vcfs/*vcf" --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
*/


params.project = "the1001genomes"
build_index = false
params.outdir = './snpcall'
params.fasta = false
params.tmpdir = "/lustre/scratch/users/rahul.pisupati/tempFiles/"

params.name = false

if ( params.fasta ){
  genome = file(params.fasta)
  reffol = genome.parent
  refid = genome.baseName
  if( !genome.exists() ) exit 1, "Reference fasta file not found: ${params.fasta}"
  bwa_indices = Channel
    .fromPath( "$reffol/${refid}.fasta.b*" )
    .ifEmpty { build_index = true }
    .subscribe onComplete: { checked_genome_index = true }
} else {
  exit 1, "Provide reference fasta file. Ex., --fasta file_path"
}

input_gvcfs = Channel
      .fromPath( "${params.reads}" )
      .map { it -> [ file("$it"), file("$it" + ".idx") ] }

process getSamples {
  tag "joinGVCF"

  input:
  set file(gvcf), file(gvcf_index) from input_gvcfs

  output:
  set file("$gvcf"), stdout into sample_names

  script:
  """
  grep -m 1 "^#CHROM" $gvcf | awk '{for(i=10;i<=NF;i++) print \$i}'
  """
}

input_names = sample_names
    .splitText(elem: 1)
    .map{ row -> [file("${row[0]}"), file("${row[0]}.idx"), "${row[1]}".replace('\n','') ]  }

process selectSNPs {
  tag "$name"
  publishDir "${params.outdir}/$out_fol", mode: 'copy'

  input:
  set file(gvcf), file(gvcf_index), val(name) from input_names

  when:
  file("${params.outdir}/${out_fol}/${name}.filter.vcf").isEmpty()

  output:
  set file("${name}.filter.vcf"), file("${name}.filter.vcf.idx") into filter_vcf

  script:
  out_fol =  file("$gvcf").baseName.replace(/.12Jan16.all.vcf/, "")
  """
  java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTGATK/GenomeAnalysisTK.jar\
  -T SelectVariants -R $reffol/${refid}.fasta\
  -nt ${task.cpus} \
  -V $gvcf \
  -o ${name}.filter.vcf \
  -se "${name}" \
  -selectType SNP -restrictAllelesTo BIALLELIC
  """
}
