/*
Provide directories as inputs



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
      .fromPath( "${params.input}", type: 'dir' )
      .map{ it -> [it.name, file("$it/*gvcf")] }

process joinGVCFs {
  tag "$fol_name"
  publishDir "$params.outdir", mode: 'copy'

  input:
  set val(fol_name), file(in_vcf) from input_gvcfs

  output:
  set file("${fol_name}.vcf.gz"), file("${fol_name}.vcf.gz.tbi") into combgVCF

  script:
  def try_vcfs = in_vcf.collect { "-V $it" }.join(' ')
  """
  java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTGATK/GenomeAnalysisTK.jar\
  -T GenotypeGVCFs -R $reffol/${refid}.fasta\
  -nt ${task.cpus} \
  $try_vcfs -o ${fol_name}.vcf.gz \
  -allSites -variant_index_type LINEAR -variant_index_parameter 128000
  """
}
