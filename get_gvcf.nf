/*
Provide directories as inputs with .gvcf files present

nextflow run ~/mygit/nf-haplocaller/get_gvcf.nf --input "0*" --fasta 100.REF_SEQ/TAIR10_wholeGenome.fasta --outdir temp

*/


params.project = "the1001genomes"
build_index = false
params.outdir = './snpcall'
params.fasta = false
params.tmpdir = "/lustre/scratch/users/rahul.pisupati/tempFiles/"


if ( params.fasta ){
  genome = file(params.fasta)
  reffol = genome.parent
  refid = genome.baseName
  if( !genome.exists() ) exit 1, "Reference fasta file not found: ${params.fasta}"
  bwa_indices = Channel
    .fromPath( "$reffol/${refid}.fasta.b*" )
    .ifEmpty { build_index = true }
} else {
  exit 1, "Provide reference fasta file. Ex., --fasta file_path"
}

/*
* 1. Create a channel for checking bwa index for genome ref
*/
if (build_index == true){
  process makeBWAindex {
      publishDir "${reffol}", mode: 'copy'
      label 'env_bwa_small'

      input:
      file genome

      output:
      file "${refid}.fasta.*" into fasta_index
      file "${refid}.dict" into fasta_dict

      script:
      """
      samtools faidx ${genome}
      bwa index $genome
      java -jar \$EBROOTPICARD/picard.jar  CreateSequenceDictionary R=$genome O=${refid}.dict
      """
  }
} else {
  fasta_index = Channel
    .fromPath( "$reffol/${refid}.fasta.*" )
  fasta_dict = Channel
    .fromPath( "$reffol/${refid}.dict" )
}

input_gvcfs = Channel
      .fromPath( "${params.input}", type: 'dir' )
      .map{ it -> [it.name, file("$it/*g.vcf.gz"), file("$it/*g.vcf.gz.tbi")] }

process joinGVCFs {
  tag "$fol_name"
  publishDir "$params.outdir", mode: 'copy'
  label 'env_gatk_large'

  input:
  set val(fol_name), file(in_vcf), file(in_vcf_idx) from input_gvcfs
  file genome
  file indices from fasta_index.collect()
  file fa_dict from fasta_dict.collect()


  output:
  set file("${fol_name}.vcf.gz"), file("${fol_name}.vcf.gz.tbi") into combgVCF

  script:
  def try_vcfs = in_vcf.collect { "-V $it" }.join(' ')
  """
  java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTGATK/GenomeAnalysisTK.jar\
  -T GenotypeGVCFs -R $genome\
  -nt ${task.cpus} \
  $try_vcfs -o ${fol_name}.vcf.gz \
  -allSites -variant_index_type LINEAR -variant_index_parameter 128000
  """
}
