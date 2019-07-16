/*
provide all the GVCFs together, the script can split the samples separately


nextflow run get_filter_vcf.nf --input "001.plate.raw.vcfs/*vcf" --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
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

input_gvcfs = Channel
      .fromPath( "${params.input}" )
      .map { it -> [ file("$it").getExtension(), file("$it"), file("${it}.*") ] }

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

process getSamples {
  tag "$gvcf"
  label 'env_gatk_small'

  input:
  set val(file_ext), file(gvcf), file(gvcf_index) from input_gvcfs

  output:
  set val(out_fol), file("$gvcf"), file("$gvcf_index"), stdout into sample_names

  script:
  out_fol =  file("$gvcf").baseName.replace( /.vcf/, '')
  if (file_ext == "gz"){
    """
    tabix -H $gvcf | grep -m 1 "^#CHROM" | awk '{for(i=10;i<=NF;i++) print \$i}'
    """
  } else if (file_ext == "vcf"){
    """
    grep -m 1 "^#CHROM" $gvcf | awk '{for(i=10;i<=NF;i++) print \$i}'
    """
  }
}

input_names = sample_names
    .splitText(elem: 3)
    .map{ row -> ["${row[0]}", file("${row[1]}"), file("${row[2]}"), "${row[3]}".replace('\n','') ]  }

process selectSNPs {
  tag "$name"
  publishDir "${params.outdir}/$out_fol", mode: 'copy'
  label 'env_gatk_small'

  input:
  set val(out_fol), file(gvcf), file(gvcf_index), val(name) from input_names

  when:
  file("${params.outdir}/${out_fol}/${name}.filter.vcf.idx").isEmpty()

  output:
  set file("${name}.filter.vcf"), file("${name}.filter.vcf.idx") into filter_vcf

  script:
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
