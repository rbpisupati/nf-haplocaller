/*
provide all the GVCFs together, the script can split the samples separately


nextflow run get_filter_vcf.nf --reads "001.plate.raw.vcfs/*vcf" --fasta ref_seq/TAIR10_wholeGenome.fasta  --outdir filter_vcfs
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
      .map { it -> [ file("$it").getExtension(), file("$it"), file("$it" + ".idx") ] }

/*
* 1. Create a channel for checking bwa index for genome ref
*/
if (build_index == true){
  process makeBWAindex {
      publishDir "${reffol}", mode: 'copy'

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
}

process getSamples {

  input:
  set val(file_ext), file(gvcf), file(gvcf_index) from input_gvcfs

  output:
  set file("$gvcf"), stdout into sample_names

  script:
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
    .splitText(elem: 1)
    .map{ row -> [file("${row[0]}"), file("${row[0]}.idx"), "${row[1]}".replace('\n','') ]  }

process selectSNPs {
  tag "$name"
  publishDir "${params.outdir}/$out_fol", mode: 'copy'

  input:
  set file(gvcf), file(gvcf_index), val(name) from input_names
  file fasta_index

  when:
  file("${params.outdir}/${out_fol}/${name}.filter.vcf").isEmpty()

  output:
  set file("${name}.filter.vcf"), file("${name}.filter.vcf.idx") into filter_vcf

  script:
  out_fol =  file("$gvcf").baseName.replace( /.vcf/, '')
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
