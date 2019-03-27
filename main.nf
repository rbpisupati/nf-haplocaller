#!/usr/bin/env nextflow
/*
========================================================================================
             GATK HaplotypeCaller B E S T - P R A C T I C E
========================================================================================
 New GATK HaplotypeCaller Best Practice Analysis Pipeline. Started December 2018.
 #### Homepage / Documentation

 #### Authors
 Rahul Pisupati <rahul.pisupati@gmi.oeaw.ac.at>
----------------------------------------------------------------------------------------
*/


/*

Simply run this

nextflow run main.nf --reads "*bam" --file_ext bam --fasta ~/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta --outdir output_folder

*/

/*
 * SET UP CONFIGURATION VARIABLES
 */
params.project = "the1001genomes"
params.outdir = './snpcall'
params.fasta = false   // reference fasta file
params.file_ext = "bam"  // please change this accordingly..
params.singleEnd = false
params.tmpdir = "/lustre/scratch/users/rahul.pisupati/tempFiles/"

params.saveTrimmed = false
build_index = false
params.name = false
params.notrim = true
params.clusterOptions = false
params.email = false
params.plaintext_email = false
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
  if( ! nextflow.version.matches(">= $params.nf_required_version") ){
    throw GroovyException('Nextflow version too old')
  }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

/*
 * Create a channel for input read files
 */
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

num_files = 1
if ( params.file_ext == 'fastq' ){
  num_files = params.singleEnd ? 1 : 2
}
read_files_processing = Channel
    .fromFilePairs( params.reads, size: num_files )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }

/// ______________________________

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Library prep presets
params.illumina = true

log.info "=================================================="
log.info " nf-haplocaller : SNP calling Best Practice v${params.version}"
log.info "=================================================="
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Data Type']      = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']         = params.fasta
if(params.notrim)       summary['Trimming Step'] = 'Skipped'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "========================================="


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
      file "${refid}.fasta.*" into bwa_index
      file "${refid}.dict" into fasta_dict

      script:
      """
      samtools faidx ${genome}
      bwa index $genome
      java -jar \$EBROOTPICARD/picard.jar  CreateSequenceDictionary R=$genome O=${refid}.dict
      """
  }
} else {
  bwa_index = Channel
    .fromPath( "$reffol/${refid}.fasta.*" )
  fasta_dict = Channel
    .fromPath( "$reffol/${refid}.dict" )
}

/*
* 2. Generate FASTQ from BAM file
*/
if (params.file_ext == "fastq"){
  read_files_processing.into { read_files_fastqc; read_files_trimming }
} else {

  process extractFastq {
    tag "$name"
    storeDir "${params.tmpdir}/rawreads"
    label 'env_qual_small'

    input:
    set val(name), file(reads) from read_files_processing

    output:
    set val(name), file("${name}.fastq") into read_files_fastqc
    set val(name), file("${name}.fastq") into read_files_trimming

    script:
    if (params.singleEnd) {
      if (reads.getExtension() == "sra") {
        """
        fastq-dump $reads
        """
      } else if (reads.getExtension() == "bam") {
        """
        java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTPICARD/picard.jar SamToFastq I=$reads FASTQ=${name}.fastq VALIDATION_STRINGENCY=LENIENT
        """
      }
    } else {
      if (reads[0].getExtension() == "sra") {
        """
        fastq-dump --split-files $reads
        """
      } else if (reads.getExtension() == "bam") {
        """
        java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTPICARD/picard.jar SamToFastq I=$reads FASTQ=${name}.fastq INTER=true VALIDATION_STRINGENCY=LENIENT
        """
      }
    }
  }

}

/*
* 3. FastQC for the input files
*/
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
* 4. Trimming the reads
*/
if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = Channel.from(false)
} else {
    process trim_galore {
        tag "$name"
        label 'env_trim'
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(name), file(reads) from read_files_trimming

        output:
        set val(name), file('*fq.gz') into trimmed_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        rrbs = params.rrbs ? "--rrbs" : ''
        illumina = params.illumina ? "--illumina" : ''
        non_directional = params.rrbs && params.non_directional ? "--non_directional" : ''
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $illumina $rrbs $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $illumina $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}


/*
* 5. Aligning the reads -- BWA-MEM
*/
process alignReads {
  tag "$name"
  label 'env_bwa_large'

  input:
  set val(name), file(reads) from trimmed_reads
  file genome
  file indices from bwa_index

  output:
  set val(name), file("${name}.sam") into aligned_sam

  script:
  """
  bwa mem -t ${task.cpus} -p $reffol/${refid}.fasta $reads > ${name}.sam
  """
}

/*
* 6. Processing sam to bam and sorting the bam
*/
process processBam {
  tag "$name"
  label 'env_picard_medium'

  input:
  set val(name), file(sam) from aligned_sam

  output:
  set val(name), file("${name}.sorted.bam") into sorted_bam


  script:
  """
  samtools view -b -o ${name}.bam -S $sam
  samtools sort -m 10G --threads ${task.cpus} -o ${name}.sorted.bam ${name}.bam
  """
}

/*
* 7. Picard tool on bam file to remove duplicates
*/
process picardBam {
  tag "$name"
  label 'env_picard_small'

  input:
  set val(name), file(bam) from sorted_bam

  output:
  set val(name), file("${name}.modified.bam") into modified_bam
  set val(name), file("${name}.modified.bam.bai") into modified_bam_index

  script:
  """
  java -Djava.io.tmpdir=$params.tmpdir -jar \$EBROOTPICARD/picard.jar MarkDuplicates INPUT=$bam OUTPUT=${name}.dedup.bam METRICS_FILE=${name}.metrics
  java -Djava.io.tmpdir=$params.tmpdir -jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups INPUT=${name}.dedup.bam OUTPUT=${name}.modified.bam ID=$name LB=$name PL=illumina PU=none SM=$name
  samtools index ${name}.modified.bam
  """
}

/*
* 8. GATK to realign the reads at the positions where there are indels
*/
process realignBam {
  tag "$name"
  publishDir "${params.outdir}/alignedBam", mode: 'copy'
  label 'env_gatk_small'

  input:
  set val(name), file(bam) from modified_bam
  set val(name), file(bam_index) from modified_bam_index

  output:
  set val(name), file("${name}.realignedBam.bam") into realigned_bam
  set val(name), file("${name}.realignedBam.bam.bai") into realigned_bam_index


  script:
  """
  java -Djava.io.tmpdir=${params.tmpdir}  -jar \$EBROOTGATK/GenomeAnalysisTK.jar\
  -T RealignerTargetCreator -R $reffol/${refid}.fasta\
  -I $bam -o ${name}.intervals

  java -Djava.io.tmpdir=${params.tmpdir}  -jar \$EBROOTGATK/GenomeAnalysisTK.jar\
  -T IndelRealigner -R $reffol/${refid}.fasta\
  -I $bam -targetIntervals ${name}.intervals\
  -o ${name}.realignedBam.bam

  samtools index ${name}.realignedBam.bam
  """
}

/*
* 9. GATK HaplotypeCaller for the SNPs
*/
process doSNPcall {
  tag "$name"
  publishDir "${params.outdir}/gvcf", mode: 'copy'
  label 'env_gatk_medium'

  input:
  set val(name), file(bam) from realigned_bam
  set val(name), file(bam_index) from realigned_bam_index

  output:
  file("${name}.g.vcf.gz") into raw_gvcf
  file("${name}.g.vcf.gz.tbi") into raw_gvcf_index

  script:
  """
  java -Djava.io.tmpdir=$params.tmpdir -jar \$EBROOTGATK/GenomeAnalysisTK.jar \
  -T HaplotypeCaller -R $reffol/${refid}.fasta \
  -I $bam -o ${name}.g.vcf.gz \
  -nct ${task.cpus} \
  --genotyping_mode DISCOVERY -ERC GVCF \
  -variant_index_type LINEAR -variant_index_parameter 128000
  """
}


/*
* 10. GenotypeGVCF for all the files
*/
input_genotypegvcf = raw_gvcf.collect()
input_genotypegvcf_index = raw_gvcf_index.collect()

process joinGVCFs {
  publishDir "$params.outdir", mode: 'copy'
  label 'env_gatk_large'

  input:
  file in_vcf from input_genotypegvcf
  file(vcf_idx) from input_genotypegvcf_index

  output:
  file("allSample_Combined.vcf.gz") into combgVCF
  file("allSample_Combined.vcf.gz") into combgVCF_name
  file("allSample_Combined.vcf.gz.tbi") into combgVCF_index
  file("allSample_Combined.vcf.gz.tbi") into combgVCF_name_index

  script:
  def try_vcfs = in_vcf.collect { "-V $it" }.join(' ')
  """
  java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTGATK/GenomeAnalysisTK.jar\
  -T GenotypeGVCFs -R $reffol/${refid}.fasta\
  -nt ${task.cpus} \
  $try_vcfs -o allSample_Combined.vcf.gz \
  -variant_index_type LINEAR -variant_index_parameter 128000
  """
}

/*
* 11. Get sample names and filter the GVCF by SelectVariants
*/
process getSamples {
  tag "joinGVCF"
  label 'env_gatk_small'

  input:
  file gvcf from combgVCF_name
  file gvcf_index from combgVCF_name_index

  output:
  stdout sample_names

  script:
  """
  tabix -H $gvcf | grep -m 1 "^#CHROM" | awk '{for(i=10;i<=NF;i++) print \$i}'
  """
}

input_names = sample_names
    .splitCsv()
    .map {row -> "${row[0]}"}


process selectSNPs {
  tag "$name"
  publishDir "${params.outdir}/vcf", mode: 'copy'
  label 'env_gatk_small'

  input:
  val name from input_names
  file gvcf from combgVCF
  file gvcf_index from combgVCF_index

  output:
  file "${name}.filter.vcf" into filter_vcf

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
