#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help                     = null
params.ref                      = null
params.input_file               = null
params.input_folder             = null
params.strelka                  = "/opt/conda/envs/strelka2-nf/share/strelka-2.9.10-0/"
params.config                   = null
params.cpu                      = "2"
params.mem                      = "20"
params.output_folder            = "strelka_output"
params.mode                     = "somatic"
params.exome                    = null
params.rna                      = null
params.outputCallableRegions    = null
params.callRegions              = "NO_FILE"
params.AF                       = null

log.info ""
log.info "--------------------------------------------------------------------"
log.info "  Strelka2 1.3a : variant calling with Strelka2 using nextflow DSL2 "
log.info "--------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""


if (params.help) {

    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "-------------------SOMATIC -----------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/strelka2.nf -profile singularity --mode somatic --ref hg38.fa --input_file pairs.txt"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--ref                  FILE                 Genome reference file"
    log.info "--input_folder         FOLDER               Folder containing BAM or CRAM files"
    log.info "--input_file           FILE                 Tab delimited text file with at least two columns called normal and tumor;"
    log.info "                                            optionally a sample column and a vcf column to use mutect's --forcedGT"
    log.info ""
    log.info "------------------GERMLINE -----------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/strelka2.nf -profile singularity --mode germline --ref hg38.fa --input_folder path/to/bam/"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--ref                  FILE                 Genome reference file"
    log.info "--input_folder         FOLDER               Folder containing BAM files"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Optional arguments:"
    log.info "--mode                 STRING               Mode for variant calling; one of somatic, germline, genotyping (default: somatic)"
    log.info "--cpu                  INTEGER              Number of cpu to use (default: 2)"
    log.info "--mem                  INTEGER              Memory in GB (default: 20)"
    log.info "--output_folder        PATH                 Output directory for vcf files (default: strelka_ouptut)"
    log.info "--strelka              PATH                 Strelka installation dir (default: path inside docker and singularity containers)"
    log.info "--config               FILE                 Path to custom strelka configuration file (default: none)"
    log.info "--callRegions          PATH                 Region bed file (default: none)"
    log.info ""
    log.info "Flags:"
    log.info "--exome                                     automatically set up parameters for exome data"
    log.info "--rna                                       automatically set up parameters for rna data"
    log.info "--outputCallableRegions                     Create a BED track containing regions which are determined to be callable"
    log.info "--AF                                        Get allelic fractions"
    log.info "--help                                      Display this message"
    log.info ""
    System.exit(0) 
}


/***************************************************************************************/
/************************ handle global parameters *************************************/
/***************************************************************************************/

if (!(params.mode in ["somatic", "genotyping", "germline", "rna"])) {
        error "Invalid value for params.mode: ${params.mode}. Supported values are somatic, genotyping, germline, rna."
}

outputCallableRegionsFlag = params.outputCallableRegions ? "--outputCallableRegions" : ""
exomeFlag = params.exome ? "--exome" : ""
rnaFlag = params.rna ? "--rna" : ""

log.info ""
log.info "ref                   = ${params.ref}"
log.info "input_file            = ${params.input_file}"
log.info "input_folder          = ${params.input_folder}"
log.info "strelka               = ${params.strelka}"
log.info "cpu                   = ${params.cpu}"
log.info "mem                   = ${params.mem}Gb"
log.info "output_folder         = ${params.output_folder}"
log.info "mode                  = ${params.mode}"
log.info "exome                 = ${params.exome}"
log.info "rna                   = ${params.rna}"
log.info "callRegions           = ${params.callRegions}"
log.info "outputCallableRegions = ${params.outputCallableRegions}"
log.info "ext                   = ${params.ext}"
log.info ""


/***************************************************************************************/
/************************  Process : run_strelka_somatic *******************************/
/***************************************************************************************/

process run_strelka_somatic {

    cpus params.cpu
    memory params.mem +'GB'

    publishDir params.output_folder+"/VCFs/raw", mode: 'copy', pattern: "*somatic*vcf.gz*"
    publishDir params.output_folder+"/CallableRegions", mode: 'copy', pattern: "*callable.regions.bed.gz*"

    input:
        tuple val(sample), path(bamT), path(baiT), path(bamN), path(baiN), path(vcfSNV), path(vcfSNVtbi), path(vcfINDEL), file(vcfINDELtbi)
        tuple path(fasta_ref), path(fasta_ref_fai)
        tuple path(bed), path(tbi)

    output:
        tuple val(sample),path('*.somatic.snvs.vcf.gz'), 
                          path('*.somatic.snvs.vcf.gz.tbi'),
                          path('*.somatic.indels.vcf.gz'), 
                          path('*.somatic.indels.vcf.gz.tbi'), emit: vcffiles
        tuple val(sample),path('*regions.bed.gz'), path('*regions.bed.gz.tbi'), optional: true

    script:
        output_prefix = sample ? "${sample}.somatic" : "${bamT}_vs_${bamN}.somatic"
        config = (params.config==null) ? params.strelka + "/bin/configureStrelkaSomaticWorkflow.py.ini" : params.config
        callRegions = (vcfSNV.name!="NO_SNPS") ? "--callRegions regions.bed.gz" : (bed != "NO_FILE") ? "--callRegions $bed" : ""
        """
        forcedGT=''
        if [ -e "${vcfSNV}" ]; then
            ${projectDir}/bin/prep_vcf_bed.sh
            for v in *.vcf.gz; do forcedGT=\$forcedGT' --forcedGT '\$v; done
        fi
        
        ${params.strelka}/bin/configureStrelkaSomaticWorkflow.py \$forcedGT \
            --tumorBam $bamT \
            --normalBam $bamN \
            --referenceFasta $fasta_ref \
            --config $config \
            --runDir strelkaAnalysis \
            $exomeFlag $callRegions $outputCallableRegionsFlag


        cd strelkaAnalysis
        ./runWorkflow.py -m local -j ${task.cpus} -g ${params.mem}

        cd ..
        mv strelkaAnalysis/results/variants/* .
        mv somatic.indels.vcf.gz ${output_prefix}.indels.vcf.gz
        mv somatic.snvs.vcf.gz ${output_prefix}.snvs.vcf.gz
        mv somatic.indels.vcf.gz.tbi ${output_prefix}.indels.vcf.gz.tbi
        mv somatic.snvs.vcf.gz.tbi ${output_prefix}.snvs.vcf.gz.tbi
        ${projectDir}/bin/fixStrelkaOutput.sh *.vcf.gz

        if [ -d strelkaAnalysis/results/regions/ ]; then
            mv strelkaAnalysis/results/regions/* .
            mv somatic.callable.regions.bed.gz ${output_prefix}.callable.regions.bed.gz
            mv somatic.callable.regions.bed.gz.tbi ${output_prefix}.callable.regions.bed.gz.tbi
        fi
        """

    stub:
        output_prefix = sample ? "${sample}.somatic" : "${bamT}_vs_${bamN}.somatic"
        """
        touch ${output_prefix}.snvs.vcf.gz
        touch ${output_prefix}.snvs.vcf.gz.tbi
        touch ${output_prefix}.indels.vcf.gz
        touch ${output_prefix}.indels.vcf.gz.tbi
        touch ${output_prefix}.indels.vcf.gz.tbi
        touch ${output_prefix}.callable.regions.bed.gz
        """
}

/***************************************************************************************/
/************************  Process : run_strelka_germline ******************************/
/***************************************************************************************/

process run_strelka_germline {

    cpus params.cpu
    memory params.mem + 'GB'

    publishDir params.output_folder+"/VCFs/raw", mode: 'copy', pattern: "*germline*vcf.gz*"

    input:
        tuple val(sample), path(bam), path(bai), path(vcfSNV), path(vcfSNVtbi), path(vcfINDEL), file(vcfINDELtbi)
        tuple path(fasta_ref), path(fasta_ref_fai)
        tuple path(bed), path(tbi)
    
    output:
        tuple val(sample),path('*.germline.snvs.vcf.gz'), 
                          path('*.germline.snvs.vcf.gz.tbi'), 
                          path('*.germline.indels.vcf.gz'), 
                          path('*.germline.indels.vcf.gz.tbi'), emit: vcffiles

    script:
        output_prefix = sample ? "${sample}.germline" : "${bam}.germline"
        config = (params.config==null) ? params.strelka + "/bin/configureStrelkaGermlineWorkflow.py.ini" : params.config
        callRegions = (vcfSNV.name!="NO_SNPS") ? "regions.bed.gz" : (bed != "NO_FILE") ? "--callRegions $bed" : ""

        """
        forcedGT=''
        if [ -e "${vcfSNV}" ]; then
            ${projectDir}/bin/prep_vcf_bed.sh
            for v in *.vcf.gz; do forcedGT=\$forcedGT' --forcedGT '\$v; done
        fi

        ${params.strelka}/bin/configureStrelkaGermlineWorkflow.py \$forcedGT \
            --bam ${bam} \
            --referenceFasta ${fasta_ref} \
            --config ${config} \
            --runDir strelkaAnalysis \
            $rnaFlag $exomeFlag $callRegions

        cd strelkaAnalysis
        ./runWorkflow.py -m local -j ${task.cpus} -g ${params.mem}
        cd ..
        mv strelkaAnalysis/results/variants/* .
        mv variants.vcf.gz ${output_prefix}.vcf.gz
        mv variants.vcf.gz.tbi ${output_prefix}.vcf.gz.tbi
        """

    stub:
        output_prefix = sample ? "${sample}.germline" : "${vcfSNV}.germline"
        """
        touch ${output_prefix}.snvs.vcf.gz
        touch ${output_prefix}.snvs.vcf.gz.tbi
        touch ${output_prefix}.indels.vcf.gz
        touch ${output_prefix}.indels.vcf.gz.tbi
        """
    
}

/***************************************************************************************/
/************************  Process : run_strelkaGenotyping *****************************/
/***************************************************************************************/

process run_strelkaGenotyping {

    cpus params.cpu
    memory params.mem + 'GB'
    
    publishDir params.output_folder+"/VCFs/GT/raw", mode: 'copy', pattern: "*vcf*"
    publishDir params.output_folder+"/CallableRegions", mode: 'copy', pattern: "*bed*"

    input:
        tuple val(sample), path(bamT), path(baiT), path(bamN), path(baiN), path(vcfSNV), path(vcfSNVtbi), path(vcfINDEL), path(vcfINDELtbi)
        tuple path(fasta_ref), path(fasta_ref_fai)
        tuple path(bed), path(tbi)

    output:
        tuple val(sample), path("*"+${bamT[0]}+"*.somaticGT.snvs.vcf.gz"),
                           path("*"+${bamT[0]}+"*.somaticGT.snvs.vcf.gz.tbi"),
                           path("*"+${bamT[0]}+"*.somaticGT.indels.vcf.gz"),
                           path("*"+${bamT[0]}+"*.somaticGT.indels.vcf.gz.tbi"),
                           path("*"+${bamT[1]}+"*.somaticGT.snvs.vcf.gz"),
                           path("*"+${bamT[1]}+"*.somaticGT.snvs.vcf.gz.tbi"),
                           path("*"+${bamT[1]}+"*.somaticGT.indels.vcf.gz"),
                           path("*"+${bamT[1]}+"*.somaticGT.indels.vcf.gz.tbi"), emit: vcffiles
        tuple val(sample),path('*regions.bed.gz'), path('*regions.bed.gz.tbi'), optional: true
     
    shell:
        config = (params.config==null) ? params.strelka + "/bin/configureStrelkaSomaticWorkflow.py.ini" : params.config
        callRegions = "regions.bed.gz"

        if( params.rna ){ 
	        files1="--bam ${bamT[0]}"
	        files2="--bam ${bamT[1]}"
        }else{ 
            files1="--normalBam $bamN --tumorBam ${bamT[0]}" 
            files2="--normalBam $bamN --tumorBam ${bamT[1]}"
        }
        
        """
        forcedGT=''
        ${projectDir}/bin/prep_vcf_bed.sh
        for v in *.vcf.gz; do forcedGT=\$forcedGT' --forcedGT '\$v; done

        ${params.strelka}/bin/configureStrelkaSomaticWorkflow.py \$forcedGT /
            --referenceFasta $fasta_ref /
            --config $config /
            --runDir strelkaAnalysis_T1 /
            $files1 $exomeFlag $callRegions $outputCallableRegionsFlag

        cd strelkaAnalysis_T1
        ./runWorkflow.py -m local -j ${task.cpus} -g ${params.mem}
        cd ..
        mv strelkaAnalysis_T1/results/variants/* .
        mv strelkaAnalysis_T1/results/regions/* .
        mv somatic.indels.vcf.gz ${sample}_${bamT[0]}.somaticGT.indels.vcf.gz
        mv somatic.snvs.vcf.gz ${sample}_${bamT[0]}.somaticGT.snvs.vcf.gz
        mv somatic.indels.vcf.gz.tbi ${sample}_${bamT[0]}.somaticGT.indels.vcf.gz.tbi
        mv somatic.snvs.vcf.gz.tbi ${sample}_${bamT[0]}.somaticGT.snvs.vcf.gz.tbi
        mv somatic.callable.regions.bed.gz ${sample}_${bamT[0]}.somaticGT.callable.regions.bed.gz
        mv somatic.callable.regions.bed.gz.tbi ${sample}_${bamT[0]}.somaticGT.callable.regions.bed.gz.tbi

        ${params.strelka}/bin/configureStrelkaSomaticWorkflow.py \$forcedGT /
            --referenceFasta $fasta_ref /
            --config $config /
            --runDir strelkaAnalysis_T2 /
            $files2 $exomeFlag $callRegions $outputCallableRegionsFlag

        cd strelkaAnalysis_T2
        ./runWorkflow.py -m local -j ${task.cpus} -g ${params.mem}
        cd ..
        mv strelkaAnalysis_T1/results/variants/* .
        mv strelkaAnalysis_T1/results/regions/* .
        mv somatic.indels.vcf.gz ${sample}_${bamT[1]}.somaticGT.indels.vcf.gz
        mv somatic.snvs.vcf.gz ${sample}_${bamT[1]}.somaticGT.snvs.vcf.gz
        mv somatic.indels.vcf.gz.tbi ${sample}_${bamT[1]}.somaticGT.indels.vcf.gz.tbi
        mv somatic.snvs.vcf.gz.tbi ${sample}_${bamT[1]}.somaticGT.snvs.vcf.gz.tbi
        mv somatic.callable.regions.bed.gz ${sample}_${bamT[1]}.somaticGT.callable.regions.bed.gz
        mv somatic.callable.regions.bed.gz.tbi ${sample}_${bamT[1]}.somaticGT.callable.regions.bed.gz.tbi
        """
        

    stub:
        """
        touch ${sample}_${bamT[0]}.somaticGT.indels.vcf.gz
        touch ${sample}_${bamT[0]}.somaticGT.snvs.vcf.gz
        touch ${sample}_${bamT[0]}.somaticGT.indels.vcf.gz.tbi
        touch ${sample}_${bamT[0]}.somaticGT.snvs.vcf.gz.tbi
        touch ${sample}_${bamT[0]}.somaticGT.callable.regions.bed.gz
        touch ${sample}_${bamT[0]}.somaticGT.callable.regions.bed.gz.tbi
        touch ${sample}_${bamT[1]}.somaticGT.indels.vcf.gz
        touch ${sample}_${bamT[1]}.somaticGT.snvs.vcf.gz
        touch ${sample}_${bamT[1]}.somaticGT.indels.vcf.gz.tbi
        touch ${sample}_${bamT[1]}.somaticGT.snvs.vcf.gz.tbi
        touch ${sample}_${bamT[1]}.somaticGT.callable.regions.bed.gz
        touch ${sample}_${bamT[1]}.somaticGT.callable.regions.bed.gz.tbi
        """

}


process filter_pass{
   
    cpus 1
    memory '1GB'

    publishDir params.output_folder+"/VCFs/filtered/", mode: 'copy', pattern: "*vcf.gz*"

    input:
        tuple val(sample), path(snvs), path(snvs_tbi), path(indels), path(indels_tbi)

    output:
        tuple val(sample),path('*snvs.pass.vcf.gz'), 
                          path('*.snvs.pass.vcf.gz.tbi'), 
                          path('*.indels.pass.vcf.gz'), 
                          path('*.indels.pass.vcf.gz.tbi'), emit: vcffiltred

    script:
        snvs_tag = snvs[0].name.replace(".vcf.gz",".pass.vcf.gz")
        indels_tag = indels[0].name.replace(".vcf.gz",".pass.vcf.gz")
        """
        bcftools view -f PASS -O z ${snvs} -o ${snvs_tag}
        bcftools index -t ${snvs_tag}
        bcftools view -f PASS -O z ${indels} -o ${indels_tag}
        bcftools index -t ${indels_tag}
        """

    stub:
        snvs_tag = snvs[0].name.replace(".vcf.gz",".pass.vcf.gz")
        indels_tag = indels[0].name.replace(".vcf.gz",".pass.vcf.gz")
        """
        touch ${snvs_tag}
        touch ${snvs_tag}.tbi
        touch ${indels_tag}
        touch ${indels_tag}.tbi
        """
    
}

/****************************************************************************************/
/************************  Workflow   ***************************************************/
/****************************************************************************************/

workflow {

    if(params.input_file && params.mode=="somatic"){
        println("input from file")
        input = Channel.fromPath(params.input_file).splitCsv(header: true, sep: '\t', strip: true) | map { 
            row -> 
            assert (row.tumor != null ) : "Error: tumor file is missing check your input_file"
            assert (row.normal != null) : "Error: normal file is missing check your input_file"
            tuple(
                row.sample, 
                file(row.tumor),file("${row.tumor}.{bai,crai}"),
                file(row.normal),file("${row.normal}.{bai,crai}"),
                row.snp ? file(row.snp) : file("NO_SNPS"), 
                row.snp ? file(row.snp + ".tbi") : file("NO_SNPS.TBI"),
                row.indels ? file(row.indels) : file("NO_INDEL"), 
                row.indels ? file(row.indels + ".tbi") : file("NO_INDEL.TBI")
            )}
    } else if(params.input_file && params.mode=="germline"){
        println("input from file")
        input = Channel.fromPath(params.input_file).splitCsv(header: true, sep: '\t', strip: true) | map { 
            row -> 
            assert (row.tumor != null ) : "Error: tumor file is missing check your input_file"
            tuple(
                row.sample, 
                file(row.tumor),file("${row.tumor}.{bai,crai}"),
                row.snp ? file(row.snp) : file("NO_SNPS"), 
                row.snp ? file(row.snp + ".tbi") : file("NO_SNPS.TBI"),
                row.indels ? file(row.indels) : file("NO_INDEL"), 
                row.indels ? file(row.indels + ".tbi") : file("NO_INDEL.TBI")
            )}
    } else if (params.input_folder && params.mode=="germline"){
        println("input from folder")
        input = Channel.fromFilePairs( params.input_folder + "/*.{bam,bam.bai,cram,cram.crai}",flat: true) | map{
            row -> tuple(
                row[0],row[1],row[2],file("NO_SNPS"), file("NO_SNPS.TBI"), file("NO_INDEL"), file("NO_INDEL.TBI")
            )}
    } else {
        error "ERROR : NO INPUTS !  You must provide --input_file (with all modes) or --input_folder (with mode germline)"
    }

    fasta_ref = tuple( file(params.ref), file( params.ref+'.fai' ) )

    bed = tuple( file( params.callRegions ), file( params.callRegions+'.tbi' ) )
    
    /************ process *****************/
    
    if(params.mode=="somatic"){
        println "Entering somatic mode"
        println(input)
        run_strelka_somatic(input,fasta_ref,bed)
        filter_pass(run_strelka_somatic.out.vcffiles)
    } else if(params.mode=="germline"){
        println "Entering germline mode"
        run_strelka_germline(input,fasta_ref,bed) | filter_pass
    }else if(params.mode=="genotyping"){
        println "Entering genotyping mode"
        run_strelka_genotyping(input,fasta_ref,bed)
    }else{
        error "ERROR : WRONG MODE !  --mode must be one of somatic, germline, genotyping"
    }
    
    
    
}

