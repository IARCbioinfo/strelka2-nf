#! /usr/bin/env nextflow

/*vim: syntax=groovy -*- mode: groovy;-*- */

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

params.help          		 = null
params.ref           		 = null
params.tn_pairs      		 = null
params.input_folder  		 = "./"
params.strelka        		= "$CONDA_PREFIX/share/strelka-2.9.10-0/"
params.config         		= null
params.cpu            		= "2"
params.mem           		 = "20"
params.output_folder  		= "strelka_output"
params.mode           		= "somatic"
params.AF             		= null
params.exome          		= null
params.rna            		= null
params.outputCallableRegions    = null
params.callRegions    		= "NO_FILE"
params.suffix                   = ".PASS"

log.info ""
log.info "----------------------------------------------------------------"
log.info "  Strelka2 1.1 : variant calling with Strelka2 using nextflow "
log.info "----------------------------------------------------------------"
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
    log.info "nextflow run iarcbioinfo/strelka2.nf --mode somatic --ref hg38.fa --tn_pairs pairs.txt --input_folder path/to/bam/ --strelka path/to/strelka/"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--ref                  FILE                 Genome reference file"
    log.info "--input_folder         FOLDER               Folder containing BAM files"
    log.info "--tn_pairs             FILE                 Tab delimited text file with two columns called normal and tumor"
    log.info ""
    log.info "------------------GERMLINE -----------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/strelka2.nf --mode germinal --ref hg38.fa --input_folder path/to/bam/ --strelka path/to/strelka/"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--ref                  FILE                 Genome reference file"
    log.info "--input_folder         FOLDER               Folder containing BAM files"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Optional arguments:"
    log.info "--cpu                  INTEGER              Number of cpu to use (default=2)"
    log.info "--mem                  INTEGER              Memory (in GB)"
    log.info "--output_folder        PATH                 Output directory for vcf files (default=strelka_ouptut)"
    log.info "--strelka              PATH                 Strelka installation dir"
    log.info "--config               FILE                 Use custom configuration file"
    log.info "--AF                                        get Allele frequency and Filter"
    log.info "--exome                                     automatically set up parameters for exome data"
    log.info "--rna                                       automatically set up parameters for rna data"
    log.info "--callRegions          PATH                 Region bed file"
    log.info "--outputCallableRegions                     Create a BED track containing regions which are determined to be callable"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    exit 0
} else {


workflow=""
if( params.mode=="somatic"||params.mode=="genotyping" ) { workflow= params.strelka + '/bin/configureStrelkaSomaticWorkflow.py' }
else if(params.mode=="germline"){ workflow= params.strelka + '/bin/configureStrelkaGermlineWorkflow.py' }
else { println "ERROR: wrong value for --mode option. Must be somatic, germline, or genotyping"; System.exit(0) }

if (params.config==null){ config = workflow + ".ini" } else {config=params.config}
config = file(config)

fasta_ref = file(params.ref)
fasta_ref_fai = file( params.ref+'.fai' )

exome="" ; if (params.exome) { exome="--exome" }
rna=""; if (params.rna) { rna="--rna" }
bed = file( params.callRegions ) 
tbi = file( params.callRegions+'.tbi' )
outputCallableRegions="" ; if (params.outputCallableRegions) { outputCallableRegions="--outputCallableRegions" }

/* Software information */
log.info ""
log.info "ref           	= ${params.ref}"
log.info "tn_pairs      	= ${params.tn_pairs}"
log.info "input_folder  	= ${params.input_folder}"
log.info "strelka       	= ${params.strelka}"
log.info "config        	= ${config}"
log.info "cpu           	= ${params.cpu}"
log.info "mem           	= ${params.mem}Gb"
log.info "output_folder 	= ${params.output_folder}"
log.info "mode          	= ${params.mode}"
log.info "Allele Freq.  	= ${params.AF}"
log.info "exome         	= ${params.exome}"
log.info "rna           	= ${params.rna}"
log.info "config     	  	= ${config}"
log.info "callRegions   	= ${params.callRegions}"
log.info "outputCallableRegions = ${outputCallableRegions}"
log.info ""
}

if(params.mode=="genotyping"){ 
	if( params.rna ) { workflow= params.strelka + '/bin/configureStrelkaGermlineWorkflow.py' }
	pairs = Channel.fromPath(params.tn_pairs).splitCsv(header: true, sep: '\t', strip: true)
			.map{ row -> [ row.sample , file(params.input_folder + "/" + row.tumor), file(params.input_folder + "/" + row.tumor+'.bai'), file(params.input_folder + "/" + row.normal), file(params.input_folder + "/" + row.normal+'.bai'), 
					file(params.input_folder + "/" + row.tumor+"_vs_" + row.normal + ".somatic.snvs${params.suffix}.vcf.gz") , file(params.input_folder + "/" + row.tumor+"_vs_" + row.normal + ".somatic.snvs${params.suffix}.vcf.gz.tbi") ,
					file(params.input_folder + "/" + row.tumor+"_vs_" + row.normal + ".somatic.indels${params.suffix}.vcf.gz") , file(params.input_folder + "/" + row.tumor+"_vs_" + row.normal + ".somatic.indels${params.suffix}.vcf.gz.tbi") ]}
			
	pairs2genotype = pairs.groupTuple(by: 0)
			      .map { row -> tuple(row[0] , row[1], row[2] , row[3][0] , row[4][0] ,  row[5],row[6],row[7],row[8]  ) }

	
 process run_strelkaSomaticGenotyping {
     cpus params.cpu
     memory params.mem+'GB'
     tag { sample }

     publishDir params.output_folder, mode: 'copy'

     input:
     set val(sample), file(bamT), file(baiT), file(bamN), file(baiN), file(vcfSNV), file(vcfSNVtbi), file(vcfINDEL), file(vcfINDELtbi) from pairs2genotype
     file bed
     file tbi
     file fasta_ref
     file fasta_ref_fai
     file config

     output:
     file 'strelkaAnalysis_T1/results/variants/*vcf.gz'
     file 'strelkaAnalysis_T1/results/variants/*.tbi'
     file 'strelkaAnalysis_T1/results/regions/*.bed.gz' optional true
     file 'strelkaAnalysis_T1/results/regions/*.tbi' optional true
     file 'strelkaAnalysis_T2/results/variants/*vcf.gz'
     file 'strelkaAnalysis_T2/results/variants/*.tbi'
     file 'strelkaAnalysis_T2/results/regions/*.bed.gz' optional true
     file 'strelkaAnalysis_T2/results/regions/*.tbi' optional true   

     shell:
     if (params.callRegions!="NO_FILE") { callRegions="--callRegions $bed" } else { callRegions="" }
     if( params.rna ){ 
	files1="--bam ${bamT[0]}"
	files2="--bam ${bamT[1]}"
     }else{ 
	files1="--normalBam $bamN --tumorBam ${bamT[0]}" 
	files2="--normalBam $bamN --tumorBam ${bamT[1]}"
     }
     '''
     !{baseDir}/bin/prep_vcf_bed.sh
     forcedGT=''
     for v in `ls *.vcf.gz`; do forcedGT=$forcedGT' --forcedGT '$v; done
     forcedGT=$forcedGT" --callRegions regions.bed.gz"
     !{workflow} $forcedGT !{files1} --referenceFasta !{fasta_ref} --config !{config} !{rna} !{exome} --runDir strelkaAnalysis !{callRegions} !{outputCallableRegions}
     cd strelkaAnalysis
     ./runWorkflow.py -m local -j !{params.cpu} -g !{params.mem}
     cd results/variants
     mv somatic.indels.vcf.gz !{sample}_!{bamT[0]}.somaticGT.indels.vcf.gz
     mv somatic.snvs.vcf.gz !{sample}_!{bamT[0]}.somaticGT.snvs.vcf.gz
     mv somatic.indels.vcf.gz.tbi !{sample}_!{bamT[0]}.somaticGT.indels.vcf.gz.tbi
     mv somatic.snvs.vcf.gz.tbi !{sample}_!{bamT[0]}.somaticGT.snvs.vcf.gz.tbi
     cd ../regions
     mv somatic.callable.regions.bed.gz !{sample}_!{bamT[0]}.somaticGT.callable.regions.bed.gz
     mv somatic.callable.regions.bed.gz.tbi !{sample}_!{bamT[0]}.somaticGT.callable.regions.bed.gz.tbi

     cd ../../..
     mv strelkaAnalysis strelkaAnalysis_T1
     !{workflow} $forcedGT !{files2} --referenceFasta !{fasta_ref} --config !{config} !{rna} !{exome} --runDir strelkaAnalysis !{callRegions} !{outputCallableRegions}
     cd strelkaAnalysis
     ./runWorkflow.py -m local -j !{params.cpu} -g !{params.mem}
     cd results/variants
     mv somatic.indels.vcf.gz !{sample}_!{bamT[1]}.somaticGT.indels.vcf.gz
     mv somatic.snvs.vcf.gz !{sample}_!{bamT[1]}.somaticGT.snvs.vcf.gz
     mv somatic.indels.vcf.gz.tbi !{sample}_!{bamT[1]}.somaticGT.indels.vcf.gz.tbi
     mv somatic.snvs.vcf.gz.tbi !{sample}_!{bamT[1]}.somaticGT.snvs.vcf.gz.tbi
     cd ../regions
     mv somatic.callable.regions.bed.gz !{sample}_!{bamT[1]}.somaticGT.callable.regions.bed.gz
     mv somatic.callable.regions.bed.gz.tbi !{sample}_!{bamT[1]}.somaticGT.callable.regions.bed.gz.tbi
     cd ../../..
     mv strelkaAnalysis strelkaAnalysis_T2
     '''
  }

}else{

if (params.mode=="somatic"){

  pairs = Channel.fromPath(params.tn_pairs).splitCsv(header: true, sep: '\t', strip: true)
  .map{ row -> [ file(params.input_folder + "/" + row.tumor), file(params.input_folder + "/" + row.tumor+'.bai'), file(params.input_folder + "/" + row.normal), file(params.input_folder + "/" + row.normal+'.bai') ] }

  process run_strelkaSomatic {

     cpus params.cpu
     memory params.mem+'GB' 
      
     publishDir params.output_folder, mode: 'copy'

     input:
     file pair from pairs
     file bed
     file tbi
     file fasta_ref
     file fasta_ref_fai
     file config

     output:
     file 'strelkaAnalysis/results/variants/*vcf.gz' into vcffiles
     file 'strelkaAnalysis/results/variants/*.tbi' into tbifiles
     file 'strelkaAnalysis/results/regions/*.bed.gz' optional true into regionfiles
     file 'strelkaAnalysis/results/regions/*.tbi' optional true into regiontbifiles

     shell:
     if (params.callRegions!="NO_FILE") { callRegions="--callRegions $bed" } else { callRegions="" }
     '''
     !{workflow} --tumorBam !{pair[0]} --normalBam !{pair[2]} --referenceFasta !{fasta_ref} --config !{config} !{rna} !{exome} --runDir strelkaAnalysis !{callRegions} !{outputCallableRegions}
     cd strelkaAnalysis
     ./runWorkflow.py -m local -j !{params.cpu} -g !{params.mem}
     cd results/variants
     mv somatic.indels.vcf.gz !{pair[0]}_vs_!{pair[2]}.somatic.indels.vcf.gz
     mv somatic.snvs.vcf.gz !{pair[0]}_vs_!{pair[2]}.somatic.snvs.vcf.gz
     mv somatic.indels.vcf.gz.tbi !{pair[0]}_vs_!{pair[2]}.somatic.indels.vcf.gz.tbi
     mv somatic.snvs.vcf.gz.tbi !{pair[0]}_vs_!{pair[2]}.somatic.snvs.vcf.gz.tbi
     cd ../regions
     mv somatic.callable.regions.bed.gz !{pair[0]}_vs_!{pair[2]}.somatic.callable.regions.bed.gz
     mv somatic.callable.regions.bed.gz.tbi !{pair[0]}_vs_!{pair[2]}.somatic.callable.regions.bed.gz.tbi
     '''
  }
    
  if (params.AF){
      process getAllelicFractionSomatic{

         publishDir params.output_folder, mode: 'copy'

         input:
         file vcf from vcffiles

         output:
         file '*.vcf' into passfiles

         shell:
         '''
         getAllelicFraction
         '''
      }
   }
  
}


if (params.mode=="germline"){
  bamFiles = Channel.fromFilePairs( params.input_folder + '/*.{bam,bam.bai}')

  process run_strelkaGermline {

    cpus params.cpu
    memory params.mem+'GB'
    tag { sample_Id }
      
    publishDir params.output_folder, mode: 'copy'

    input:
    set sample_Id, file(bam) from bamFiles
    file config
    file bed
    file tbi 
    file fasta_ref
    file fasta_ref_fai
    
    output:
    file 'strelkaAnalysis/results/variants/*.germline.vcf.gz*' into vcffiles

    shell:
    if (params.callRegions!="NO_FILE") { callRegions="--callRegions $bed" } else { callRegions="" }
    '''
    runDir="results/variants/"
    !{workflow} --bam !{sample_Id}.bam --referenceFasta !{fasta_ref} --config !{config} !{rna} !{exome} --runDir strelkaAnalysis !{callRegions}
    cd strelkaAnalysis
    ./runWorkflow.py -m local -j !{params.cpu} -g !{params.mem}
    mv $runDir/variants.vcf.gz $runDir/!{sample_Id}.germline.vcf.gz
    mv $runDir/variants.vcf.gz.tbi $runDir/!{sample_Id}.germline.vcf.gz.tbi
    '''
  }

  vcffiles.into{ vcffiles1 ; vcffiles2 }

  if (params.AF){
      process getAllelicFractionGermline{
         publishDir params.output_folder, mode: 'copy'

         input:
         file vcf from vcffiles1

         output:
         file '*.vcf' into passfiles

         shell:
         '''
         getAllelicFraction 
         '''
      }
  }
}

}

process filter_pass{
   cpus 1
   memory '1GB'

   publishDir params.output_folder+"/filtered/", mode: 'copy'

   input:
   file vcf from vcffiles2

   output:
   file '*_PASS.vcf.gz' into filtered

   shell:
   file_tag = vcf[0].name.replace(".vcf.gz","")
   '''
   bcftools view -f PASS -O z !{vcf[0]} -o !{file_tag}_PASS.vcf.gz
   '''
}
