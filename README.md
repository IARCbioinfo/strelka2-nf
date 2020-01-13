# strelka2-nf
### Strelka v2 pipeline with Nextflow

#### Dependencies
1. Install [Strelka v2](https://github.com/Illumina/strelka).
2. Install [nextflow](http://www.nextflow.io/).

	```bash
	curl -fsSL get.nextflow.io | bash
	```
	And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
	```bash
	sudo mv nextflow /usr/local/bin
	```

#### Execution

mode somatic
`nextflow run iarcbioinfo/strelka2.nf --mode somatic --ref hg38.fa --tn_pairs pairs.txt --input_folder path/to/bam/ --strelka path/to/strelka/`

mode germline
`nextflow run iarcbioinfo/strelka2.nf --mode germline --ref hg38.fa --input_folder path/to/bam/ --strelka path/to/strelka/`

#### Options
--rna
--exome
--callRegions
--outputCallableRegions

#### Help section
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/strelka2.nf --help
```
This shows details about optional and mandatory parameters provided by the user.  

#### pairs.txt format
The pairs.txt file is where you can define pairs of bam to analyse with strelka. It's a tabular file with 2 columns normal and tumor.

| normal | tumor |
| ----------- | ---------- |
| normal1.bam | tumor2.bam |
| normal2.bam | tumor2.bam |
| normal3.bam | tumor3.bam |

#### Global parameters
```--strelka``` and ```--ref``` are mandatory parameters but can be defined in your nextflow config file (```~/.nextflow/config``` or ```config``` in the working directory) and so not set as inputs.

The following is an example of config part defining this:
```bash
profiles {

        standard {
                params {
                   ref = '~/Documents/Data/references/hg38.fasta'
                   strelka = '~/bin/strelka/1.0.15/bin/'
                }
        }
```

variants.vcf.gz

## Output
  | Type      | Description     |
  |-----------|---------------|
  | strelkaAnalysis/results/variants/\*.vcf.gz    | VCF files |
  | filtered/\*PASS.vcf.gz    | VCF files with only variants with PASS flag |
  
  All vcf files have companion tabix index files (.tbi). Note that in germline mode, the VCF outputted corresponds to variants only (file variants.vcf.gz from strelka). 
