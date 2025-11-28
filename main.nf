nextflow.enable.dsl=2

//parameters

//input reads
params.query = "data/*.fastq"

//reference reads
params.reference = "target/*.fna"

//output directory
params.outdir = "results"

//Processes 
process fastQC {

 publishDir "${params.outdir}/fastqc", mode: "copy"

 input:
 path query 
 
 output:
 path "fastqc_out"

 script:
 """
 mkdir -p fastqc_out
 fastqc -o fastqc_out ${query}

 """
}

process fastplong {

 publishDir "${params.outdir}/fastplong", mode: "copy"

 input:
 path query
 
 output:
 path "fastp_long.fq"

 script:
 """
 fastplong --disable_adapter_trimming -i ${query} -o fastp_long.fq
 """
}

process hifiasm {

 publishDir "${params.outdir}/hifiasm", mode: "copy"
 
 input:
 path  fastp_output

 output:
 path "hifiasm_out"
 
 script:
 """
 mkdir -p hifiasm_out
 hifiasm -t ${task.cpus} -o hifiasm_out/assembly ${fastp_output} 2>&1 | tee hifiasm_out/hifiasm.log
 """

}

process gfa_tools {

 publishDir "${params.outdir}/hifiasm_conversion", mode: "copy"

 input:
 path hifiasm_out

 output:
 path "hifiasm_conversion"

 script:
 """
 mkdir -p hifiasm_conversion

 # find gfa primary file and take the first one.
 gfa_file=\$(find ${hifiasm_out}/ -name "*.p_ctg.gfa" | head -1)

 # if gfa primary file found then 
 if [ -n "\$gfa_file" ]; then
    gfatools gfa2fa "\$gfa_file" > hifiasm_conversion/hifiasm_assembly.fa
 # else look for any .gfa file
 else
     gfa_file=\$(find ${hifiasm_out}/ -name "*.gfa" | head -1)
     if [ -n "\$gfa_file" ]; then
         gfatools gfa2fa "\$gfa_file" > hifiasm_conversion/hifiasm_assembly.fa
	# if none found then exit 1
     else
         echo "No GFA files found in ${hifiasm_out}"
         exit 1
     fi
 fi
"""
}

process minimap2 {

 publishDir "${params.outdir}/minimap2_output", mode: "copy"

 input:
 path hifiasm_conversion
 path query

 output:
 path "minimap2_output"

 script:
 """
 mkdir -p minimap2_output
 
 #create reference index 
 minimap2 -t ${task.cpus} -d hifiasm_assembly.mmi ${hifiasm_conversion}/hifiasm_assembly.fa

 #use index to align and pipe the output straight to bam file
 minimap2 -t${task.cpus} -a -x map-hifi hifiasm_assembly.mmi ${query} | \\
     samtools view -b - > minimap2_output/alignment.bam
 
 #generate stats for MultiQC
 samtools stats minimap2_output/alignment.bam > minimap2_output/alignment_stats.txt
 samtools flagstat minimap2_output/alignment.bam > minimap2_output/alignment_flagstat.txt
 """
}

process merqury {

 publishDir "${params.outdir}/merqury_output", mode: "copy"

 input:
 path hifiasm_conversion


 output:
 path "merqury_out_prefix*"

 script:
 """
 export MERQURY="/usr/local/share/merqury"
 #count k-mers in reads using meryl and save as a .meryl file
 meryl count k=21 ${hifiasm_conversion}/hifiasm_assembly.fa output merqury_out.meryl

 #run merqury
 merqury.sh merqury_out.meryl ${hifiasm_conversion}/hifiasm_assembly.fa merqury_out_prefix

 """
}

process busco {

 publishDir "${params.outdir}/busco_output", mode: "copy"

 input:
 path hifiasm_conversion

 output:
 path "busco_output"

 script:
 """
 mkdir -p busco_output
 busco -i ${hifiasm_conversion}/hifiasm_assembly.fa -m genome -l saccharomycetes_odb10 -o busco_out --out_path busco_output

 """
}

process quast {

 publishDir "${params.outdir}/quast_output", mode: "copy"

 input:
 path hifiasm_conversion
 path reference

 output:
 path "quast_output"

 script:
 """
 quast.py ${hifiasm_conversion}/hifiasm_assembly.fa -r ${reference} -o quast_output
  
 """
}

process multiqc {

 publishDir "${params.outdir}/multiqc_output", mode: "copy"

 input:
 path multiqc_files

 output:
 path "multiqc_report.html", emit: report
 path "multiqc_report_data", emit: data

 script:
 """
 multiqc . --filename multiqc_report.html

 """
}

workflow {
  query_ch = Channel.fromPath(params.query)
  reference_ch = Channel.fromPath(params.reference)
   
  fastQC(query_ch)

  fastplong(query_ch) | hifiasm | gfa_tools 

  minimap2(gfa_tools.out, query_ch )
  merqury(gfa_tools.out) 
  busco(gfa_tools.out)
  quast(gfa_tools.out, reference_ch)

  //collect outputs for multiqc - not including merqury as it is not supported
    multiqc_files = fastQC.out.mix(
    hifiasm.out,
    minimap2.out,
    busco.out,
    quast.out
)
  multiqc(multiqc_files.collect())
}



