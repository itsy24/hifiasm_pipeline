## Nextflow pipeline incorporating hifiasm

### Folder Structure
1. Sequencing fastq reads must be saved in a folder named "data"
2. Reference reads must be saved in a folder named "target"

# Hifiasm pipeline that pre-processes raw reads and builds a phased assembly
By: Isabel Sy

## 1. Introduction
Next generation sequencing techniques alongside the rapid development of technology has led to a plethora of bioinformatic tools available for downstream data analysis. The influx of these tools has led to an overwhelming sea of possibilities. A subset of these tools are assemblers which as the name suggests, assembles a whole genome or DNA sequence from DNA fragments. The assembly of a high quality phased genome is critical for studying genetic variation, structural rearrangements, and haplotype-specific effects on traits and diseases. Their downstream applications in functional genomics and clinical research would enable more accurate discovery of gene variants and regulatory mechanisms. Of all the assemblers currently available, HiFiasm, a long read assembly specific for PacBio data has emerged to be a powerful haplotype de novo assembler. Hifiasm was created due to current limitations in algorithms to create consensus genomic copies<sup>1</sup>. The program itself is great if the goal is to create phased assemblies, but that's not usually the case. This is where workflows or pipelines come into play. 

Bioinformatic pipelines are a sequence of tools and/or algorithms designed to automate a process. These tools are incorporated in academia and industry and are utilized in almost all sectors. While Hifiasm has emerged as a powerful assembler, there is no standardized pipeline that incorporates Hifiasm and therefore allowing reproducible workflows. This project builds a NextFlow<sup>2</sup> pipeline that preprocesses HiFi data to assemble a phase-resolved genome assembly with Hifiasm. Utilizing Nextflow as workflow management allows for checkpointing, reproducibility, and scalability. The pipeline will also utilize other tools for read quality control checks, assembly quality check and analysis (Figure 1)<sup>3-9</sup>. The final output is reproducible nextflow pipeline that preprocesses HiFi data, assembles phased genomes using Hifiasm, validates assemblies with multiple QC tools, and generates an aggregated report. This pipeline will serve as a standardized, user-friendly framework for long-read genome assembly.

![Screenshot 2025-10-13 135132.png](attachment:88440710-d563-4932-a93e-94ac9da450cc.png)
<p><em>Figure 1. Pipeline Workflow</em></p>


## 2. Materials & Methods
    
### Data Acquisition
Validation pipeline testing was performed using Candida albicans SC5314 data from a study by Hoyer et al. This dataset was chosen because yeast organisms are typically smaller (~12 million base pairs), allowing for rapid assembly and reduced computational load during iterative pipeline testing. Raw fastq files and reference assembly were downloaded from the National Library of Medicine - National Center for Biotechnology Information (NCBI) Sequence Read Archive(SRA) under accession SRR23724250 <sup>10</sup>. 

### Computational Environment
Analysis were executed on the New York University (NYU) high-performance computing (HPC) cluster. Jobs were submitted using SLURM with max 32 CPU cores and max 64GB RAM. The workflow was implemented in Nextflow DSL2 (version 25.04.3) with tools execution managed through the NYU HPC module list or Docker/Singularity containers. The full pipeline script and parameters can be found here: https://github.com/itsy24/hifiasm_pipeline.
    
### Pipeline Workflow

The pipeline consists of Four major stages:
    1. Quality Control and Read preprocessing
    2. Genome Assembly and file conversion
    3. Structural validation
    4. Assembly quality assessment
    5. Reporting
All the tools were ran using default parameters. 

#### 1. Quality Control and Read preprocessing
Quality control (QC) of raw HiFi sequencing reads was conducted with FASTQC (v.0.11.9). FastQC provides base-quality distribution, GC content, adapter content, and other read metrics that could be used to identify sequencing anomalies. Read preprocessing for preliminary filtering and per-read quality correction was then performed with FastPLong (v.0.4.1), which conducts quality filtering, low-quality base trimming, and k-mer profiling. These two steps ensure that downstream analysis is performed on high quality reads.
    
#### 2. Genome Assembly and file conversion
Genome assembly was performed using Hifiasm (v.0.20.0), a de novo assembler for PacBio HiFi reads. The assembler default parameters were utilized and produced a phased assembly in GFA format. The GFA primary output files were then converted to FASTA format using GFAtools (v.0.5.5) to generate the contig sequences.

#### 3. Structural validation 
Structural assembly correctness was evaluated by aligning the converted FASTA file to the original HiFi reads using Minimap2 (v.2.22) with the parameter -x map-hifi preset. SAM/BAM conversion, sorting, and indexing were done with Samtools (v.1.22.1). The alignment serves as the basis for measuring mapping rates and identifying potential structural inconsistencies.

#### 4. Assembly quality assessment
Three complementary tools were used to assess overall assembly quality:
- Merqury (v.1.3): Estimated assembly quality (QV score) and completeness using k-mer comparison of the raw HiFi reads against the assembled genome.
- BUSCO (v.6.0.0): Estimated gene content completeness of assemblies using orthologs by measuring recovery of conserved single-copy orthologs. In the test    validation, the completeness was based off of the lineage: saccharomycetes_odb10.
- QUAST (v.5.3): Provided structural metrics reporting such as N50, assembly size, and misassemblies relative to the reference genome

Overall, these three tools, give per-base accuracy, gene content completeness, and structural contiguity of the final assembly.

#### 5. Reporting
Outputs from all steps except Merqury were aggregated using MultiQC (v.1.32) to generate a unified HTML report. All intermediate and final files were stored in organized output directories for further downstream analysis. 

## 3. Results & Discussions
![fastqc_sequence_counts_plot.png](attachment:c209485d-68da-42ff-a8e3-981dff4314d7.png)
<p><em>Figure 2. Bar plot of total and duplicate reads for sample: SRR23724250. Majority of the reads are unique, with a samll percentage of duplicate reads. The low duplication rate is typical for long-read Hifi sequencing and reflects high library complexity and minimal data redundancy.</em></p>

![fastqc_per_sequence_quality_scores_plot.png](attachment:e6d194c8-80d3-423e-94ce-75a05da8288f.png)
<p><em>Figure 3. Histogram of mean per-sequence Phred quality scores. Most of the reads exhibit very high mean quality score (Phred > 90). The stron right-skewed peak in the green region confirms excellent overall read accuracy with no substantial population of low quality reads.</em></p>

![fastqc_adapter_content_plot.png](attachment:b3996d22-ba3f-4321-a271-de6b12c87311.png)
<p><em>Figure 4. Adapter content graph showing the percent of reads containing adapter sequences across reads positions. Both samples had nearly zero adapter contamination, indicating that the HiFi reads are already well-trimmed and free of residual adapter sequences. The green background reflects no adapter-related quality issues detected.</em></p>

![fastqc_per_sequence_gc_content_plot.png](attachment:42a2834f-f8f4-4ccb-961c-a19e615feb9b.png)
<p><em>Figure 5. Distribution of GC content across all sequencing reads. A unimodal GC peak is centered at 34%, which matches the expected composition for Candia albicans<sup>10</sup>. The sharp, narrow peak and low frequency of extreme GC values indicate high quality sequencing with minimal GC bias.</em></p>

![hifiasm-kmer-graph.png](attachment:36038358-8abe-45c2-954b-4d119465ac5c.png)
<p><em>Figure 6. K-mer distribution generated by Hifiasm. The plots shows distinct k-mers observed in the raw HiFi reads. The two broad peaks represent the heterozygous and homozygous k-mer populations typical of diploid genomes. The steep decline at low multiplicities show remove of sequencing errors, while the smooth curve at higher multiplicities indicate uniform coverage and high read accuracy. Overall shape of the distribution shows suffcient depth and data quality to support a high confidence assembly.</em></p>


![samtools_alignment_plot.png](attachment:06409fc0-cf67-4b52-899c-eebc68b60636.png)
<p><em>Figure 7. Alignment statistics for HiFi input reads mapped to the Candida albicans assembly. Reads were first aligned with Minimap2 using the -x map-hifi preset, followed by alignment metrics computed using Samtools stats. A majority of the reads map with a positive mapping quality (MQ > 0), which indicates a strong and confident alignment to the assembled genome. A small portion of the reads had MQ = 0, which may be from repetitive or low-complexity regions. Overall, the alignment profile demonstrates excellent agreement between the raw input HiFi reads and the assembled contigs.</em></p>

![busco_plot_saccharomycetes_odb10.png](attachment:5991caad-1fb9-44a5-b3d1-9c291b139f1a.png)
<p><em>Figure 8. BUSCO completeness for the assembled Candida albicans genome using saccharomycetes_odb10 as the lineage. The bar labeled "short" corresponds to the main general BUSCO summary, while the bar labeled "short_summary.specific.saccharomycetes_odb10.busco_out" is the reference-specifc summary. This reference summary is only created when BUSCO is ran with a lineage dataset. The assembly shows a very high proprtion of complete and single-copy BUSCOs (green portion), with only a small portion of duplicated, fragmented, or missing orthologs. These results indicate excellent gene content completeness and suggests minimal assembly fragmentation or gene loss. Consistency between the two summarys relfect stable results across the pipeline outputs.</em></p>


![Screenshot 2025-12-01 150247.png](attachment:20e75ce6-4475-4e70-83b6-63ede3972ee4.png)
<p><em>Figure 9. QUAST summary of genome assembly quality metrics for Candida albicans. The assembly has an N50 of 1247.6Kbp, largest contig is 3233.6Kbp in size, covers 99.4% of the reference genome, and has a average nucleotide identity (ANI) of 99.7%. The results from this table reveal a highly complete assembly with strong structural accuracy and high sequence identity.</em></p>


![Screenshot 2025-12-01 152007.png](attachment:6af62ca0-d1c0-4c22-a254-eb5ef277be8f.png)
<p><em>Figure 10. QUAST distribution of contig sizes. 193 total contigs were created, with 45 of them >50,000bp, 133 between 25,000-50,000bp, and 15 between 10,000-25,000bp. Presence of many medium and large size contigs with an N50 >1.2Mb reflects assembly contiguity and effective reconstruction of major genomic regions.</em></p>


![merqury_out_prefix.hifiasm_assembly.spectra-cn.fl.png](attachment:4471c3d7-45c2-4e05-8fe4-401a8f22a507.png)
<p><em>Figure 11. Merqury k-mer copy number spectrum (spectra-cn) showing k-mer multiplicities in the raw reads compared to their copy number in the Hifiasm assembly. The read peak at ~300-350x is the homozygous k-mers, while the the smaller peak at ~175-180x are the heterozygous k-mers. These peaks are expected in a diploid organisms such as Candida albicans. The tall orange peak at the lower multiplicities represents sequencing error k-mers present in the reads, but not in the assembly. Minimal signal in the higher-copy bins (>4) indicates repetitive elementers were assembled without excessive collapse. Demonstrates appropriate diploid structure, low sequencing noise, and high agreement between reads and final assembly.</em></p>

![merqury_out_prefix.spectra-asm.fl.png](attachment:a7ea87f5-456d-4b54-8d3e-f4a86c4ae88c.png)
<p><em>Figure 12. Merqury assembly spectrum (sepctra-asm) comparing k-mers present in the Hifiasm assembly (gray) to the raw reads k-mers absent from the assembly (red). Assembly curve has a dominant peak at ~300-350x, which is consistent with correctly assembled homozygous regions. Small heterozygous peak at ~175-180x is the preserved diploid variation. The very low proportion of read only k-mers indicated that nearly all valid k-mers from the reads were incorporated in the assembly. This graph shows that the assembly was highly accurate and complete at the k-mer level.</em></p>

This project successfully implemented the NextFlow pipeline creating a complete genome assembly with hifiasm and evaluation using Candida Albicans PacBio HiFi data. The resulting assembly had high contiguity (N50 was 1.25 Mb), high completeness (BUSCO showed ~99% complete genes), and a strong accuracy (Merqury QV was high, values ranging from 33.7 to >60, with most large contigs QV >40-50). Together, these metrics indicate that the pipeline reliably created a high-quality Candida albicans assembly. FastQC metrics revealed that the HiFi sequencing reads had uniformly high per-base quality scores with almost all the reads having a mean Phred score of 90 (Fig. 3). The minimal adapter contamination and low duplication rates further support the high integrity of the data (Fig. 4). These quality metrics are essential for Hifiasm to perform accurately because it relies on long range information to resolve repetitive regions and phase haplotypes. Hifiasm itself performed as expected. Producing an assembly that was 22.3Mb in length, which is consistent with Candid albicans. An N50 score of 1.25Mb (Fig. 9) and several contigs >3Mb (Fig. 10) emphasize the assembler's ability to reconstruct long genomic regions without fragmentation. The Hifiasm k-mer spectrum (Fig.6) showed homozygous peaks with minimal error associated k-mers, further supporting sufficient coverage and sequence accuracy.

Assembly alignment with raw reads was then performed to determine assembly completeness and accuracy. Minimap2 alignment statistics (performed with Samtools) revealed that an overwhelming amount of the reads mapped to the assembly (MQ>0), with a small amount of MQ=0 and unmapped reads observed. BUSCO analysis with saccharomycetes_odb10 set as lineage identified ~99% complete single-copy orthologs with very few missing/fragmented BUSCOs. QUAST showed a reported genome fraction of 99.4% and an ANI of 99.7%. Although 151 misassembles were detected, these do not necessarily represent assembly errors since genomic rearrangements could occur. Merqury's k-mer spectra cn plot showed clear peaks consistent with homozygous and heterozygous k-mer distributions, while the spectra asm plot had nearly all k-mers from the raw reads in the assembly. Overall, the assembly and alignment metrics performed reveal that the Hifiasm assembly is of excellent quality with high level gene completeness and demonstrated that essential genes were fully recovered. 

Despite the successful assembly and evaluation of Candida albicans, several challenges were encountered during pipeline development. First, the configuration file had problems, due to NYU HPC requiring a specific NextFlow configuration. Second, Merqury initially failed due to a missing envinroment variable (/usr/local/bin/merqury.sh: line 16: /util/util.sh: No such file or directory), which was resolved by defining export MERQURY="/usr/local/share/merqury" within the workflow. Third, out of memory failures occurred with some of the processes such as Hifiasm, Merqury, and Minimap2/Samtools and required resource adjustment. A technical limitation of this pipeline is the absence of scaffolding, which may be why fragmentation was observed in mid-sized contigs. Additionally, the workflow's memory allocation and parameters were tuned specificly for Candida albicans and have only been validated on this relatively small fungal genome. 

Future work should incorporate Hi-C scaffolding and variant calling modules to produce chromosome level assembly. Integrating containerization and expanding the script for modularization would improve reproducibility and broaden applicability across diverse organisms and sequencing datasets.

Overall, this workflow provides a robust, reproducible, and efficient framework for de novo assembly using long-read PacBio Hifi sequencing. This pipeline has broad applications both in academia and industry, enabling researchers to address key biological and genomic questions. These questions could include assessing assembly completeness, identifying species-specific genomic features, finding evidence of genome duplication or heterozygosity, and reconstruction of high-quality genomes from metagenomic or environmental samples. Furthermore, high confidence assemblies produced through this workflow can directly improve downstream analysis such as variant detection, resistance profiling, and comparative genomics.

## 4. Conclusion

This project was able to successfully develope and implement a reproducible Nextflow workflow for de nove genome assembly and quality evaluation of PacBio HiFi sequencing data. Applied to Candida albicans, the pipeline created a high quality assembly characterized by strong contiguity, high completeness and base level accuracy. Integrating multiple tools enabled a comprehensive assessment of genome structure, content, and k-mer derived accuracy. Outputs demonstrated that the pipeline was able to reliably construct Candida albicans with minimal fragmentation and robust agreement between reads and assembly. 

The significance pipeline extends beyond a single organism. Its modular design and reliance on long read sequencing makes it broadly applicable for addressing key genomics questions. This framework provides an accessible, transparent solution that can be adopted in both academic and industrial settings where accurate genome assemblies are essential. 

Future improvements will focus on expanding the workflow with Hi-C scaffolding, polishing, and variant calling modules, and enhanced resource configurability to support larger and more complex genomes. Additional development may include containerizing to house the entire pipeline to improve portability and reproducibility across computing environments. With these enhancements, the workflow has the potential to produce chromosome-level assemblies and serve as a flexible, generalizable platform for a wide range of genomic applications.

## 5. References
1. Cheng, H., Concepcion, G.T., Feng, X. et al. Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nat Methods 18, 170–175 (2021). https://doi.org/10.1038/s41592-020-01056-5
2. P. Di Tommaso, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319 (2017) doi:10.1038/nbt.3820
3.Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
4.Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu, fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
5. Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191
6.Rhie A, Walenz BP, Koren S, Phillippy AM. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. Genome Biol. 2020 Sep 14;21(1):245. doi: 10.1186/s13059-020-02134-9. PMID: 32928274; PMCID: PMC7488777.
7.Tegenfeldt F., Kuznetsov D., Manni M., Berkeley M., Zdobnov E.M., Kriventseva E.V. OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes. Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D516–D522, https://doi.org/10.1093/nar/gkae987
8.Alla Mikheenko, Vladislav Saveliev, Pascal Hirsch, Alexey Gurevich, WebQUAST: online evaluation of genome assemblies, Nucleic Acids Research (2023) 51 (W1): W601–W606. doi: 10.1093/nar/gkad406 First published online: May 17, 2023
9. Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis results for multiple tools and samples in a single report, Bioinformatics, Volume 32, Issue 19, October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354
10. Hoyer LL, Freeman BA, Hogan EK and Hernandez AG (2024) Use of a Candida albicans SC5314 PacBio HiFi reads dataset to close gaps in the reference genome assembly, reveal a subtelomeric gene family, and produce accurate phased allelic sequences. Front. Cell. Infect. Microbiol. 14:1329438. doi: 10.3389/fcimb.2024.1329438 
      dataset: https://www.ncbi.nlm.nih.gov/sra/SRR23724250
      New C. Albicans assembly: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032688725.1 
