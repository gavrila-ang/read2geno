#!/usr/bin/env python

""" 
Palmer Lab Genotyping and Validation Pipeline
Step 3: Quality control of metadata file
Gracie Ang, Ben Johnson, Denghui Chen (est. 2021)

Inputs:
    - Sample processing input file
    - Sample processing array file
    - Reference genome
    - Raw, demultiplexed sample-level FASTQs 

Resources
    - Linux High-Performance Clusters
    - Slurm Workload Manager 

Outputs:
    - Cleaned sample-level FASTQs
    - Alignemnt files (.bam)
    - Alignment quality control files (.csv, .png, .fastqc, .flagstat, etc.)

Assumptions:
    - Custom input file directory structure

"""

##############################################################################################################################################################################

# === Load python packages ===
import user
import subprocess
import pandas as pd
import sys
import os
from time import time

# === Export user arguments as global variables ===
user_vars_dict = user.user_args.import_user_config()
globals().update(user_vars_dict)
softwares_dict = user.user_args.import_softwares()
globals().update(softwares_dict)

# === Load required modules (CPU, GCC, Java) for picard ===
subprocess.run("module load cpu/0.17.3 && module load gcc/10.2.0-2ml3m2l && module load openjdk/1.8.0_265-b01-7lpdkk2", shell=True)

# === Load system arguments ===
array_taskid = int(sys.argv[1]) 
ncpu = int(sys.argv[2])
mem_per_cpu = int(sys.argv[3])
javamem = int(ncpu * mem_per_cpu)
javamem = f"{javamem}G"

print(f"array_taskid : {array_taskid}")
print(f"ncpu : {ncpu}")
print(f"mem_per_cpu : {mem_per_cpu}")
print(f"javamem : {javamem}")

# === Declare data type of special columns ===
dtype_mapping={'rfid':str}

# === Extract info for the specific sample_id being processed ===
input_filepath = f"{round_directory}/output/03_sample_processing/{current_round}_sample_processing_input.csv"
sample_processing_input = pd.read_csv(input_filepath, dtype=dtype_mapping).reset_index(drop=True)
row = sample_processing_input.loc[sample_processing_input.idx == array_taskid].iloc[0]

##############################################################################################################################################################################

# === Runner ===

STEP = 0
elapsed_total = 0

def run_step(name, func):
    global STEP
    STEP += 1
    print(f"\n[{STEP:02d}] {name}")
    start = time()
    func()

    # Calculate time in hours and minutes
    elapsed = time() - start
    print(elapsed)

def sh(cmd, stdout_path=None):
    """
    cmd: list[str]
    """    
    print(cmd)
    if stdout_path:
        with open(stdout_path, "w") as out:
            subprocess.run(cmd, check=True, stdout=out)
    else:
        subprocess.run(cmd, check=True)

##############################################################################################################################################################################

def riptide_block():

    # ===  Runs BBDuk to remove adapter sequences from paired-end reads using k-mer trimming ===
    # === k=23: Use 23‑nt k-mers to detect adapters.
    # === mink=11: Allow shorter k-mers down to 11 nt at read ends for partial matches.
    # === hdist=1: Permit 1 mismatch when matching k-mers to adapters.
    # === Trims poly-G tails ≥50 nt, and ensures both reads in a pair are trimmed consistently ===

    run_step(
        "BBDuk removal of adapter sequences and polyg tails",
        lambda: sh([
            f"{bbmap}/bbduk.sh",
            f"ref={bbmap}/resources/adapters.fa",
            f"in1={row.demuxed_fq1_filepath}",
            f"in2={row.demuxed_fq2_filepath}",
            f"out1={row.output_directory}/bbduk/{row.Sample_ID}_adapter_trimmed_R1.fastq.gz",
            f"out2={row.output_directory}/bbduk/{row.Sample_ID}_adapter_trimmed_R2.fastq.gz",
            "ktrim=r", "k=23", "mink=11", "hdist=1",
            "trimpolyg=50", "tpe", "tbo",
        ])
    )

    # === Runs cutadapt to quality-trim paired-end reads at Q5, discard reads shorter than 70 bp === 
    # === Filter pairs if either read fails, use all CPU cores ===

    run_step(
        "Cutadapt removal of low quality reads and short reads",
        lambda: sh([
            "cutadapt",
            "-q", "5",
            "--minimum-length", "70",
            "--pair-filter=any",
            "--cores=0",
            "-o", f"{row.output_directory}/cutadapt/{row.Sample_ID}_qualen_trimmed_R1.fastq.gz",
            "-p", f"{row.output_directory}/cutadapt/{row.Sample_ID}_qualen_trimmed_R2.fastq.gz",
            f"{row.output_directory}/bbduk/{row.Sample_ID}_adapter_trimmed_R1.fastq.gz",
            f"{row.output_directory}/bbduk/{row.Sample_ID}_adapter_trimmed_R2.fastq.gz",
        ])
    )

    # === Runs FastQC with multiple threads to quietly generate QC reports for the raw, adapter-trimmed ===
    # === and quality/length-trimmed FASTQ files ===
    run_step(
        "FastQC paired reads before and after trimming",
        lambda: sh([
            fastqc,
            "-t", str(ncpu),
            "-q",
            "-o", f"{row.output_directory}/fastqc_trim/",
            row.demuxed_fq1_filepath,
            row.demuxed_fq2_filepath,
            f"{row.output_directory}/cutadapt/{row.Sample_ID}_qualen_trimmed_R1.fastq.gz",
            f"{row.output_directory}/cutadapt/{row.Sample_ID}_qualen_trimmed_R2.fastq.gz",
            f"{row.output_directory}/bbduk/{row.Sample_ID}_adapter_trimmed_R1.fastq.gz",
            f"{row.output_directory}/bbduk/{row.Sample_ID}_adapter_trimmed_R2.fastq.gz",
        ])
    )

    # === Runs bwa mem to align paired-end, quality- and length-trimmed reads to the reference genome using multiple threads, marks secondary/split alignments === 
    # === Adds detailed read-group metadata, and writes the output to a SAM file. ===
    run_step(
        "BWA mem paired read alignment to ref genome",
        lambda: sh([
            bwa, "mem",
            "-aM",
            "-t", str(ncpu),
            "-R", f"@RG\\tID:{row.instrument_name}.{row.run_id}.{row.flowcell_id}.{row.flowcell_lane}\\tLB:{row.Library_ID}\\tPL:ILLUMINA\\tSM:{row.Sample_ID}\\tPU:{row.flowcell_id}.{row.flowcell_lane}.{row.Sample_Barcode}",
            reference_genome,
            f"{row.output_directory}/cutadapt/{row.Sample_ID}_qualen_trimmed_R1.fastq.gz",
            f"{row.output_directory}/cutadapt/{row.Sample_ID}_qualen_trimmed_R2.fastq.gz",
        ], stdout_path=f"{row.output_directory}/bwamem/{row.Sample_ID}.sam")
    )


def gbs_block():

    # === Runs cutadapt to remove a 5′ sample barcode from reads (allowing 25% mismatches) ===
    run_step(
        "Cutadapt single reads 5 prime end",
        lambda:sh([
            "cutadapt",
            row.demuxed_fq1_filepath,
            "-o", f"{row.output_directory}/cutadapt/{row.Sample_ID}_adapter_trimmed.fq.gz",
            "-g", f"^{row.Sample_Barcode}",
            "-e", "0.25",
            "--untrimmed-output",
            row.demuxed_fq1_filepath, 
        ])
    )

    # === Runs cutadapt to trim a 3′ adapter sequence, requiring ≥8 bp overlap === 
    # === Discarding reads <25 bp, quality-trimming at Q20, allowing 25% mismatches ===
    run_step(
        "Cutadapt single reads 3 prime end",
        lambda:sh([
            "cutadapt",
            f"{row.output_directory}/cutadapt/{row.Sample_ID}_adapter_trimmed.fq.gz",
            "-a", "AGATCGGAAGAGCGGTTCAGCAGGAATGCCG",
            "-O", "8",
            "-m", "25",
            "-q", "20",
            "-e", "0.25",
            "-o", f"{row.output_directory}/cutadapt/{row.Sample_ID}_qualen_trimmed.fq.gz",            
        ])
    )

    # === Runs FastQC with multiple threads to quietly generate QC reports for the raw, adapter-trimmed ===
    # === and quality/length-trimmed FASTQ files ===
    run_step(
        "FastQC single reads before and after trimming",
        lambda:sh([
            fastqc,
            "-t", str(ncpu),
            "-q",
            "-o", f"{row.output_directory}/fastqc_trim/",
            row.demuxed_fq1_filepath,
            f"{row.output_directory}/cutadapt/{row.Sample_ID}_adapter_trimmed.fq.gz",
            f"{row.output_directory}/cutadapt/{row.Sample_ID}_qualen_trimmed.fq.gz",
        ])
    )

    # === Runs bwa mem to align reads to the reference genome using multiple threads, marking shorter split hits, filtering low-quality alignments (MAPQ ≥20) ===
    # === Adding read-group metadata, and writing the output SAM file. ===
    run_step(
        "BWA mem single read alignment to ref genome",
        lambda:sh([
            bwa, "mem",
            "-t", str(ncpu),
            "-M",
            "-T", "20",
            "-R", f"@RG\tID:{row.Sample_ID}\tLB:{row.Library_ID}\tPL:ILLUMINA\tSM:{row.Sample_ID}\tPU:unit1",
            reference_genome,
            row.demuxed_fq1_filepath,
        ], stdout_path=f"{row.output_directory}/bwamem/{row.Sample_ID}.sam")
    )



def alignqc_block():

    # === Sorts the BAM/SAM file by genomic coordinates ===
    run_step(
        "Samtools sort alignment file",
        lambda: sh([
            samtools, "sort",
            "-@", str(ncpu),
            "-o", f"{row.output_directory}/samtools_sort/{row.Sample_ID}_sorted.bam",
            f"{row.output_directory}/bwamem/{row.Sample_ID}.sam"
        ])
    )

    # === Identifies PCR and optical duplicates in the alignment file ===
    # === Does not remove duplicates (--REMOVE_DUPLICATES=false), only flags them ===
    # === Requires the BAM to be sorted, NEEDS to come after samtools sort ===
    run_step(
        "Mark duplicates with Picard",
        lambda: sh([
            "java", f"-Xmx{javamem}", "-XX:+AggressiveOpts", "-XX:+AggressiveHeap",
            "-jar", picard, "MarkDuplicates",
            f"--INPUT", f"{row.output_directory}/samtools_sort/{row.Sample_ID}_sorted.bam",
            "--REMOVE_DUPLICATES", "false",
            "--ASSUME_SORTED", "true",
            f"--METRICS_FILE", f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.metrics",
            f"--OUTPUT", f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam"
        ])
    )

    # === Creates a BAM index (.bai) for fast random access by position ===
    run_step(
        "Samtools index for marked-duplicate bam",
        lambda: sh([
            samtools, "index",
            f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam",
            f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bai"
        ])
    )

    # === Produces basic alignment statistics: total reads, mapped reads, properly paired reads, duplicates, etc. ===
    run_step(
        "Samtools flagstat to produce basic alignment statistics",
        lambda: sh([
            samtools, "flagstat",
            "-@", str(ncpu),
            f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam",
        ], stdout_path=f"{row.output_directory}/samtools_flagstat/{row.Sample_ID}.flagstat")
    )

    # === Reports the number of bases in mapped reads across genome ===
    # === Excludes duplicated, unmapped reads, secondary alignments, and failed quality reads ===
    # === --no-per-base disables per-base depth (faster for large genomes) ===
    # === --fast-mode speeds up processing for high-coverage data ===
    run_step(
        "Mosdepth excluding duplicates",
        lambda: sh([
            "mosdepth",
            "--threads", str(ncpu),
            "--no-per-base",
            "--fasta", reference_genome,
            "--fast-mode",
            "-F", "1796",
            f"{row.output_directory}/mosdepth/{row.Sample_ID}_excluded_duplicates",
            f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam"
        ])
    )

    # === Reports the number of bases in mapped reads across genome ===
    # === Excludes, unmapped reads, secondary alignments, and failed quality reads ===
    # === --no-per-base disables per-base depth (faster for large genomes) ===
    # === --fast-mode speeds up processing for high-coverage data ===
    run_step(
        "Mosdepth excluding unmapped reads",
        lambda: sh([
            "mosdepth",
            "--threads", str(ncpu),
            "--no-per-base",
            "--fasta", reference_genome,
            "--fast-mode",
            "-F", "772",
            f"{row.output_directory}/mosdepth/{row.Sample_ID}_excluded_unmapped",
            f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam"
        ])
    )

    # === Calculates depth only at variant positions defined in the positions bed file for each chromosome ===
    # === Ignores duplicate reads automatically because MarkDuplicates was run ===
    chrom_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14",
                "15","16","17","18","19","20","X","Y","M"]

    for chrom in chrom_list:
        run_step(
            f"Samtools depth for chromosome {chrom}",
            lambda chrom=chrom: sh([
                samtools, "depth",
                "-b", f"{round_directory}/output/03_sample_processing/chr{chrom}_variant.bed",
                f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam",
                "-o", f"{row.output_directory}/samtools_depth/{row.Sample_ID}_chr{chrom}.depth"
            ])
        )

    # === Counts the number of reads mapped to each chromosome at mapping quality ≥20 ===
    # === Useful to assess sex chromosome coverage (chrX, chrY) and total coverage distribution ===
    
    # === chrX ===
    run_step(
        "Samtools view to count reads mapped to chrX",
        lambda: sh([
            samtools, "view", "-c", "-q", "20",
            f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam",
            "chrX",
            "-o", f"{row.output_directory}/samtools_view/{row.Sample_ID}_chrX_count.csv"
        ])
    )

    # === chrY ===
    run_step(
        "Samtools view to count reads mapped to chrY",
        lambda: sh([
            samtools, "view", "-c", "-q", "20",
            f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam",
            "chrY",
            "-o", f"{row.output_directory}/samtools_view/{row.Sample_ID}_chrY_count.csv"
        ])
    )

    # === All chromosomes ===
    run_step(
        "Samtools view to count reads mapped across all chromosomes",
        lambda: sh([
            samtools, "view", "-c", "-q", "20",
            f"{row.output_directory}/picard_markdup/{row.Sample_ID}_sorted_mkDup.bam",
            "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chrX","chrY","chrM",
            "-o", f"{row.output_directory}/samtools_view/{row.Sample_ID}_allchrom_count.csv"
        ])
    )

    # === Redundant Riptide and GBS files removal ===
    files_to_remove = [
        # Riptide-specific files
        f"{row.output_directory}/bbduk/{row.Sample_ID}_adapter_trimmed_R1.fastq.gz",
        f"{row.output_directory}/bbduk/{row.Sample_ID}_adapter_trimmed_R2.fastq.gz",
        # Common files for both Riptide and GBS
        f"{row.output_directory}/bwamem/{row.Sample_ID}.sam",
        f"{row.output_directory}/samtools_sort/{row.Sample_ID}_sorted.bam",
    ]

    def remove_files(files):
        """Remove files if they exist, skipping missing files safely."""
        for file in files:
            if os.path.exists(file):
                try:
                    os.remove(file)
                    print(f"Removed: {file}")
                except Exception as e:
                    print(f"Error removing {file}: {e}")
            else:
                print(f"File not found (skipped): {file}")

    # Wrap in run_step for logging and timing
    run_step(
        "Remove redundant intermediate files",
        lambda: remove_files(files_to_remove)
    )

##############################################################################################################################################################################
            

def main():

    print(f"Running sample processing on {row.Sample_ID}....")    

    # === Pre-processing and alignment is sequencing method specific ===
    if row.seq_method in ['riptide', 'LSW_lcwgs']:
        riptide_block()
    elif row.seq_method == 'gbs':
        print(f"{WARNING}: All gbs files should have been processed.")
        # gbs_block()
    else:
        print(f"{FAIL}: Invalid sequencing method.")
        sys.exit(1)

    # === Post-processing and quality-control of alignment is common to all methods ===
    alignqc_block()


##############################################################################################################################################################################

if __name__ == "__main__":
    main()