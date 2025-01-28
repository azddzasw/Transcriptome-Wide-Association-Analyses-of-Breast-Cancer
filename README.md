# Transcriptome-Wide-Association-Analyses-of-Breast-Cancer

### Project Directory Structure

Here is the directory structure:

```bash
./wgs_practice/
├── bin
├── input
└── output
```
- **input**: Stores all input data.
- **output**: Stores all output data.
- **bin**: Contains all executable programs and scripts.

The **output** directory only stores result data, which is generated by combining the data and scripts in the **input** and **bin** directories.

### Sequence comparison

The first step in the analysis is sequence comparison, which involves mapping sequencing reads to a reference genome to determine the genomic location of each read. For this, we use **BWA**, one of the most widely used tools for this purpose. Before performing the comparison, we need to build the FM-index (comparison index) for the reference sequence required by BWA.

```bash
$ /Tools/common/bin/bwa index BreastCancerGenome.fa
```
 After the process, five index files will be generated with the prefix `BreastCancerGenome.fa`:

```bash
BreastCancerGenome.fa.amb
BreastCancerGenome.fa.ann
BreastCancerGenome.fa.bwt
BreastCancerGenome.fa.pac
BreastCancerGenome.fa.sa
```

Now we proceed with sequence comparison using BWA, followed by format conversion, sorting, and marking of PCR duplicates using SAMtools and GATK. The steps are as follows:

```bash
#1 Sequence comparison
time /Tools/common/bin/bwa mem -t 4 -R '@RG\tID:foo\tPL:illumina\tSM:BreastCancer' /Project/wgs_practice/input/fasta/BreastCancerGenome.fa /Project/wgs_practice/input/fastq/BreastCancerGenome.fastq.gz /Project/wgs_practice/input/BreastCancerGenome.fastq.gz | /Tools/common/bin/samtools view -Sb - > /Project/wgs_practice/output/BreastCancerGenome.bam && echo "** bwa mapping done **"

#2 Sorting
time /Tools/common/bin/samtools sort -@ 4 -m 4G -O bam -o /Project/wgs_practice/output/BreastCancerGenome.sorted.bam /Project/wgs_practice/output/BreastCancerGenome.bam && echo "** BAM sort done"

rm -f /Project/wgs_practice/output/BreastCancerGenome.bam

#3 Marking PCR duplicates
time /Tools/common/bin/gatk/4.0.1.2/gatk MarkDuplicates -I /Project/wgs_practice/output/BreastCancerGenome.sorted.bam -O /Project/wgs_practice/output/BreastCancerGenome.sorted.markdup.bam -M /Project/wgs_practice/output/BreastCancerGenome.markdup_metrics.txt && echo "** markdup done **"

#4 Deleting unnecessary files (optional)
rm -f /Project/wgs_practice/output/BreastCancerGenome.bam
rm -f /Project/wgs_practice/output/BreastCancerGenome.bam

#5 Creating index
time /Tools/common/bin/samtools index /Project/wgs_practice/output/BreastCancerGenome.sorted.markdup.bam && echo "** index done **"
```

### Explanation 

The script, named `bwa_and_markdup.sh`, is stored in the `bin` directory for better organization. Here’s a brief breakdown of its operations:

1. **Sequence comparison**: The `bwa mem` command maps the quality-controlled reads of Breast Cancer Genome to its reference genome. The `-t` option sets the number of threads (4 in this case). The data is piped (`|`) to SAMtools for conversion to the binary BAM format, and the result is saved to a file.

   > **Note**: The `-R` parameter specifies the read group (RG) information, which is crucial for downstream processes like GATK. It includes `ID` (identifier), `PL` (platform), and `SM` (sample name). For compatibility with GATK, the platform (`PL`) must be one of the following: `ILLUMINA`, `SLX`, `SOLEXA`, `SOLID`, `454`, `LS454`, `COMPLETE`, `PACBIO`, `IONTORRENT`, `CAPILLARY`, `HELICOS`, or `UNKNOWN`.

2. **Sorting**: The BAM file is sorted by genomic position using SAMtools. Sorted files are required for most downstream analyses.

3. **Marking PCR Duplicates**: PCR duplicates are identified and marked using GATK's `MarkDuplicates` tool.

4. **Indexing**: An index is generated for the BAM file using SAMtools to allow efficient access to specific regions of the genome during analysis.

The `time` command at the beginning of each step records the runtime, which can be useful for workflow optimization.

---

### Variant detection

Next, we perform variant detection using GATK. Before starting, we need to generate a `.dict` file for the reference genome using the `CreateSequenceDictionary` module:

```bash
$ /Tools/common/bin/gatk/4.0.1.2/gatk CreateSequenceDictionary -R BreastCancerGenome.fa -O BreastCancerGenome.dict && echo "** dict done **"
```

The `.dict` file must have the same prefix as the FASTA file and be located in the same directory for GATK to recognize it.

#### Steps for Variant Calling

1. **Generate GVCF**: The first step is to generate a GVCF (Genomic VCF) file using GATK’s `HaplotypeCaller`:

   ```bash
   time /Tools/common/bin/gatk/4.0.1.2/gatk HaplotypeCaller \
     -R /Project/wgs_practice/input/fasta/BreastCancerGenome.fa \
     --emit-ref-confidence GVCF \
     -I /Project/wgs_practice/output/BreastCancerGenome.sorted.markdup.bam \
     -O /Project/wgs_practice/output/BreastCancerGenome.g.vcf && echo "** gvcf done **"
   ```

2. **Joint Genotyping**: To improve the accuracy of variant detection, multiple GVCF files can be jointly genotyped using `GenotypeGVCFs`. Here, we demonstrate it for a single sample:

   ```bash
   time /Tools/common/bin/gatk/4.0.1.2/gatk GenotypeGVCFs \
     -R /Project/wgs_practice/input/BreastCancerGenome.fa \
     -V /Project/wgs_practice/output/BreastCancerGenome.g.vcf \
     -O /Project/wgs_practice/output/BreastCancerGenome.vcf && echo "** vcf done **"
   ```

Finally, the resulting VCF file is compressed and indexed using `bgzip` and `tabix` for easier downstream analysis:

```bash
# Compress VCF
time /Tools/common/bin/bgzip -f /Project/wgs_practice/output/BreastCancerGenome.vcf

# Index VCF
time /Tools/common/bin/tabix -p vcf /Project/wgs_practice/output/BreastCancerGenome.vcf.gz
```

---

#### Hard Filtering  

The following is  hard filtering, which  conducted using the latest version of GATK. In GATK 4.0, there is a dedicated `VariantFiltration` module (inherited from GATK 3.x) that makes this process straightforward. However, when applying filtering, it is essential to separate SNPs and Indels as they require different thresholds for filtering.  The detailed execution commands are as follows:  

```bash
# Select SNPs using SelectVariants  
time /Tools/common/bin/gatk/4.0.1.2/gatk SelectVariants \
    -select-type SNP \
    -V ../BreastCancerGenome.vcf.gz \
    -O ../BreastCancerGenome.snp.vcf.gz  

# Apply hard filtering for SNPs  
time /Tools/common/bin/gatk/4.0.1.2/gatk VariantFiltration \
    -V ../BreastCancerGenome.snp.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O ../BreastCancerGenome.snp.filter.vcf.gz  

# Select Indels using SelectVariants  
time /Tools/common/bin/gatk/4.0.1.2/gatk SelectVariants \
    -select-type INDEL \
    -V ../BreastCancerGenome.vcf.gz \
    -O ../BreastCancerGenome.indel.vcf.gz  

# Apply filtering for Indels  
time /Tools/common/bin/gatk/4.0.1.2/gatk VariantFiltration \
    -V ../BreastCancerGenome.indel.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "Filter" \
    -O ../BreastCancerGenome.indel.filter.vcf.gz  

# Merge the filtered SNP and Indel variants  
time /Tools/common/bin/gatk/4.0.1.2/gatk MergeVcfs \
    -I ../BreastCancerGenome.snp.filter.vcf.gz \
    -I ../BreastCancerGenome.indel.filter.vcf.gz \
    -O ../BreastCancerGenome.filter.vcf.gz  

# Remove unnecessary intermediate files  
rm -f ../output/BreastCancerGenome.snp.vcf.gz* ../output/BreastCancerGenome.snp.filter.vcf.gz* ../output/BreastCancerGenome.indel.vcf.gz* ../output/BreastCancerGenome.indel.filter.vcf.gz*  
```

