This file URL: https://github.com/JetBrains-Research/sysbio-workshop/blob/master/2019_bioinf_summer_school/chipseq/commands.md

Monocytes CD14 session: https://artyomovlab.wustl.edu/publications/supp_materials/SBW-epigenetics-StL/sessions/human_cd14_reads_cov.yaml




Run pipeline:
```bash
conda activate snakemake
snakemake --use-conda -pr --ri --cores 2
```

# Explore results

Change directory:
```bash
cd ~/chipseq/workdir
```

Activate `bio` conda environment
```bash
cd ~/chipseq/workdir
```

View files in current directory and subdirectories
```bash
tree
```
# Explore aligned files (first 4 lines)
```bash
samtools view bams/GSM1102797_CD14_H3K4me3_hg19.chr15.bam | head -n 4
```

# Explore MACS2 broad results

**Find macs2 results:**
```bash
find . -name "*.broadPeak"
tree | grep "broadPeak"
```

**Preview peaks file content (first 4 lines):**
```bash
head -n 4 ./macs2/GSM1102797_CD14_H3K4me3_hg19.chr15_broad0.1_peaks.broadPeak
```

**Most confident peaks (by 9-th column):**
```bash
cat ./macs2/GSM1102797_CD14_H3K4me3_hg19.chr15_broad0.1_peaks.broadPeak | sort -k9,9r | head -n 4
```

**Max peak length**:
```bash
cat ./macs2/GSM1102797_CD14_H3K4me3_hg19.chr15_broad0.1_peaks.broadPeak | awk 'BEGIN { max_len=0 }; { len = $3-$2; if (len > max_len) max_len = len } END { print "Max peak length:", max_len }'
```

**Min peak length**:
```bash
cat ./macs2/GSM1102797_CD14_H3K4me3_hg19.chr15_broad0.1_peaks.broadPeak | awk 'BEGIN { min_len=99999999 }; {if (NF > 0) { len = $3-$2; if (len < min_len) min_len = len } }  END { print "Min peak length",  min_len }'
```

**Avg peak length**:
```bash
cat ./macs2/GSM1102797_CD14_H3K4me3_hg19.chr15_broad0.1_peaks.broadPeak | awk 'BEGIN { sum_len=0 }; { if (NF > 0) { sum_len += $3-$2 } }; END { print "Avg peak length:", sum_len / NR} '
```

# Explore SICER results

**Preview peaks file content (first 4 lines):**
```bash
head -n 4 sicer/GSM1102797_CD14_H3K4me3_hg19.chr15-W200-G600-islands-summary-FDR0.01
```

# MACS2 vs SICER
**Peaks number in some file:**
```bash
cat macs2/GSM1102807_CD14_input_hg19.chr15_broad0.1_peaks.broadPeak | wc -l
```

**File with biggest number of peaks:**
```bash
ls macs2/*.broadPeak sicer/* | xargs wc -l | grep -v total | sort -k1,1nr
```

# Explore SPAN results

```bash
head -n 4 span/GSM1102797_CD14_H3K4me3_hg19.chr15_200_1E-6_5.peak
```

# Interpretation & Playground

**Create directory for results:**
```bash
mkdir -p downstream
cd downstream
```

**Switch to 'bio' conda environment:**
```bash
conda activate bio
```

**Example:**

Sort BED file by chromosome and by start position (sometimes could be useful)
```bash
cat ../sicer/GSM1102797_CD14_H3K4me3_hg19.chr15-W200-G600-islands-summary-FDR0.01 | sort -k1,1 -k2,2n > sorted.bed
```

**Download genes annotations GTF file:**

Already downloaded from https://www.gencodegenes.org/human/release_31lift37.html to
`~/chipseq/tracks/gencode.v31lift37.annotation.gtf`

**GTF to BED file:**

* FILTER out all the non-chromosome positions
* ONLY 15 chromosome is taken into account
* FILTER out the transcripts
* Convert GTF file to TAB-separated BED file
* SORT the result

```bash
cat ~/chipseq/tracks/gencode.v31lift37.annotation.gtf | grep -v "#" | grep "^chr15" | awk -v OFS='\t' '($3=="gene") {print $1,$4-1,$5,$10}' | sort -k1,1 -k2,2n > gencode.v31lift37.annotation.bed
```
```bash
head gencode.v31lift37.annotation.bed 
```

**Find the closest genes for the peak file:**
```bash
bedtools closest -a ../macs2/GSM1102797_CD14_H3K4me3_hg19.chr15_broad0.1_peaks.broadPeak -b gencode.v31lift37.annotation.bed -D ref | head -n 1
```

**Plot signal profile (deeptools):**

* Step 1 - matrix computation:
    ```bash
    computeMatrix reference-point -S ../bw/GSM1102797_CD14_H3K4me3_hg19.chr15.bw -R  gencode.v31lift37.annotation.bed -a 3000 -b 3000 -out matrix.tsv.gz
    ```
* Step2 - plot profile
    ```bash
    plotProfile -m matrix.tsv.gz -out TssProfile.png --plotTitle "TSS K4me3 profile"
    ```
* Step3 - plot heatmap
    ```bash
    plotHeatmap -m matrix.tsv.gz --outFileName TssHeatmap.png --plotTitle "TSS k4me3 coverage heatmap"
    ```

**Function annotation:**

ChIPSeeqer R code: open R studion, see `commands.r` file

**Make peaks file GREAT again:**

Preparing my peaks file for uploading to GREAT, picking first 3 columns:

H3K4me3 (promoter), first & last 200 peaks :
```bash
cat ../macs2/GSM1102797_CD14_H3K4me3_hg19.chr15_broad0.1_peaks.broadPeak | awk -v OFS='\t' '{print($1,$2,$3)}' > h3k4me3.bed
```

Pay attention to `>` vs `>>`:
```bash
head -n 100 h3k4me3.bed > h3k4me3_200.bed
tail -n 100 h3k4me3.bed >> h3k4me3_200.bed
```

Or Take 100 top significant H3K27ac peaks (~active enhancers)
```bash
cat ../macs2/GSM1102782_CD14_H3K27ac_hg19.chr15_broad0.1_peaks.broadPeak |  sort -k9,9r | awk -v OFS='\t' '{print($1,$2,$3)}' > h3k27ac.bed
head -n 100 h3k27ac.bed > h3k27ac_100.bed
```

Try filtered peaks in GREAT w/o background or all peaks as background
