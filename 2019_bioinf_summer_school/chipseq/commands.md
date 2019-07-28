Change directory:
```bash
cd ~/chipseq/workdir
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

# Explore SICER results

**Preview peaks file content (first 4 lines):**
```bash
head -n 4 sicer/GSM1102797_CD14_H3K4me3_hg19.chr15-W200-G600-islands-summary-FDR0.01
```

# Interpretation & Playground

**Peaks number in some file:**
```bash
cat macs2/GSM1102807_CD14_input_hg19.chr15_broad0.1_peaks.broadPeak | wc -l
```

**File with biggest number of peaks:**
```bash
ls macs2/*.broadPeak sicer/* | xargs wc -l | grep -v total | sort -k1,1nr
```

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

**Plot signal profile (deeptools):**

* Step 1 - matrix computation:
    ```bash
    computeMatrix reference-point -S ../bw/GSM1102797_CD14_H3K4me3_hg19.chr15.bw -R  gencode.v31lift37.annotation.bed -a 3000 -b 3000 -out matrix.tsv.gz
    ```
* Step2 - plot profile
    ```bash
    plotProfile -m matrix.tsv.gz -out TssProfile.png --plotTitle "TSS K4me3 profile"
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
cat ../macs2/GSM1102782_CD14_H3K27ac_hg19.chr15_broad0.1_peaks.broadPeak |  sort -k9,9r | awk -v OFS='\t' '{print($1,$2,$3)}' > k27ac.bed
head -n 100 h3k27ac.bed > h3k27ac_100.bed
```

Try filtered peaks in GREAT w/o background or all peaks as background