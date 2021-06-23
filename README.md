# EPITRAN workshop -RNA Editing 24/06/21 
<a href="mailto:ernesto.picardi@uniba.it ">Ernesto Picardi</a><br>
<a href="mailto:claudio.logiudice@uniba.it ">Claudio Lo Giudice</a>
<p align="justify"> Epitranscriptome refers to all biochemical RNA modifications affecting cellular RNAs. 
To date, >150 different RNA modifications have been discovered, and, among them, RNA editing represents the most prominent post-transcriptional mechanism.<br> 
In mammals, it involves the deamination of cytosine to uridine or, more commonly, the conversion of adenosine to inosine (A-to-I) and, concurrently with alternative splicing, increases the diversity of the eukaryotic transcriptome and proteome.<br>
Nowadays, RNAseq is the de facto standard approach to discover RNA editing candidates in whole eukaryotic genomes.
Although the identification of editing sites is, in principle, quite simple, it represents a computational challenge because true RNA editing events have to be discriminated from genome-encoded SNPs and technical artefacts caused by sequencing or read-mapping errors.<br>
The use of genomic reads from whole genome sequencing (WGS) or whole exome sequencing (WXS) experiments in single individuals, annotations in dbSNP
and several stringent filters can minimize the detection of false RNA editing candidates
To promote and facilitate the investigation of RNA editing at genomic scale, we developed and released the first software package devoted to this purpose, named <b>REDItools</b>. <br>
It comprises a suite of python scripts has been conceived to handle massive transcriptome sequencing data through a variety of filters to provide accurate RNA editing calls.<br>
Reditools have been applied in thousands of human RNAseq experiments, leading to the discovery of >16 million events in several human body sites, collected in our specialized <b>REDIportal</b> database.<br>
REDItools and REDIportal represent two relevant and efficient computational resources to investigate RNA editing in a variety of organisms and in different physiological and pathological conditions.<br> 
Today, we are going to reproduce step-by-step a bioinformatics pipeline aimed to study A-to-I editing in human samples, combininig both Reditools and Rediportal information.
During the morning we will start with exploring the REDIportal database and we will conclude our practice by obtaining the annotated editing events for two genes
GRIA2 and FLNA. These data will by used as input for one of the script embedded with the REDItools package to explore the RNA editing potential of two different RNA-Seq data sets using known editing events.
Finally, during the afternoon we will illustrate a computational workflow to discover novel RNA editing events using RNAseq and WGS reads from the same sample.
</p>

<p> 
The main steps described during the practice are reported below and can be easily copy/pasted in your terminal.<br>
<b>Note.</b> Assuming you're traineeX, please change X according to your workspace.<br>
<b>Note2.</b> Choose a GTEX sample from Cerebellum or Lung and copy DNAseq/RNAseq accordingly in your home folder.<br>
<b>IMPORTANT!</b> REDItoolDnaRna.py outTable (eg. outTable_892028847) contains 8digit random number, so it usually varies among users and different script launches on the same machine.
  
<table>
<thead>
<th>Sample</th>
<th>Tissue</th>
<th>Gender</th>
<th>Path</th>
<tr>
<td>SRR1307770.bam</td>
<td>Brain - Cerebellum</td>
<td>male normal</td> 
<td>/usr/share/course_data/rnaediting/cerebellum/RNAseq_fq</td>
</tr>
<tr>
<td>SRR1434867.bam</td>
<td>Brain - Cerebellum</td>
<td>male normal</td>
<td>/usr/share/course_data/rnaediting/cerebellum/RNAseq_fq</td>
</tr>
<tr>
<td>SRR1361736.bam</td>
<td>Brain - Cerebellum</td>
<td>male normal</td>
<td>/usr/share/course_data/rnaediting/cerebellum/RNAseq_fq</td>
</tr>
<tr>
<td>SRR1310520.bam</td>
<td>Lung</td>
<td>male normal</td>
<td>/usr/share/course_data/rnaediting/lung/RNAseq_fq</td>
</tr>
<tr>
<td>SRR1334490.bam</td>
<td>Lung</td>
<td>male normal</td>
<td>/usr/share/course_data/rnaediting/lung/RNAseq_fq</td>
</tr>
<tr>
<td>SRR1475524.bam</td>
<td>Lung</td>
<td>male normal</td>
<td>/usr/share/course_data/rnaediting/lung/RNAseq_fq</td>
</tr>
</thead>
</table>
</p>

<pre>
<b>1) Create a folder for your RNAseq data (eg. RNAseq)</b>
$ mkdir RNAseq

<b>2) Copy the RNAseq data inside it</b>
$ cp ../data/rnaediting/lung/RNAseq_fq/SRR1310520_r* .

<b>3) Align RNASeq reads to the reference genome with STAR:</b>
$ STAR --runThreadN 4   --genomeDir /usr/share/course_data/STAR_genome_index_ucsc/   --genomeLoad NoSharedMemory   --readFilesIn ./SRR1310520_r1.fastq.gz   ./SRR1310520_r2.fastq.gz  --readFilesCommand zcat  --outFileNamePrefix SRR1310520_chr21_   --outReadsUnmapped Fastx   --outSAMtype BAM   SortedByCoordinate      --outSAMstrandField intronMotif   --outSAMattributes All      --outFilterType BySJout   --outFilterMultimapNmax 1  --outFilterMismatchNmax 999   --outFilterMismatchNoverLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000  --alignMatesGapMax 1000000   --alignSJoverhangMin 8   --alignSJDBoverhangMin 1

$cd ..

<b>4) Create a folder for your DNAseq data (eg. DNAseq)</b>
$ mkdir DNAseq

>b>5) Copy the DNAseq data inside it</b>
$ cp ../data/rnaediting/DNAseq/Lung* .

<b>6) Detect all potential DNA–RNA variants in your data (limited to chromosome 21) using the REDItoolDnaRNA.py script:</b>

$ REDItoolDnaRna.py -i ./RNAseq/SRR1310520_chr21_Aligned.sortedByCoord.out.bam -j ./DNAseq/Lung_sorted.bam -o editing -f /usr/share/course_data/rnaediting/hg19ref/GRCh37.primary_assembly.genome.fa  -c1,1 -m30,255 -v1 -q30,30 -e -n0.0 -N0.0 -u -l -p -s2 -g2  -S -Y chr21:1-48129895
For detailed REDItoolDnaRnaGTEX.py options <a href="https://github.com/BioinfoUNIBA/REDItools/blob/master/README_1.md#reditooldnarna-py">click here</a>

<b>7) Exclude invariant positions as well as positions not supported by ≥10 WGS reads:</b>

$ awk 'FS="\t" {if ($8!="-" && $10>=10 && $13=="-") print}' editing/DnaRna_892028847/outTable_892028847 > outTable_892028847_chr21.out

<b>8) Annotate positions using RepeatMasker and dbSNP annotations:</b>

$ AnnotateTable.py -a /usr/share/course_data/rnaediting/rptmsk/rmsk_chr21.sorted.gtf.gz -n rmsk -i outTable_892028847_chr21.out -o outTable_892028847_chr21.out.rmsk -u

$ AnnotateTable.py -a /usr/share/course_data/rnaediting/dbsnp/snp151_chr21.sorted.gtf.gz -n snp151 -n snp151 -i outTable_892028847_chr21.out.rmsk -o outTable_892028847_chr21.out.rmsk.snp -u
For detailed AnnotateTable.py options <a href="https://github.com/BioinfoUNIBA/REDItools/blob/master/README_1.md#annotatetable-py">click here</a>

<b>9) Create a first set of positions selecting sites supported by at least five RNAseq reads and a single mismatch:</b>

$ selectPositions.py -i outTable_892028847_chr21.out.rmsk.snp -c 5 -v 1 -f 0.0 -o outTable_892028847_chr21.out.rmsk.snp.sel1
For detailed selectPositions.py options <a href="https://github.com/BioinfoUNIBA/REDItools/blob/master/README_1.md#selectpositions-py">click here</a>

<b>10) Create a second set of positions selecting sites supported by ≥10 RNAseq reads, three mismatches and minimum editing frequency of 0.1: </b>

$ selectPositions.py -i outTable_892028847_chr21.out.rmsk.snp -c 10 -v 3 -f 0.1  -o outTable_892028847_chr21.out.rmsk.snp.sel2

<b>11) Select ALU sites from the first set of positions:</b>

$ awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)=="Alu" && $17=="-" && $8!="-" && $10>=10 && $13=="-") print}' outTable_892028847_chr21.out.rmsk.snp.sel1 > outTable_892028847_chr21.out.rmsk.snp.alu

<b>12) Select REP NON ALU sites from the second set of positions, excluding sites in Simple repeats or Low complexity regions:</b>

$ awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15!="-" && $15!="Simple_repeat" && $15!="Low_complexity" && $17=="-" && $8!="-" && $10>=10 && $14<=0.05 && $9>=0.1) print}' outTable_892028847_chr21.out.rmsk.snp.sel2 > outTable_892028847_chr21.out.rmsk.snp.nonalu

<b>13) Select NON REP sites from the second set of positions:</b>

$ awk 'FS="\t" {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15=="-" && $17=="-" && $8!="-" && $10>=10 && $14<=0.05 && $9>=0.1) print}' outTable_892028847_chr21.out.rmsk.snp.sel2 > outTable_892028847_chr21.out.rmsk.snp.nonrep

<b>14) Annotate ALU, REP NON ALU and NON REP sites using known editing events from REDIportal:</b>

$ AnnotateTable.py -a /usr/share/course_data/rnaediting/rediportal/atlas.gtf.gz -n ed -k R  -c 1 -i outTable_892028847_chr21.out.rmsk.snp.alu -o outTable_892028847_chr21.out.rmsk.snp.alu.ed -u

$ AnnotateTable.py -a /usr/share/course_data/rnaediting/rediportal/atlas.gtf.gz -n ed -k R  -c 1 -i outTable_892028847_chr21.out.rmsk.snp.nonalu -o outTable_892028847_chr21.out.rmsk.snp.nonalu.ed -u

$ AnnotateTable.py -a /usr/share/course_data/rnaediting/rediportal/atlas.gtf.gz -n ed -k R  -c 1 -i outTable_892028847_chr21.out.rmsk.snp.nonrep -o outTable_892028847_chr21.out.rmsk.snp.nonrep.ed -u

<b>15) Extract known editing events from ALU, REP NON ALU and NON REP sites:</b>

$ mv outTable_892028847_chr21.out.rmsk.snp.alu.ed alu

$ mv outTable_892028847_chr21.out.rmsk.snp.nonalu.ed nonalu

$ mv outTable_892028847_chr21.out.rmsk.snp.nonrep.ed nonrep

$ cat alu nonalu nonrep > alu-nonalu-nonrep

$ awk 'FS="\t" {if ($19=="ed") print}' alu-nonalu-nonrep > knownEditing 

<b>16) Convert editing candidates ($19!="ed" selects novel RNA editing events.) in REP NON ALU and NON REP sites in GFF format for further filtering:</b>

$ cat nonalu nonrep > nonalu-nonrep

$ awk 'FS="\t" {if ($19!="ed") print}' nonalu-nonrep > pos.txt

$TableToGFF.py -i pos.txt -s -t -o pos.gff
For detailed TableToGFF.py options <a href="https://github.com/BioinfoUNIBA/REDItools/blob/master/README_1.md#tabletogff-py-new-in-version-1-0-3">click here</a>

<b>17) Convert editing candidates in ALU sites in GFF format for further filtering:</b>

$ awk 'FS="\t" {if ($19!="ed") print}' alu > posalu.txt

$ TableToGFF.py -i posalu.txt -s -t -o posalu.gff

<b>18) Launch REDItoolDnaRna.py on ALU sites using stringent criteria to recover potential editing candidates:</b>

$ REDItoolDnaRna.py -s 2 -g 2 -S -t 4 -i ./RNAseq/SRR1310520_chr21_Aligned.sortedByCoord.out.bam -f /usr/share/course_data/rnaediting/hg19ref/GRCh37.primary_assembly.genome.fa -c 5,5 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 1 -n 0.0 -e -T posalu.sorted.gff.gz -w /usr/share/course_data/rnaediting/Gencode_annotation/gencode.v30lift37.chr21.splicesites.txt -k /usr/share/course_data/rnaediting/hg19ref/nochr -R -o firstalu

<b>19) Launch REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria to recover RNAseq reads harboring reference mismatches:</b>

$ REDItoolDnaRnaGTEX.py -s 2 -g 2 -S -t 4 -i ./RNAseq/SRR1310520_chr21_Aligned.sortedByCoord.out.bam -f /usr/share/course_data/rnaediting/hg19ref/GRCh37.primary_assembly.genome.fa -c 10,10 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T pos.sorted.gff.gz -w /usr/share/course_data/rnaediting/Gencode_annotation/gencode.v30lift37.chr21.splicesites.txt -k /usr/share/course_data/rnaediting/hg19ref/nochr --reads -R --addP -o first

<b>20) Launch pblat on RNAseq reads harboring reference mismatches from previous step and select multimapping reads:</b>

$ pblat -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 /usr/share/course_data/rnaediting/hg19ref/GRCh37.primary_assembly.genome.fa first/DnaRna_304977045/outReads_304977045 reads.psl

$ readPsl.py reads.psl badreads.txt

<b>21) Extract RNAseq reads harboring reference mismatches from Step 19 and remove duplicates:</b>
  
$ sort -k1,1 -k2,2n -k3,3n first/DnaRna_595685947/outReads_595685947 | mergeBed > bed 

$ samtools view -@ 4 -L bed -h -b ./RNAseq/SRR1310520_chr21_Aligned.sortedByCoord.out.bam > SRR1310520_chr21_bed.bam

$ samtools sort -@ 4 -n SRR1310520_chr21_bed.bam -o SRR1310520_chr21_bed_ns.bam 

$ samtools fixmate -@ 4 -m SRR1310520_chr21_bed_ns.bam  SRR1310520_chr21_bed_ns_fx.bam 

$ samtools sort -@ 4 SRR1310520_chr21_bed_ns_fx.bam -o SRR1310520_chr21_bed_ns_fx_st.bam

$ samtools markdup -r -@ 4 SRR1310520_chr21_bed_ns_fx_st.bam SRR1310520_chr21_bed_dedup.bam

$ samtools index SRR1310520_chr21_bed_dedup.bam

<b>22) Re-run REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria, deduplicated reads and mis-mapping info:</b>

$ REDItoolDnaRna.py -s 2 -g 2 -S -t 4 -i SRR1310520_chr21_bed_dedup.bam -f /usr/share/course_data/rnaediting/hg19ref/GRCh37.primary_assembly.genome.fa -c 10,10 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T pos.sorted.gff.gz -w /usr/share/course_data/rnaediting/Gencode_annotation/gencode.v30lift37.chr21.splicesites.txt -R -k /usr/share/course_data/rnaediting/hg19ref/nochr -b badreads.txt --rmIndels -o second
  
<b>23) Collect filtered ALU, REP NON ALU and NON REP sites:</b>

$ collect_editing_candidates.py

$ sort -k1,1 -k2,2n editing.txt > editing_sorted.txt

<b>24) Inspect the distribution of editing candidates to look at A-to-I enrichment: </b>

$ get_Statistics.py
  
