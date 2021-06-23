# EPITRAN workshop -RNA Editing 24/06/21 MORNING SESSION
<p class="text-justify">
Searching into <a href="http://srv00.recas.ba.infn.it/atlas/search.html">REDIportal</a> is quite straightforward and also users with no bioinformatics skills can perform accurate searches across the database. RNA editing sites are stored according to their genomic positions and can be retrieved providing a genomic locus (“Genomic Region” field) or a known gene symbol (“Gene Name” field). Both fields are mutually exclusive.
Genomic loci can be interrogated entering chromosome coordinates in the format Chr:start-end (for example chr4:158101247-158308846).<br>
<b>Note.</b> Today we are going to extract all the editing events related to the genes GRIA2 and FLNA.<br> 
The protein encoded by the <b>FLNA</b> gene is an actin-binding protein that is involved in remodeling the cytoskeleton to effect changes in cell shape and migration.<br>
<b>GRIA2</b> gene product belongs to a family of glutamate receptors.Human and animal studies suggest that pre-mRNA editing is essential for brain function, and defective GRIA2 RNA editing may be relevant to amyotrophic lateral sclerosis (ALS) etiology.

FLNA and GRIA2 annotations will be adapted to be used with the REDItoolKnown.py script (see below).<br>
This script has been developed to explore the RNA editing potential of RNA-Seq data sets using <b>only known</b> editing events.    
</p>
<br>
<div align="center"><img src="img/sc1.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
RNA editing events in known genes can be retrieved entering the gene symbol in the “Gene Name” field. To avoid editing sites in intergenic regions surrounding the entered gene name, the “Extact Match” check box must be selected. The “Gene Name” field allows the autocomplete function to facilitate the selection of right gene.
</p>
<br>
<div align="center"><img src="img/sc2.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
Once the genomic region or gene name has been entered, the search can be refined using additional select menus. The following options are admitted:
</p>

<table class="table table-bordered">
<tr>
<th>Menu</th>
<th>Name</th>
<th>Option</th>
</tr>
<tr>
<td width="30%"><img src="img/sc3.png" width="100%" height="100%"/></td>
<td style='text-align:center' width="20%"><b>Location</b></td>
<td>
<ul>
<li>ALU
<li>NONREP
<li>REP
</ul>
<p class="text-justify">
Location menu allows the selection of RNA editing sites residing in Alu elements (ALU) or repetitive elements non-Alu (REP) or non repetitive regions (NONREP).
</p>
</td>
</tr>
<tr>
<td width="30%"><img src="img/sc4.png" width="100%" height="100%"/></td>
<td style='text-align:center' width="20%"><b>Genic Region</b></td>
<td>
<ul>
<li>5'UTR
<li>3'UTR
<li>Intronic
<li>Intergenic
<li>Exonic
</ul>
<p class="text-justify">
This menu allows the selection of RNA editing sites residing in specific genic regions such as: untranslated regions (UTR) or intronic regions or coding/non-coding exons or intergenic regions.
Reported classification has been carried out by <a href="http://annovar.openbioinformatics.org/en/latest/" target="blank">ANNOVAR</a>.
</p>
</td>
</tr>
<tr>
<td width="30%"><img src="img/sc5.png" width="100%" height="100%"/></td>
<td style='text-align:center' width="20%"><b>AA Change</b></td>
<td>
<ul>
<li>Synonymous
<li>Nonsynonymous
<li>Stop Loss
<li>Unknown
</ul>
<p class="text-justify">
This menu allows the selection of RNA editing sites residing in protein coding regions and affecting codon integrity.
Reported classification has been carried out by <a href="http://annovar.openbioinformatics.org/en/latest/" target="blank">ANNOVAR</a>. 
</p>
</td>
</tr>
<tr>
<td width="30%"><img src="img/sc6.png" width="100%" height="100%"/></td>
<td style='text-align:center' width="20%"><b>Tissue</b></td>
<td>
<ul>
<li>Adipose Tissue
<li>Adrenal Gland
<li>Blood
<li>...
<li>Thyroid
</ul>
<p class="text-justify">
This menu allows the selection of RNA editing sites residing in specific human tissues. More than one tissue can be selected per each search.
Tissue names are according to <a href="http://www.gtexportal.org/home/tissueSummaryPage" target="blank">GTEx</a>. 
</p>
</td>
</tr>
<tr>
<td width="30%"><img src="img/sc7.png" width="100%" height="100%"/></td>
<td style='text-align:center' width="20%"><b>Body Site</b></td>
<td>
<ul>
<li>Brain - Hypothalamus
<li>Brain - Substantia nigra
<li>...
<li>Whole Blood 
</ul>
<p class="text-justify">
This menu allows the selection of RNA editing sites residing in specific human body sites. More than one body site can be selected per each search.
Body site names are according to <a href="http://www.gtexportal.org/home/tissueSummaryPage" target="blank">GTEx</a>. 
</p>
</td>
</tr>
</table>
<p class="text-justify">
A search example can be performed clicking on "Example" button.
All searches, instead, are activated by clicking the "Search" button. The search form can also be reset by clicking the "Clean" button.
</p>
<p class="text-justify">
Once a search has been performed, results will be displayed in a table including the following columns:
</p>
<br>
<div align="center"><img src="img/sc8.png" width="100%" height="100%"/></div>
<br><br>

<table class="table table-bordered">
<tr>
<th width="20%">Column Name</th>
<th>Meaning</th>
</tr>
<tr>
<td><b>Chr</b></td>
<td><p class="text-justify">Chromosome Name</p></td>
</tr>
<tr>
<td><b>Position</b></td>
<td><p class="text-justify">Chromosome Coordinate</p></td>
</tr>
<tr>
<td><b>Ref</b></td>
<td><p class="text-justify">Reference Nucleodite</p></td>
</tr>
<tr>
<td><b>Ed</b></td>
<td><p class="text-justify">Edited Nucleotide</p></td>
</tr>
<tr>
<td><b>Strand</b></td>
<td><p class="text-justify">Strand (+ or -)</p></td>
</tr>
<tr>
<td><b>dbSNP</b></td>
<td><p class="text-justify">a colored flag indicating the presence of a SNP in dbSNP. Only SNPs classified as "genomic" are taken into account. A green flag indicates a match with dnSNP and provides also an external link to NCBI</p>
<div align="center"><img src="img/sc10.png" width="100%" height="100%"/></div>
</td>
<tr>
<td><b>Location</b></td>
<td><p class="text-justify">Location of RNA Editing in repetitive or non-repetitive regions.</p></td>
</tr>
<tr>
<td><b>Repeats</b></td>
<td><p class="text-justify">Class and family of repeat including the RNA editing position.</p></td>
</tr>
<tr>
<td><b>Gene</b></td>
<td><p class="text-justify">Gene Symbol</p></td>
</tr>
<tr>
<td><b>Region</b></td>
<td><p class="text-justify">Genic Region according to ANNOVAR</p></td>
</tr>
<tr>
<td><b>EditedIn</b></td>
<td><p class="text-justify">The number of Samples in which the specific position appears to be edited. It is showed by a progression bar.</p>
<div align="center"><img src="img/sc11.png" width="100%" height="100%"/></div>
</td>
</tr>
<tr>
<td><b>ExFun</b></td>
<td><p class="text-justify">Exonic function limited to synonymous and non-synonymous positions. A colored flag is used to indicate if a site is synonymous (green) or non-synonymous (red). Click on to open a pop-up with details.</p>
<div align="center"><img src="img/sc9.png" width="100%" height="100%"/></div>
</td>
</tr>
<tr>
<td><b>Phast</b></td>
<td><p class="text-justify">PhastCons conservation scores calculated for multiple alignments
of 45 vertebrate genomes to the human genome. It ranges from 0 (no conservation) to 1000 (max conservation). Values derive from UCSC phastCons46way table.</p></td>
</tr>
<tr>
<td><b>KnownIn</b></td>
<td><p class="text-justify">A colored flag indicating the presence of a site in other available database (A: ATLAS, R: RADAR, D: DARNED). Click on R or D to open an external link to RADAR or DARNED databases, respectively.</p></td>
</tr>
</table>

<p class="text-justify">
For each position, REDIportal provides additional info by clicking on blue arrow in the first column.
This will cause the opening of four tabs. The first tab named "Heat-Map" displays an RNA Editing heat-map in which mean editing level per body site is reported.
Mouse over each body site to open a tooltip showing the average editing level.  
</p>
<br>
<div align="center"><img src="img/sc12.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
The second tab named "Box Plot" displays RNA Editing levels per each body site by means of box plots. Relevant values are available by mousing over each box plot.
</p>
<br>
<div align="center"><img src="img/sc13.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
The third tab named "Alternative Annotations" displays a table with gene/transcript annotations from RefSeq database and UCSC KnownGene table. 
</p>
<br>
<div align="center"><img src="img/sc14.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
The last tab named "Editing Details" displays the number of samples, tissues and body sites in which the position appears to be edited.
Clicking on "View Editing Details" button will cause the opening of a new windows with a table including editing levels per each experiment.
</p>
<br>
<div align="center"><img src="img/sc15.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
The "View Editing Details" button enables the opening of a new windows including relevant editing info described in the table below.
</p>
<br>
<div align="center"><img src="img/sc16.png" width="100%" height="100%"/></div>
<br><br>


<table class="table table-bordered">
<tr>
<th width="20%">Column Name</th>
<th>Meaning</th>
</tr>
<tr>
<td><b>RNAseq Run</b></td>
<td><p class="text-justify">RNAseq Run accession number according to SRA database.</p></td>
</tr>
<tr>
<td><b>WGS Run</b></td>
<td><p class="text-justify">Whole Genome Sequencing Run accession number according to SRA database.</p></td>
</tr>
<tr>
<td><b>Tissue</b></td>
<td><p class="text-justify">Tissue Name according to GTEx project.</p></td>
</tr>
<tr>
<td><b>BodySite</b></td>
<td><p class="text-justify">Body Site Name according to GTEx project.</p></td>
</tr>
<tr>
<td><b>n.As</b></td>
<td><p class="text-justify">Number of RNAseq reads supporting Adenosine</p></td>
</tr>
<tr>
<td><b>n.Gs</b></td>
<td><p class="text-justify">Number of RNAseq reads supporting Guanosine</p></td>
</tr>
<tr>
<td><b>EditingFreq</b></td>
<td><p class="text-justify">RNA Editing Frequecy</p></td>
</tr>
<tr>
<td><b>gCoverage</b></td>
<td><p class="text-justify">Number of supporting genomic reads</p></td>
</tr>
<tr>
<td><b>gFreq</b></td>
<td><p class="text-justify">Max Frequency of AG change at genomic level.</p></td>
</tr>
</table>

<!--p class="text-justify">
Individual run or tissue or bosy sites can be selected by using the "Select" button below each column.
Numerical columns can be sorted by clicking on each column title.
</p>
<br>
<div align="center"><img src="img/sc17.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
Result table can be downloaded and exported in Excel or PDF format for further analyses. 
</p>
<br>
<div align="center"><img src="img/sc18.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
Result table can also be filtered by clicking on the "Filter Editing Levels" button. 
This will cause the opening of a pop-up in which the user can insert numeric values to filter RNA editing levels as well as reads supporting adenosines or guanosines. 
</p>
<br>
<div align="center"><img src="img/sc19.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
Specific columns of result table can be hided by clicking on "Column visibility" button. 
</p>
<br>
<div align="center"><img src="img/sc20.png" width="100%" height="100%"/></div>
<br><br-->

<p class="text-justify">
Users can increase the number of visible rows by using the "Show" button. 
</p>
<br>
<div align="center"><img src="img/sc21.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
Also in the main result table, specific columns can be hided by clicking on "Column visibility" button.
</p>
<br>
<div align="center"><img src="img/sc22.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
Search results can be downloaded using the "Download" button. This will cause the opening of a pop-up in which users can select columns to download.
</p>
<br>
<div align="center"><img src="img/sc23.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
Columns of each result table can be exchanged or moved in order to customize the aspect and column order. 
</p>
<br>
<div align="center"><img src="img/sc24.png" width="100%" height="100%"/></div>
<br><br>

<p class="text-justify">
Columns with gray arrows are sortable in ascending or descending order. 
</p>
<br>
<div align="center"><img src="img/sc25.png" width="100%" height="100%"/></div>
<br><br>


<!--hr>
<h4 id="S4">Search RNAseq samples in REDIportal</h4>
<p class="text-justify">
REDIportal allows also RNA editing searches at sample level. Users can browse RNA editing statistics detected in each RNAseq experiment by selecting the “Search Sample” page from the main menu and providing specific options. Indeed, users can interrogate the database introducing the name of a sample in the “Sample name” form using the run accession number provided by SRA (or dbGAP or GTEx) (for example SRR1069188) or can select samples by source (GTEx) or status (normal, tumor) or data type (bulk tissue or single cell) or tissue or body site. In addition, samples can be selected according to the expression of ADAR genes or Alu Editing index values.
</p>
<br>
<div align="center"><img src="img/sm1.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
Once a search has been performed, results will be displayed in a table including the following columns:
</p>
<br>
<div align="center"><img src="img/sm2.png" width="100%" height="100%"/></div>
<br><br>

<table class="table table-bordered">
<tr>
<th width="20%">Column Name</th><th>Meaning</th></tr>
<tr>
<td><b>Sample</b></td>
<td><p class="text-justify">Sample Name (RNAseq accession number)</p></td>
</tr>
<tr>
<td><b>WGS/WES</b></td>
<td><p class="text-justify">WGS/WES Name (DNAseq accession number) from the same individual, if available</p></td>
</tr>
<tr>
<td><b>Source</b></td>
<td><p class="text-justify">Project source name</p></td>
</tr>
<tr>
<td><b>Organism</b></td>
<td><p class="text-justify">Organism name</p></td>
</tr>
<tr>
<td><b>Events</b></td>
<td><p class="text-justify">Number of RNA editing events detected in the sample</p></td>
</tr>
<td><b>Hyper</b></td>
<td><p class="text-justify">Number of hyper-edited events detected in the sample</p></td>
</tr>
<tr>
<td><b>Body Site</b></td>
<td><p class="text-justify">Name of the body site</p></td>
</tr>
<tr>
<td><b>Status</b></td>
<td><p class="text-justify">Disease Status</p></td>
</tr>
<tr>
<td><b>Type</b></td>
<td><p class="text-justify">Tissue type: bulk or single cell</p></td>
</tr>
<tr>
<td><b>AEI</b></td>
<td><p class="text-justify">Alu Editing Index</p>
</td>
</tr>
<tr>
<td><b>REI</b></td>
<td><p class="text-justify">Recoding Editing Index</p>
</td>
</tr>
<tr>
<td><b>ADAR</b></td>
<td><p class="text-justify">Expression of ADAR gene (in TPM)</p></td>
</tr>
<tr>
<td><b>ADARB1</b></td>
<td><p class="text-justify">Expression of ADARB1 gene (in TPM)</p></td>
</tr>
<tr>
<td><b>ADARB2</b></td>
<td><p class="text-justify">Expression of ADARB2 gene (in TPM)</p></td>
</tr>
</table>
<p class="text-justify">
For each sample, REDIportal provides additional info by clicking on blue arrow in the first column.
This will cause the opening of five tabs.
</p>
<br>
<table class="table table-bordered">
<tr><th width="20%">Tab name</th><th>Content</th></tr>

<tr><th width="20%">Genomics Facts</th>
<th><p class="text-justify">Main statistics about the genomics location of detected RNA editing events</p>
<div align="center"><img src="img/sm3.png" width="100%" height="100%"/></div>
</th></tr>

<tr><th width="20%">Base Distribution</th>
<th><p class="text-justify">It shows the distribution of detected variants by our HPC REDItools pipeline</p>
<div align="center"><img src="img/sm4.png" width="100%" height="100%"/></div>
</th></tr>

<tr><th width="20%">RNA Editing Indices</th>
<th><p class="text-justify">Box plots of AEI and REI indices for the specific body site. Details of recoding events per sample are available by clicking on the "REI details"</p>
<div align="center"><img src="img/sm5.png" width="100%" height="100%"/></div>
</th></tr>

<tr><th width="20%">RNA Editing Levels</th>
<th><p class="text-justify">Distribution of RNA editing levels from detected sites</p>
<div align="center"><img src="img/sm6.png" width="100%" height="100%"/></div>
</th></tr>

<tr><th width="20%">Transcriptome Coverage</th>
<th><p class="text-justify">Fraction of edited genes over the entire annotation. Details for each edited gene per sample are available by clicking on the "Gene details"</p>
<div align="center"><img src="img/sm7.png" width="100%" height="100%"/></div>
</th></tr>
</table>
<hr>
<h4 id="S5">Browse RNA Editing sites per gene locus</h4>
<p class="text-justify">
RNA editing events stored in REDIportal are also visible in their genic context through our novel Gene View functionality.
Users can explore known events per gene by selecting the "Gene View" page from the Search Menu and provide the name of the favourite gene (according to available organisms and genome assemblies).
</p>
<div align="center"><img src="img/sm8.png" width="50%" height="50%"/></div>
<p class="text-justify">
Once a gene has been selected, the user will be able to see the structure of the gene locus organised in transcripts and a panel containing all known editing events for the specific locus. Users can also zoom on specific gene locations.
</p>
<div align="center"><img src="img/sm9.png" width="100%" height="100%"/></div>
<hr>
<br>
<h4 id="S2">Browse RNA Editing sites in JBrowse</h4>
<p class="text-justify">
All RNA editing events stored in REDIportal are visible in their genomic context through <a href="http://jbrowse.org/" target="blank">JBrowse</a>, a fast genome browser based on JavaScript and HTML5. It is embedded in REDIportal by default, allowing the browsing of basic tracks such as individual RNA editing sites, SNPs, RefSeq gene annotations, Alu elements and LINEs.
</p>
<br>
<div align="center"><img src="img/img1.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
Genomic intervals can be inspected entering chromosome coordinates in the JBrowse search box (inside the red ovale below) using the format Chr:start..end or Chr:start-end (for example chr4:158101247..158308846). Commas can be used as thousands separators (like in UCSC) in the start and stop nucleotide position numbers but they are not required.
</p>
<br>
<div align="center"><img src="img/img2.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
Alternatively, the JBrowse search box accepts gene symbols and allows the autocomplete function to easily suggest gene names during the typing (an example is shown in the red circle below).
</p>
<br>
<div align="center"><img src="img/img3.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
If the gene symbol is present in multiple JBrowse tracks implemented in REDIportal, a dialog window will be opened allowing the selection of the correct track.
In the example below, GRIA1 gene is in both Gencode Basic V19 track and RefSeq track. The dialog window will allow the selection of needed track using "Go" buttons.  
</p>
<br>
<div align="center"><img src="img/img4.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
Navigation buttons are located to the left of the search box in the consolidated header region. 
The arrow buttons move the view about the distance of one screen left or right.
The larger zoom buttons zoom in or out about twice as far as the smaller buttons.
</p>
<br>
<div align="center"><img src="img/img5.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
In addition to these buttons, JBrowse supports click and drag selection of regions in both the chromosome-level and detail-level position bars.
</p>
<br>
<div align="center"><img src="img/img6.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
Each JBrowse data track has a context-specific menu, hidden by default. The down arrow on the title bar of the track allows the visualization of the following options:
</p>
<ul>
<li><p class="text-justify"> <b>About this track</b>: provides some additional information about a particular track such as the track type, category and legend.</p>
<li><p class="text-justify"> <b>Pin to top</b>: causes that track to always be displayed directly beneath the header area at the top of the browser window. </p>
<li><p class="text-justify"> <b>Edit Config</b>: allows the user to directly edit the configuration script for a particular track, even though it is not recommended for most users.</p>
<li><p class="text-justify"> <b>Delete track</b>: turns off individual tracks. </p>
<li><p class="text-justify"> <b>Save Track Data</b>: allows the viewing and saving of track data in gff3, bed or sequin format. Reference sequences can be exported in fasta format.</p>
<li><p class="text-justify"> <b>Display mode</b>: enables three display modes 1) "Normal" view, 2) "Compact" with reduced height of each object in the track and 3) "Collapse" (default for RNA editing and SNP tracks) moves all objects to a single line on the track.</p>
<li><p class="text-justify"> <b>Show labels</b>: displays labels when the view is zoomed in sufficiently. "Show labels" box is turned off by default in RNA editing and SNP tracks.</p>
</ul>
<br>
<div align="center"><img src="img/img7.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
Details of each track can be explored by clicking on each annotation, since it will cause the opening of a specific pop-up window. 
In addition, left-clicks on features will open an embedded popup window showing further options:
</p>
<ul>
<li><p class="text-justify"> <b>Zoom</b>: allows the zoom on the specific feature;</p>
<li><p class="text-justify"> <b>Highlight</b>: enables the highlighting of a feature (this behaviour can be disabled clicking the highlight button on navigation bar);</p>
<li><p class="text-justify"> <b>Link</b> to a specific web page to recover additional info. In case of gene annotations, they are linked to GeneCards database. SNPs are connected to NCBI dbSNP while RNA editing sites are linked to REDIportal details including info about edited tissues, body sites and samples with correlated frequencies values.</p>
<li><p class="text-justify"> <b>View details</b>: features with no specific links have a “View details“ option to open a pop-up window as above.</p>
</ul>
<br>
<div align="center"><img src="img/img8.png" width="50%" height="50%"/></div>
<br><br>
<p class="text-justify">
JBrowse embedded in REDIportal includes also further tracks, available as a list in the left side bar, visible only in the full-screen modality. Such full-screen view can be enabled clicking the “Full-screen view” link in the upper right corner.
</p>
<br>
<div align="center"><img src="img/img9.png" width="50%" height="50%"/></div>
<br><br><br>
<div align="center"><img src="img/img10.png" width="50%" height="50%"/></div>
<br><br>
<hr>
<h4 id="S3">Search Cell Lines in the CLAIRE (Cell Line A-to-I Rna Editing) database</h4>
<p class="text-justify">
Users can now perform searches across the CLAIRE database (Cell Line A-to-I Rna Editing) enabling the identification of cell lines
suitable for investigating specific RNA editing sites.
Cell lines can be retrieved according to the AEI (Alu Editing Index) and the expression levels of ADAR and ADARB1 genes.
</p>
<br>
<div align="center"><img src="img/img11.png" width="50%" height="50%"/></div>
<br><br>

<p class="text-justify">
Up to 5 target sites (in hg19 coordinates) can be selected.
</p>
<br>
<div align="center"><img src="img/img12.png" width="50%" height="50%"/></div>
<br><br>


<p class="text-justify">
For each site, specific search parameters can be tuned such as the expression of the target gene, the coverage depth (number of reads supporting the site) and the RNA editing level.
</p>
<br>
<div align="center"><img src="img/img13.png" width="50%" height="50%"/></div>
<br><br>

<p class="text-justify">
Once a search has been performed, results will be displayed in a table including the following columns:
</p>
<br>
<div align="center"><img src="img/img14.png" width="100%" height="100%"/></div>
<br><br>
<table class="table table-bordered">
<tr>
<th width="20%">Column Name</th>
<th>Meaning</th>
</tr>
<tr>
<td><b>Cell Line</b></td>
<td><p class="text-justify">Cell Line name</p></td>
</tr>
<tr>
<td><b>Sample</b></td>
<td><p class="text-justify">Sample name</p></td>
</tr>
<tr>
<td><b>AEI</b></td>
<td><p class="text-justify">Alu Editing Index</p></td>
</tr>
<tr>
<td><b>ADAR expr</b></td>
<td><p class="text-justify">ADAR expression (TPM)</p></td>
</tr>
<tr>
<td><b>ADARB1 expr</b></td>
<td><p class="text-justify">ADARB1 expression (TPM)</p></td>
</tr>
<tr>
<td><b>Tissue</b></td>
<td><p class="text-justify">Human tissue of origin</p></td>
</tr>
</table>

<p class="text-justify">
Cell lines or samples or tissues can be selected by using the "Select" button below each column and all columns can be sorted (in ascending or descending order) by clicking on each column title.
Results can be downloaded and exported in Excel or PDF format for further analyses.<br>
For each cell line, additional info are provided in tabular format by clicking on blue arrow in the first column. They include the expression of the target gene (in TPM), the coverage depth (number of reads) of the specific target site and its editing level.
</p>
<br>
<div align="center"><img src="img/img15.png" width="100%" height="100%"/></div>
<br><br>
<br><br>
    </div>
    <div class="col-sm-2 sidenav"-->

<h3>Practical part</h3>
<p>After searching for GRIA2 and FLNA genes use the "Download" button embedded with the Results table and select the right columns to obtain a table separated file compatible with the REDItoolKnown.py script<br> See below for further details on this file format.<p>
The main steps described during the practice are reported below and can be easily copy/pasted in your terminal.<br>
Note. Assuming you're traineeX, please change X according to your workspace.<br>
IMPORTANT! REDItoolKnown.py outTable (eg. outTable_892028847) contains 9digit random number, so it usually varies among users and different script launches on the same machine.<br>
Due to multiple available versions of the core module pysam, it is possibile that some commands will return you a pysam error.
In those cases just type:
    
<pre>$ conda activate rnaediting2</pre>

*rnaediting2 environment contains pysam=0.15.2
Type again the command that returned errors and revert to your main environment with:

<pre>$ conda activate rnaediting</pre>

*rnaediting envirnment contains pysam==0.7.7     
 
<p>1) Log into your area and create two separate folders for strand-oriented and unstranded RNAseq data (eg. RNAseq_strnd)</p>
<pre> 
$ mkdir RNAseq_strnd
$ cd ..
$ mkdir RNAseq_unstrnd
$ cd ..
</pre>
<p>2) According to each folder copy the RNAseq data from Editing_knwn folder</p>
<pre>
$ cd RNAseq_unstrnd
$ cp /usr/share/course_data/rnaediting/Editing_knwn/Unstrnd/*.bam* .
$ cd ..
$ cd RNAseq_strnd/
$ cp /usr/share/course_data/rnaediting/Editing_knwn/Strndd/*.bam* .
$ cd ..
</pre>
<b>Note.</b>If unable to upload them from your computer to your home folder, you can recover the tab files from /usr/share/course_data/rnaediting/Editing_knwn/, by entering your workspace and giving the command:
<pre>$ cp /usr/share/course_data/rnaediting/Editing_knwn/*.gz* .</pre>

<p>3) Launch REDItoolKnown.py on RNAseq unstranded data and check them ONLY for those known positions extracted previously from REDIportal (e.g GRIA2)</p>
<pre>
$ cd RNAseq_unstrnd
$ REDItoolKnown.py -i Cerebellum_unstrnd.bam -f  /usr/share/course_data/rnaediting/hg19ref/GRCh37.primary_assembly.genome.fa -l ../GRIA2nrptAtlasTable.txt.gz  
$ REDItoolKnown.py -i Lung_unstrnd.bam -f  /usr/share/course_data/rnaediting/hg19ref/GRCh37.primary_assembly.genome.fa -l ../FLNAAtlasTable.txt.gz
$ cd ..
</pre>

<p>4) Launch REDItoolKnown.py on RNAseq strand-oriented data and check them ONLY for those known positions extracted previously from REDIportal (e.g GRIA2)</p>
<pre>
$ cd RNAseq_strnd
$ REDItoolKnown.py -i SRR-6H_HD.bam -f  /usr/share/course_data/rnaediting/hg19ref/GRCh37.primary_assembly.genome.fa -l ../GRIA2nrptAtlasTable.txt.gz -s2 -S
$ REDItoolKnown.py -i SRR-6H_HD.bam -f  /usr/share/course_data/rnaediting/hg19ref/GRCh37.primary_assembly.genome.fa -l ../AGBL4AtlasTable.txt.gz -s2 -S
$ cd ..
</pre>

<h2>TAB<a class="headerlink" href="#tab" title="Permalink to this headline"></a></h2>
<p>TAB files are simple textual files with at least three tabulated columns including:</p>
<ul class="simple">
<li>genomic region (generally the chromosome name according to the reference genome)</li>
<li>coordinate of the position (1-based)</li>
<li>strand (+ or -). You can also indicate strand by 0 (strand -), 1 (strand +) or 2 (+ and - or unknown)</li>
</ul>
<table border="0" class="docutils">
<thead valign="bottom">
<tr class="row-odd"><th class="head">genomic region</th>
<th class="head">coordinate</th>
<th class="head">strand</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td>chr21</td>
<td>10205589</td>
<td>-</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>10205629</td>
<td>-</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15411496</td>
<td>+</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15412990</td>
<td>+</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15414553</td>
<td>+</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15415901</td>
<td>+</td>
</tr>
<tr class="row-even"><td>chr21</td>
<td>15417667</td>
<td>+</td>
</tr>
<tr class="row-odd"><td>chr21</td>
<td>15423330</td>
<td>+</td>
</tr>
</tbody>
</table>
<p>TAB files must be coordinate sorted. In unix/linux environment they can be sorted by the sort command:</p>
<div class="highlight-python"><pre>sort -k1,1 -k2,2n mytable.txt &gt; mytable.sorted.txt</pre>
</div>
    
<div class="section" id="reditoolknown-py">
<h2>REDItoolKnown.py<a class="headerlink" href="#reditoolknown-py" title="Permalink to this headline"></a></h2>
<p>REDItoolKnown.py has been developed to explore the RNA editing potential of RNA-Seq data
sets using known editing events. Such events can be downloaded from REDIportal database or
generated from supplementary materials of a variety of publications. Known RNA editing events
have to be stored in TAB files (see above for details).</p>
<dl class="docutils">
<dt>Options:</dt>
<dd><table class="first last docutils option-list" frame="void" rules="none">
<col class="option" />
<col class="description" />
<tbody valign="top">
<tr><td class="option-group">
<kbd><span class="option">-i</span></kbd></td>
<td>BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-I</span></kbd></td>
<td>Sort input BAM file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-f</span></kbd></td>
<td>Reference in fasta file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-l</span></kbd></td>
<td>List of known RNA editing events</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-C</span></kbd></td>
<td>Base interval to explore [100000]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-k</span></kbd></td>
<td>List of chromosomes to skip separated by comma or file</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-t</span></kbd></td>
<td>Number of threads [1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-o</span></kbd></td>
<td>Output folder [rediFolder_XXXX] in which all results will be stored. XXXX is a random number generated at each run.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-F</span></kbd></td>
<td>Internal folder name [null] is the main folder containing output tables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-c</span></kbd></td>
<td>Min. read coverage [10]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-q</span></kbd></td>
<td>Minimum quality score [25]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-m</span></kbd></td>
<td>Minimum mapping quality score [25]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-O</span></kbd></td>
<td>Minimum homoplymeric length [5]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-s</span></kbd></td>
<td>Infer strand (for strand oriented reads) [1]. It indicates which read is in line with RNA. Available values are: 1:read1 as RNA,read2 not as RNA;
2:read1 not as RNA,read2 as RNA; 12:read1 as RNA,read2 as RNA; 0:read1 not as RNA,read2 not as RNA.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-g</span></kbd></td>
<td>Strand inference type 1:maxValue 2:useConfidence [1]; maxValue: the most prominent strand count will be used; useConfidence: strand is assigned if over a prefixed frequency confidence (-x option)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-x</span></kbd></td>
<td>Strand confidence [0.70]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-S</span></kbd></td>
<td>Strand correction. Once the strand has been inferred, only bases according to this strand will be selected.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-G</span></kbd></td>
<td>Infer strand by GFF annotation (must be sorted, otherwise use -X). Sorting requires grep and sort unix executables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-X</span></kbd></td>
<td>Sort annotation files. It requires grep and sort unix executables.</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-K</span></kbd></td>
<td>File with positions to exclude (chromosome_name coordinate)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-e</span></kbd></td>
<td>Exclude multi hits</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-d</span></kbd></td>
<td>Exclude duplicates</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-p</span></kbd></td>
<td>Use paired concardant reads only</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-u</span></kbd></td>
<td>Consider mapping quality</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-T</span></kbd></td>
<td>Trim x bases up and y bases down per read [0-0]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-B</span></kbd></td>
<td>Blat folder for correction</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-U</span></kbd></td>
<td>Remove substitutions in homopolymeric regions</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-v</span></kbd></td>
<td>Minimum number of reads supporting the variation [3]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-n</span></kbd></td>
<td>Minimum editing frequency [0.1]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-E</span></kbd></td>
<td>Exclude positions with multiple changes</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-P</span></kbd></td>
<td>File containing splice sites annotations (SpliceSite file format see above for details)</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-r</span></kbd></td>
<td>Num. of bases near splice sites to explore [4]</td></tr>
<tr><td class="option-group">
<kbd><span class="option">-h</span></kbd></td>
<td>Print the help</td></tr>
</tbody>
</table>
</dd>
</dl>
<p>Example:</p>
<div class="highlight-python"><pre>REDItoolKnown.py -i rnaseq.bam -f reference.fa -l knownEditingSites.tab</pre>
</div>
</div>

