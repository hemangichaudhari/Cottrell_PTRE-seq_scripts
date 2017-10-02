# Cottrell_PTRE-seq_scripts
Scripts used in Cottrell et al, 2017 manuscript<br /><br />
Steps and Scripts
1. Library design 
<br /> Library was designed manually. This GitHub repo includes the excel file
2. Raw sequencing data to Barcode counts<br />
All replicates were indexed with a different P1 index and pooled and sequenced. FastQ files were first demultiplexed, and then counts were obtained per barcode. <br />Script: KC_barcode_counts.sh
3. Preliminary analysis: Counts to RNA expression and TE efficiency
<br />Script: Analysis_RNA_TE.R
4. Linear Modeling to RNA expression and TE efficiency 
<br />Scripts: <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1) modeling_linear.R <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2) Exploratory analysis : modeling_FOLD.R
5. Example scripts:<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1) Analysis of Let7 targets : Let7_analysis.R<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2) Analysis of combined effects of Let7 and Pumilio : Pum_Let7_analysis.R
6. Natural target analyis scripts : Natural_Target_Analysis.R

