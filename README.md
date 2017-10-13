# Cottrell_PTRE-seq_scripts
Scripts used in Cottrell et al, 2017 manuscript<br /><br />
Steps and Scripts
1. Library design 
<br /> Library was designed manually. This GitHub repo includes the excel file: Library_design.xlsx
<br /> Grouping of regulatory elements into specific categories: Grouped_Barcode_identities.txt

2. Raw sequencing data to Barcode counts<br />
All replicates were indexed with a different P1 index and pooled and sequenced. FastQ files were first demultiplexed, and then counts were obtained per barcode. <br />Script: KC_barcode_counts.sh 

3. Preliminary analysis: Counts to RNA expression and TE efficiency
<br />Script: Analysis_RNA_TE.R<br />
Files: 
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1) RNA_plasmid.txt (RNA Expression)<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2) PARNA_RNA.txt (TE efficiency)<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3) Analysis_cell_type.R (Counts to RNA expression for different cell types)<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4) Analysis_40S.R (Counts to 40S association)

4. Linear Modeling to RNA expression and TE efficiency 
<br />Scripts: <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1) Modeling_linear.R<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2) Exploratory analysis: modeling_FOLD.R<br />

5. Scripts for specific data analysis
<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1) Analysis of Let7 targets: Let-7_Analysis.R<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2) Analysis of combined effects of Let7 and Pumilio: Pum_Let-7_Analysis.R<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3) Natural target analysis scripts: Natural_Targets_Analysis.R<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4) Analysis of AU-rich element containing reporters: HuR_Analysis<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;5) Analysis of reporters with combinations of Let-7 binding sites and AREs: HuR_Let_7_Analysis.R<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6) Analysis of reporters with combinations of PREs and AREs: HuR_Pum_Analysis.R<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;7) Analysis of let-7 binding site variants: Seed_PC_Analysis.R

