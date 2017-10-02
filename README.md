# Cottrell_PTRE-seq_scripts
Scripts used in Cottrell et al, 2017 manuscript<br /><br />
Steps and Scripts
1. Library design
2. Raw sequencing data to Barcode counts<br />
All replicates were indexed with a different P1 index and pooled and sequenced. FastQ files were first demultiplexed, and then counts were obtained per barcode. <br />Script: KC_barcode_counts.sh
3. Preliminary analysis: Counts to RNA expression and TE efficiency
4. Linear Modeling to RNA expression and TE efficiency <br />Script: modeling_FOLD.R

