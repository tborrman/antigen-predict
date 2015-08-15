# antigen-predict
Code for T cell receptor cross-reactive antigen prediction

### Data processing

Sequencing data downloaded from Sequence Read Archive (SRA) under project code [SRP040021](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP040021)

- Data converted to fastq files using `fastq-dump` with `--split-3` option for mate pairs
- Extract reads from individual selection rounds with `[PDB_id]_process_step1.py` 

### Validate recovered peptide sequences
- `[PDB]_seq_check.py`
- `[PDB]_seq_check_heatmap.R`



