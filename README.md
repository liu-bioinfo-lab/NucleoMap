# NucleoMap

NucleoMap is a computational approach to identify nucleosomes from high-resolution chromatin contact maps.

## Required environment
Python 3.7.3\
numpy 1.16.2\
scipy 1.4.1\
pandas 0.24.2\
itertools 6.0.0\
bedtools v2.28.0

## Usage
### Data processing
Before calling nucleosomes, we need to prepare 1) pair files of high-resolution contact map, and 2) reference genome in fasta format.\
The required pair file format is a tab separated file with at least 7 columns in the following format:
```
Column1: contact ID, STR
Column2: chromosome 1, STR
Column3: position 1, INT
Column4: chromosome 2, STR
Column5: position 2, INT
Column6: strand 1, written in +/-
Column7: strand 2, written in +/-
```
Due to high computational burden, it is impossible to feed the entire contact map into the memory.
Therefore, the pair file needs to be separated according to chromosomes. Inter-chromosome contacts are discarded.\
After that, nucleosome centers and the adjusted alignments are calculated from the separated pair files using 
```
python infer_center.py <pair file> <out dir>
```
The output of this step is a bed file ``reads.bed`` containing adjusted alignments and a numpy array ``read_centers.npy`` containing inferred read centers.

Next, DNA sequence of the adjusted alignments are extracted using 
```
bedtools getfasta -fi ref_genome.fa -bed reads.bed > output.fa
```

At the end, the binding preference is calculated from the extracted DNA sequence.

### Nucleosome calling
Nucleosomes are called using
```
python nucleomap.py <read centers> <out dir> <PWM> <reference genome> 
```

### Nucleosome contact map generating
To be updated
