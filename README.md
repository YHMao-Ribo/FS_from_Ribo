This is a Perl Script used for frameshifting analysis based on Ribo-seq.
Two text files are required for this running.
1) A gzip compressed fasta file that contains all cDNA sequences. The format of the header for the fasta file can be found in the demo file.
2) A gzip compressed text file that contains read count of Ribo-seq on individual positions of mRNAs. Each row represents a dataset on an individual mRNA. All data are seperated by ",".

No additional package required.

By using command line: Perl codon_inframe.pl cDNA reads
where cDNA is the gziped fasta file, and reads is the gziped read count file.
The script will calculated in-frame rate of all 61*61 individual pairs, which was used to indicate the frequency of out-of-frame translation.
