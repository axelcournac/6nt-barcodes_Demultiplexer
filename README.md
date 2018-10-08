## The program

Demultiplexes reads produces by Illumina's NextSeq500 and barcoded with Bioo NEXTflex 6nt barcodes.

Input: 1 or 2 fastq files  
Output: a folder with a fastq file for each barcode


## Compile:
Demultiplexv2.cpp must be compiled with a compiler supporting the C++11 standard. For example:

```g++ -std=c++11 Demultiplexv2.cpp -o <executable file name>```

## Usage:

```./Demultiplex <fastq file (end1)> [<fastq file (end2)>]```

1 (for single-end reads) or 2 (for pair-end reads) files can be passed as argument.

The executable file must be in the same folder as the file *barcodes.txt*.

*barcodes.txt* is a file with 2 comma-separated columns. 
The 1st column corresponds to the 4 characters name of the barcodes.
The 2nd column corresponds to the 6 nucleotides barcode.

A folder `Demultiplexed_Reads` is created, where the demultiplexed fastq files are saved.
