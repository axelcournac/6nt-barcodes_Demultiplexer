# The program

Demultiplex reads produces by Illumina's NextSeq500 and barcoded with Bioo NEXTflex 6nt barcodes.

# Compile:
Demultiplexv2.cpp must be compied with compiler supporting the C++11 standard. For example:

```g++ --std=c++11 Demultiplexv2.cpp -o <executable file name>```

# Usage:

```./Demultiplex <fastq file (end1)> <fastq file (end2)>```

The executable file must be in the same folder as the file barcodes.txt.

barcodes.txt is a file with 2 comma-separated columns. 
The 1st column corresponds to the 4 character name of the barcodes.
The 2nd columns corresponds to the 6 nucleatides barcode.

A folder `Demultiplexed_Reads` is created, where the demultiplexed fastq files are saved.
