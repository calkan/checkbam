# verifybam

verifybam is a BAM file integrity verification tool.

Read mapping is verified by applying CIGAR and MD operations to both original read and corresponding sequence in reference genome. Also reads are cross-checked against original FastQ files to see if any reads are missing. It is crucial to use the same reference genome and FastQ files that are used in alignment process.

## Download

If htslib is not installed on your system, it is important to clone the project with htslib for building.
To understand whether htslib is install, you can run ldconfig command as follows:

```ldconfig -p | grep libhts```

If htslib is not installed, there will be no output. Use following clone command to also download htslib.

```git clone --recursive https://github.com/calkan/verifybam```

Otherwise, simply clone the project without htslib.

```git clone https://github.com/calkan/verifybam```

## Build

Makefile is already configured whether you have htslib or not. If it is not already installed on your system, make will try to build the htslib from local repository. This is going to *fail* if you have not downloaded the htslib submodule as mentioned.

```make all```

This command should create the necessary executable 'verifybam'.

```make install```

To install verifybam.

## Usage

Available options can be obtained by running verifybam without any parameters.

```
VERIFYBAM: BAM validity checking tool.
Version 0.2-alpha
	Last update: August 01, 2016, build date: Fri Aug 26 08:19:31 EEST 2016

	--input [BAM file]   : Input file in sorted and indexed BAM format.
	--ref   [reference]  : Reference genome in FASTA format.
	--fastq [Fastq file] : A fastq file that contains original reads of input BAM file. Can be given multiple times.
	--version            : Print version and exit.
	--help               : Print this help screen and exit.
```
