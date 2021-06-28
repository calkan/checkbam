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

### Troubleshooting

#### zlib.h no such file or directory

Some linux distributions are not shipped with zlib header files. For debian systems such as Ubuntu, Linux Mint, you can download ```libz-dev``` package from package manager.
```
sudo apt-get install libz-dev
```

## Usage

Note that the current version assumes the chromosome order in reference FASTA and the BAM header are exactly the same.

Available options can be obtained by running verifybam without any parameters.

```
VERIFYBAM: BAM validity checking tool.
Version 0.0.4
	Last update: June 17, 2021, build date: Mon Jun 28 16:18:00 +03 2021

	--input  [BAM file]  : Input file in sorted and indexed BAM format.
	--ref    [reference] : Reference genome in FASTA format.
	--mode               : Running mode. Default is sequential
		server       : Start verifybam and wait for tasks to process. Reference is required
		client       : Send a task to running verifybam server. Bam file is required
		sequential   : Run verifybam once. Reference and Bam file is required
	--threads            : Number of threads to run while processing bam file.
	--client             : Run in client mode. Only bam file is required
	--version            : Print version and exit.
	--help               : Print this help screen and exit.
```
