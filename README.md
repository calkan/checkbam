# verifybam

verifybam is a BAM file integrity verification tool.

Read mapping is verified by applying CIGAR and MD operations to both original read and corresponding sequence in reference genome. That's why it is crucial to use the same reference genome that is used in alignment process.

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

To install verifybam and use from shell.

## Usage

Available options can be obtained by running verifybam without any parameter.

```
VERIFYBAM: BAM validity checking tool.
Version 0.1-alpha
	Last update: September 21, 2015, build date: Mon Aug  1 21:01:34 EEST 2016

	--bamlist   [bamlist file] : A text file that lists input BAM files one file per line.
	--input [BAM files]        : Input files in sorted and indexed BAM format. You can pass multiple BAMs using multiple --input parameters.
	--ref   [reference genome] : Reference genome in FASTA format.
	--version                  : Print version and exit.
	--help                     : Print this help screen and exit.
```
