**HiCEnvironmentFinder** is a bash script to build a genomic environment around target gene using paired Hi-C reads.

### Running *HiCEnvironmentFinder*

To run HiCEnvironmentFinder use the following syntax:
`HiCEnvironmentFinder.sh [<Launch options>]`

Prerequisites:
* java (>= 1.8 version)
* python3
* [pandas](https://pandas.pydata.org/)
* [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)
* [samtools](http://www.htslib.org/doc/samtools-view.html)

### Usage example

Here is a bash script showing a typical usage of HiCEnvironmentFinder:

~~~
./HiCEnvironmentFinder.sh --reads "example/wgs_reads/*.fastq" \
		--seq "example/seq.fasta" \
		--hi-c-r1 "example/hic_R1.fastq" \
		--hi-c-r2 "example/hic_R2.fastq" \
		--work-dir "example_work_dir" \
		--metacherchant "../out/metacherchant.jar" \
		--k 31 \
		--coverage 5 \
		--maxradius 100000
~~~


* `--reads` list of all input files with metagenomic reads separated by space. FASTA and FASTQ formats are supported.
* `--seq` a FASTA file with the target nucleotide sequences, for each of which a genomic environment will be built.
* `--hi-c-r1` and `--hi-c-r2` two input files with paired Hi-C reads. This parameters assumes the i-th read in hic_R1.fastq and the i-th read in hic_R2.fastq constitute a read pair. FASTA and FASTQ formats are supported.
* `--work-dir` working directory with intermediate files, logs and output folder.
* `--metacherchant` path to metacherchant jar file.
* `--k` --- the size of k-mer used in de Bruijn graph. It is optional parameter. Default value is 31.
* `--coverage` the minimum coverage threshold for a k-mer to be included in the graph. It is optional parameter. Default value is 5.
* `--maxradius` maximum allowed distance between every k-mer and taret gen. It is optional parameter. Default value is 100000.

### Output description

After the end of the analysis, the results can be found in the folder specified in `--work-dir` parameter.

* `output/2/merged/graph.gfa` - de Bruijn graph in GFA format. This environment was built using pared Hi-C reads. It is recommended to be viewed using [Bandage](http://rrwick.github.io/Bandage/). Name(s) of the node(s) corresponding to the target sequence ends with a suffix `_start`. Green color corresponds to the starting sequence. 

* `output/2/merged/env.txt` contains de Bruijn graph in a simple text format for later use: each line contains a k-mer and its coverage.

* `output/2/merged/seqs.fasta` - a FASTA file containing all sufficiently long unitigs (non-branching sequences in de Bruijn graph) for later analysis.

* `2/hic_map.txt` - mapping Hi-C reads to the contigs in the format: <contig id 1> <contig id 2> <hi-c weight> where hi-c weight is equal to count of Hi-C crosslinks between two contigs.


### Results visualisation

The obtained genomic environment can be visualized in [Bandage](http://rrwick.github.io/Bandage/). 

![example](https://user-images.githubusercontent.com/17966048/154980704-5d5b06cf-bedc-4a3c-b2f7-5e5691b7b340.png)
