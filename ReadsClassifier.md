Reads classifier tools are provided to compare two metagenomes.
Algorithm splits reads from one metagenome based on their presence in de Bruijn graph built from other metagenome.

## Usage example

### Simple reads classifier

Simple reads classifier splits reads into two categories.

Here is a bash script showing a typical usage of simple reads classifier:

~~~
java -jar metacherchant.jar --tool reads-classifier \
    -k 31 \
    -i <input_1.fasta input_2.fasta> \
    -r <reads_1.fasta reads_2.fasta> \
    -found 90 \
    -w <workDir> \
    -o <outDir>
~~~

* `-k` --- the size of k-mer used in de Bruijn graph.
* `-i` --- two files with paired reads for de Bruijn graph. FASTA and FASTQ formats are supported, as well as compressed files *.gz or *.bz2.
* `-r` --- two files with paired reads to be classified. FASTA and FASTQ formats are supported, as well as compressed files *.gz or *.bz2.
* `-found` --- Minimum coverage breadth for reads from class *found* [0 - 100 %].
* `-w` --- directory with intermediate working files
* `-o` --- directory for final categories of reads

#### Output description

After the end of the analysis, the results can be found in the folder specified in `-o` parameter

* `found_[1|2|s].fastq` --- reads which were found in de Bruijn graph built from input reads

* `not_found_[1|2|s].fastq` --- reads which were not found in de Bruijn graph built from input reads

### Triple reads classifier

Triple reads classifier splits reads into three categories. It utilizes two values of `k` and several
heuristic thresholds (can be user defined) for more accurate classification.

Here is a bash script showing a typical usage of triple reads classifier:

~~~
java -jar metacherchant.jar --tool triple-reads-classifier \
    -k 31 \
    -k2 61 \
    -i <input_1.fasta input_2.fasta>
    -r <reads_1.fasta reads_2.fasta>
    -found 90 \
    -half 40 \
    -w <workDir> \
    -o <outDir>
~~~

* `-k` --- the size of k-mer used in de Bruijn graph.
* `-k2` --- the second size of k-mer used in de Bruijn graph. k2 > k
* `-i` --- two files with paired reads for de Bruijn graph. FASTA and FASTQ formats are supported, as well as compressed files *.gz or *.bz2.
* `-r` --- two files with paired reads to be classified. FASTA and FASTQ formats are supported, as well as compressed files *.gz or *.bz2.
* `-found` --- Minimum coverage breadth for reads from class *found* [0 - 100 %].
* `-half` --- Minimum coverage breadth for reads from class *half-found* [0 - 100 %].
* `-w` --- directory with intermediate working files
* `-o` --- directory for final categories of reads

#### Output description

After the end of the analysis, the results can be found in the folder specified in `-o` parameter

* `found_[1|2|s].fastq` --- reads which were found in de Bruijn graph built from input reads

* `half_found_[1|2|s].fastq` --- reads which were partially found in de Bruijn graph built from input reads

* `not_found_[1|2|s].fastq` --- reads which were not found in de Bruijn graph built from input reads
