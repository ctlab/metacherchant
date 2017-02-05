**MetaFast-Env** is a tool for analysing a genomic environment of a nucleotide sequence
within a metagenome, which is based on [MetaFast](https://github.com/ctlab/metafast/wiki) source code.

========

### Installation

To run MetaFast-Env you need to have JRE 1.6 or higher installed and have either of
these three files: `metafast.sh` for Linux/MacOS, `metafast.bat` for Windows or `metafast.jar` for any OS.

### Running ***metafast-env***

To run MetaFast-Env use the following syntax:
* `metafast.sh [<Launch options>]`
* `metafast.bat [<Launch options>]`
* `java -jar metafast.jar [<Launch options>]`

### Usage example

Here's the bash script showing the typical usage of MetaFast-Env:

~~~
./metafast.sh --tool environment-finder \
	--k 31 \
	--coverage=5 \
	--reads $READS_DIR/*.fasta \
	--seq $GENE_FILE.fasta \
	--output $OUTPUT_DIR/output/ \
	--work-dir $OUTPUT_DIR/workDir \
	--maxkmers=100000 \
	--bothdirs=False \
	--chunklength=100
~~~

* `--k` parameter sets the size of k-mer used in de Bruijn graph.
* `--coverage` parameter sets the lower threshold for k-mers to be included in the graph.
* `--reads` should list all files with reads, separated by space. Most modern formats are supported.
* `--seq` contains .fasta file with all nucleotide sequences, for which the genomic environment will be built.
* `--output` specifies an output folder, where all the results will be put in.
* `--work-dir` specifies a working directory with intermediate files and logs.
* `--maxkmers` sets the maximum number of k-mers present in resulting genomic environment.
* `--bothdirs` chooses between 2 modes of the BFS algorithm: False makes 2 one-directional passes and True makes 1 bidirectional pass.
* `--chunklength` determines the minimum length for a contracted graph node to be included in output .fasta file for further analysis.

There is also a _multi_ mode, where you can take multiple graph built as described above and join them into a single graph. Here's how:

~~~
./metafast.sh \
	--tool environment-finder-multi \
	--seq OXA-347.fasta \
	--work-dir "k31/TUE-S2_3_4/workDir" \
	--output "k31/TUE-S2_3_4" \
	--env "k31/TUE-S2_3/output/env.txt" "k31/TUE-S2_4/output/env.txt" # список из двух файлов env.txt, полученных другим запуском
~~~

* `--env` parameter lists all `env.txt` files to be joined into a single graph.

### Output description

After you run the tool, you'll find the result in folder specified in `--output` parameter (if there were multiple sequences in file in --seq,
there will be separate folder for each one).

* `graph.gfa` - de Bruijn graph in GFA format, which is adviced to be viewed using latest [Bandage](http://rrwick.github.io/Bandage/) version. 
Starting sequence corresponds to a node with a suffix `_start` (you can use Bandage's 'Find Nodes' feature to select them all). For a single graph, there is no special coloring. In graph for exactly 2 environments, green color corresponds to the starting sequence, red - nodes that are only in the first graph, blue - nodes that are only in the second graph, black nodes are in both. For 3 or more graph, nodes are in greyscale, darker shade of the node corresponds to more graphs that node is contained in. 

* `env.txt` contains de Bruijn graph in text format for later use. Format is a set of lines containing each k-mer and its read depth.

* `seqs.fasta` - .fasta file with all long enough unitigs (non-branching sequences in de Bruijn graph) for later analysis.

* `tsvs/*` - graph description in .tsv format for use in [Cytoscape](http://www.cytoscape.org/) tool.
