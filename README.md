# ROAGUE: **R**econstruction **o**f **A**ncestral **G**ene Blocks **U**sing **E**vents
## Purpose

ROAGUE is a tool to reconstruct ancestors of gene blocks in prokaryotic genomes. Gene blocks are genes co-located on the chromosome. In many cases, gene blocks are conserved between bacterial  pecies, sometimes as operons, when genes are co-transcribed. The conservation is rarely absolute: gene loss, gain, duplication, block splitting and block fusion are frequently observed. 

ROAGUE accepts a set of species and a gene block in a reference species. It then finds all gene blocks, orthologous to the reference gene blocks, and reconsructs their ancestral states.

## Dependencies
* [Conda](https://conda.io/miniconda.html) 
* [Python 3+](https://www.python.org/download/releases/3.0/)
* [Biopython 1.63+](http://biopython.org/wiki/Download)
* [ClustalW](http://www.clustal.org/clustal2/#Download)
* [Muscle](https://www.drive5.com/muscle/downloads.htm)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [ETE3](http://etetoolkit.org/download/) (python framework for tree)
* [PDA](http://www.cibiv.at/software/pda/#download) (optional if you want to debias your tree base on Phylogenetic Diversity)

## Installation
Type the following command line in your Terminal to copy the relevant scripts and examples:
```bash
git clone https://github.com/seafoxxsci/Ancestral-Blocks-Reconstruction
```

Install Miniconda and wget (if you don't have these tools already), and create a new environment for ROAGUE:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O Miniconda-latest-MacOSX-x86_64.sh
bash Miniconda-latest-MacOSX-x86_64.sh
conda create --name roague
conda activate roague
```

Add the following to your `.bashrc` file by typing `nano ~/.bashrc` and pasting the following line:
```bash
export PATH=~/miniconda3/envs/roague/bin:$PATH
```

Before installing any dependencies, you should set channel priorities (`conda-forge` > `bioconda` > `defaults`) in your `.condarc` file by typing `nano ~/.condarc` or by pasting the following 3 lines into your Terminal:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Go ahead and install the dependencies listed above using conda, and check your build for `ete3` which is needed for Newick tree generation and visualization:
```bash
conda install -c bioconda biopython blast clustalw 
conda install -c etetoolkit ete3 ete_toolchain
ete3 build check
```

Sometimes `ete3 build check` will return a message showing that some packages such as `clustalo` are missing. In these cases, download the package from source. Below is an example for Mac OSX:
```bash
wget http://www.clustal.org/omega/clustal-omega-1.2.3-macosx -O $PATH/ete3_apps/bin/clustalo
chmod u+x $PATH/ete3_apps/bin/clulstalo
```

In other cases, `ete3 build check` may return a message showing that all packages present, but you may encounter an error in `create_newick_tree.py` stating `sh: muscle not found`. This is a broad work-around:
```bash
cp $PATH/ete3_apps/bin/muscle $PATH/muscle
```

For the installation of PDA, please follow the instructions on this website: http://www.cibiv.at/software/pda/#download

## Example datasets and usage
The easiest way to run the project is to execute the script `roague.py` found inside the `Ancestral-Blocks-Reconstruction` directory. Running `roague.py --help` will output the following message which details the required inputs that you will need:
```
usage: roague.py [-h] [--genomes_directory GENOMES_DIRECTORY]
                 [--gene_blocks GENE_BLOCKS] [--reference REFERENCE]
                 [--filter FILTER] [--method METHOD] [--output OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  --genomes_directory GENOMES_DIRECTORY, -g GENOMES_DIRECTORY
                        The directory that store all the genomes file
                        (E_Coli/genomes)
  --gene_blocks GENE_BLOCKS, -b GENE_BLOCKS
                        The gene_block_names_and_genes.txt file, this file
                        stores the operon name and its set of genes
  --reference REFERENCE, -r REFERENCE
                        The ncbi accession number for the reference genome
                        (NC_000913 for E_Coli and NC_000964 for B_Sub)
  --filter FILTER, -f FILTER
                        The filter file for creating the tree
                        (E_Coli/phylo_order.txt for E_Coli or
                        B_Sub/phylo_order.txt for B-Sub)
  --method METHOD, -m METHOD
                        The method to reconstruc ancestral gene block, we
                        support either global or local
  --output OUTPUT, -o OUTPUT
                        Output directory to store the result
```
  
You can double-check your installation of ROAGUE by running `roague.py` on the following datasets, which were downloaded into the `Ancestral-Blocks-Reconstruction/E_Coli` and `Ancestral-Blocks-Reconstruction/B_Sub` directories. The finalized results, which are phylogenetic trees reconstructed from the gene block information that you provide, are by default stored in `result/E_Coli/visualization` and `result/B_Sub/visualization`:

### Escherichia coli K-12 MG1655
```bash
./roague.py --genomes_directory E_Coli/genomes/ --gene_blocks E_Coli/gene_block_names_and_genes.txt --reference NC_000913 --filter E_Coli/phylo_order.txt --method global
```

### Bacillus subtilis
```bash
./roague.py --genomes_directory B_Sub/genomes/ --gene_blocks B_Sub/gene_block_names_and_genes.txt --reference NC_000964 --filter B_Sub/phylo_order.txt --method global
```

### Bacterial gibberellin analysis
If you'd like to replicate the figures in our mSphere publication (Nett et al. 2020), simply download the example datasets and follow the instructions in the README at https://github.com/nguyenngochuy91/Gibberellin-Operon.

## Running ROAGUE on custom datasets
If you want to run `roague.py` on your own datasets, please keep in mind that the following inputs are required in a unique project directory within `Ancestral-Blocks-Reconstruction`:
  1. A `genomes` directory that contains all of the genomes you want to compare in GenBank format, including your reference species 
  2. A tab-delimited `gene_block_names_and_genes.txt` file that contains all of the gene blocks of interest in your reference species. The first column is the gene block name, followed by the names of each gene within the block. For example, here is the `gene_block_names_and_genes.txt` file from the `E_Coli/gene_block_names_and_genes.txt` file:
```bash
astCADBE	astA	astB	astC	astD	astE
atpIBEFHAGDC	atpI	atpH	atpC	atpB	atpA	atpG	atpF	atpE	atpD
caiTABCDE	caiA	caiE	caiD	caiC	caiB	caiT
casABCDE12	casE	casD	casA	casC	casB	cas1	cas2
chbBCARFG	chbG	chbF	chbC	chbB	chbA	chbR
``` 

Prior to running `./roague.py` on your own dataset, you should be aware that there may be families of species that are overrepresented in your genomes directory. This will reduce phylogenetic diversity and cause bias in the ancestral reconstruction. We recommend that you run [PDA](http://www.cibiv.at/software/pda/#download) on your generated tree before proceeding by doing the following:

1. Generate a phylogenetic tree from the genomes directory
```bash
./create_newick_tree.py -G genomes_directory -o tree_directory -f NONE -r ref_accession
usage: create_newick_tree.py [-h] [-G DIRECTORY] [-o DIRECTORY] [-f FILE] [-m STRING] [-t FILE] [-r REF] [-q]

optional arguments:
  -h, --help            show this help message and exit
  -G DIRECTORY, --genbank_directory DIRECTORY
                        Folder containing all genbank files for use by the
                        program.
  -o DIRECTORY, --outfolder DIRECTORY
                        Directory where the results of this program will be
                        stored.
  -f FILE, --filter FILE
                        File restrictiong which accession numbers this script
                        will process. If no file is provided, filtering is not
                        performed.
  -r REF, --ref REF     The reference genome number, such as NC_000913 for E_Coli
  -q, --quiet           Suppresses most program text outputs.
```

2. Debias the phylogenetic tree using PDA in the `./debias.py` script:
```bash
./debias.py -i tree_directory/out_tree.nwk -o pda_result.txt -s num -r ref_accession
usage: debias.py [-h] [-i INPUT_TREE] [-o PDA_OUT] [-s TREE_SIZE] [-r REF]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_TREE, --input_tree INPUT_TREE
                        Input tree that we want to debias
  -o PDA_OUT, --pda_out PDA_OUT
                        Output of pda to be store.
  -s TREE_SIZE, --tree_size TREE_SIZE
                        Reduce the size of the tree to this size
  -r REF, --ref REF     Force to include the following species, here I force
                        to include the reference species
```

3. Run ROAGUE,  the output is stored in directory `result`. 
```bash
./roague.py -g genomes_directory -b gene_block_names_and_genes.txt -r ref_accession -f phylo_order.txt -m global -o result
```

## Example outputs
Here are two gene blocks that were generated through our program. 
1. Gene block paaABCDEFGHIJK:

This gene block codes for genes involved in the catabolism of phenylacetate and is not conserved between the groups of studied bacteria.

![paaABCDEFGHIJK](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/images/paa_global_edit.png "Gene block paaABCDEFGHIJK")
2. Gene block atpIBEFHAGDC:

This gene block catalyzes the synthesis of ATP from ADP and inorganic phosphate and is very conserved between the groups of studied bacteria.

![atpIBEFHAGDC](https://github.com/nguyenngochuy91/Ancestral-Blocks-Reconstruction/blob/master/images/atp_global_edit.png "Gene block atpIBEFHAGDC")

## Credits
1. http://bioinformatics.oxfordjournals.org/content/early/2015/04/13/bioinformatics.btv128.full 

