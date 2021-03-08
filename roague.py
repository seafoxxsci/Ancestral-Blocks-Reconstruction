#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : ROAGUE pipeline, which includes finding orthoblocks in a set of given taxa and reconstructing the ancestral 
              gene blocks. This is basically a program that calls several other programs.
    Start   : 09/08/2016
    End     : 09/1/2016
    
    Iddo Friedberg:
    bcolors for printing
    added createdb: creating a BLAST database as an option
'''
import sys, argparse, os
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

## Traverses the genome information directory
def traverseAll(path):
    res = []
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+'/'+f)
    return res    
    

def parser_code():
    parser = argparse.ArgumentParser()            
    parser.add_argument("--genomes_directory", "-g", 
        help="The directory that stores the genome files (E_Coli/genomes)")  
    parser.add_argument("--gene_blocks", "-b", 
        help="The gene_block_names_and_genes.txt file, this file stores the operon name and its set of genes") 
    parser.add_argument("--reference", "-r", 
        help="The NCBI accession number for the reference genome (e.g. NC_000913 for E_Coli and NC_000964 for B_Sub)")
    parser.add_argument("--name", "-n",
        help="The common name of the reference organism (e.g. Escherichia_coli or Bacillus_subtilis")  
    parser.add_argument("--filter", "-f", 
        help="The filter file for creating the tree (E_Coli/phylo_order.txt for E_Coli or B_Sub/phylo_order.txt for B-Sub)")  
    parser.add_argument("--method", "-m", 
        help="The method to reconstruct the ancestral gene blocks, we support either global or local")       
    parser.add_argument("--output", "-o", 
        help="Output directory to store the results",default = "results")
    parser.add_argument("--approx", "-a", 
        help="Using approx method (T,F)", default = "F")
    parser.add_argument("--createdb", "-d", 
        help="Create a BLAST database (T,F)",default = "T")
    return parser.parse_args()


if __name__ == '__main__':
    args = parser_code()
    genomes_directory = args.genomes_directory
    reference = args.reference
    filter_file = args.filter
    name = args.name
    method = args.method
    gene_block_names_and_genes = args.gene_blocks
    outdir = args.output
    approx = args.approx
    createdb = args.createdb
    # check if we are going to output results in the current directory
    try:
        dirs = genomes_directory.split('/')
    except AttributeError:
        print ("Genomes directory argument can't not be None, please refer to the manual of the program (typing ./roague -h)")
    outdir += '/'
    try:
        os.mkdir(outdir+'/')
    except:
        print (bcolors.WARNING+"output directory has already been created"+bcolors.ENDC)
    if len(dirs) >= 3: # means that we have to go to subdirectory
        parent_dir = outdir+dirs[0]+"/"
    else:
        parent_dir = outdir
    print("PARENT DIR "+parent_dir)
    
    ### format a database for faster blasting, output in the db subdirectory
    db_dir = parent_dir+'db'
    if createdb.upper() == "T":
        print (bcolors.OKGREEN+"Creating BLAST db"+bcolors.ENDC)
        cmd1 = './format_db.py -i {} -o {}'.format(genomes_directory, db_dir)
        print (bcolors.OKCYAN+'Command 1:',cmd1+bcolors.ENDC)
        rc = os.system(cmd1)
        if rc != 0:
            sys.exit("{}Error in creating BLAST db{}".format(bcolors.FAIL, bcolors.ENDC))
    else:
        if os.path.isdir(db_dir):
            print (bcolors.WARNING+"BLAST db is already created"+bcolors.ENDC)
        else:
            raise IOError("db directory not found")
    
    ### Given the gene_block_names_and_genes.txt, create a gene_block_query.fa using the reference gene bank file. output in file gene_block_query.fa
    gene_block_names_and_genes = dirs[0]+"/"+'gene_block_names_and_genes.txt'
    gene_block_query = parent_dir+'gene_block_query.fa'
    cmd2 ='./make_operon_query.py -i {} -b {} -r {} -o {}'.format(genomes_directory, gene_block_names_and_genes, reference, gene_block_query)
    print (bcolors.OKGREEN+"Querying operons"+bcolors.ENDC)
    print (bcolors.OKCYAN+'Command 2:',cmd2+bcolors.ENDC)
    rc = os.system(cmd2)
    if rc != 0:
        sys.exit("{}Error in make_operon_query{}".format(bcolors.FAIL, bcolors.ENDC))
    
    ### blasting using db vs the gene_block_query.fa above. output in blast_result
    blast_result = parent_dir+'blast_result/'
    cmd3 = './blast_script.py -u {} -d {} -o {}'.format(gene_block_query, db_dir, blast_result)
    print (bcolors.OKCYAN+'Command 3:',cmd3)
    print ("{}BLASTing against genomes{}".format(bcolors.OKGREEN,bcolors.ENDC))
    rc = os.system(cmd3)
    if rc != 0:
        sys.exit("{}Error in blast_script{}".format(bcolors.FAIL, bcolors.ENDC))
    
    ### parsing the blast result directory into files that group by operon names, output in blast_parse
    blast_parse = parent_dir+'blast_parse/'
    cmd4 = './blast_parse.py -b {} -i {} -o {}'.format(gene_block_names_and_genes, blast_result, blast_parse)
    print (bcolors.OKGREEN+"Parsing BLAST results"+bcolors.ENDC)
    print (bcolors.OKCYAN+'Command 4:',cmd4+bcolors.ENDC)    
    rc = os.system(cmd4)
    if rc != 0:
        sys.exit("{}Error in blast_parse{}".format(bcolors.FAIL, bcolors.ENDC))
 
    ### filtering the gene blocks so that we have the most optimal gene blocks given the blast parse directory, ouput to optimized_gene_block
    optimized_gene_block = parent_dir+'optimized_gene_block/'
    cmd5 = './filter_operon_blast_results.py -i {} -o {} -e {}'.format(blast_parse, optimized_gene_block, 1e-3)
    os.system(cmd5)
    print (bcolors.OKGREEN+"Filtering gene blocks"+bcolors.ENDC)
    print (bcolors.OKCYAN+'Command 5:',cmd5+bcolors.ENDC)

    ### create newick tree file. output in tree directory
    tree_dir = parent_dir+'tree/'
    cmd6 = './create_newick_tree.py -G {} -f {} -r {} -o {}'.format(genomes_directory, filter_file, reference, tree_dir)
    print (bcolors.OKGREEN+"Create NEWICK tree file"+bcolors.ENDC)
    os.system(cmd6)
    print (bcolors.OKCYAN+'Command 6:',cmd6+bcolors.ENDC)    

    ### from the filter_operon_blast_results, create a result directory. output in result directory
    accession = tree_dir + 'accession_to_common.csv'
    result = parent_dir + 'result/'
    cmd7 = './get_result.py -g {} -a {} -o {} -i {}'.format(gene_block_names_and_genes, accession,result, optimized_gene_block)
    os.system(cmd7)
    print (bcolors.OKCYAN+'Command 7:',cmd7+bcolors.ENDC)
 
    ### generate group.txt file to color taxa based on class
    group = parent_dir+'group.txt'
    cmd8 ='./group.py -i {} -o {} -a {}'.format(genomes_directory, group, filter_file) 
    os.system(cmd8)
    print (bcolors.OKCYAN+'Command 8:',cmd8+bcolors.ENDC)    
    
    ### convert the file in directory result and store it in new_result
    new_result = parent_dir+'new_result/'
    cmd9 = './convert.py -i {} -o {}'.format(result, new_result)
    os.system(cmd9)
    print (bcolors.OKCYAN+'Command 9:',cmd9+bcolors.ENDC)    
    if approx.upper() == "F": 
        ### reconstructing the gene block using global method
        reconstruction = parent_dir+'reconstruction_'+method
        tree = tree_dir+'out_tree.nwk'
        cmd10 = './reconstruction.py -i {} -t {} -o {} -m {}'.format(new_result, tree, reconstruction, method)
        os.system(cmd10)
        print (bcolors.OKCYAN+'Command 10:',cmd10+bcolors.ENDC)
        
        ### create a visualization of the reconstruction and store it in directory visualization
        visualization = parent_dir+'/visualization/'
        try:
            os.mkdir(visualization)
        except:
            print (bcolors.WARNING+"Visualization directory is already created"+bcolors.ENDC)
        res = traverseAll(reconstruction)
        # go through all the reconstruction ancestral file and create a pdf visualization of them, using the group text to color group
        for file in res:
            operon = file.split('/')[-1]
            if 'mapping' in operon:
                continue
            outfile = visualization+operon
            cmd11 = './show_tree.py -i {} -g {} -o {} -r {}'.format(file, group, outfile, name)
            os.system(cmd11)
            print (bcolors.OKCYAN+'Command 11:',cmd11+bcolors.ENDC)
    else: 
        ### recosntruct the gene block using approximations
        new_result = new_result[:-1]+"_approx/"
        ### reconstructing the gene block using the global method
        reconstruction = parent_dir+'reconstruction_'+method+"_approx/"
        tree = tree_dir+'out_tree.nwk'
        cmd10 = './reconstruction.py -i {} -t {} -o {} -m {}'.format(new_result, tree, reconstruction, method)
        os.system(cmd10)
        print (bcolors.OKCYAN+'Command 10:',cmd10+bcolors.ENDC)
        
        ### create a visualization of the reconstruction and store it in directory visualization
        visualization = parent_dir+'/visualization_approx/'
        try:
            os.mkdir(visualization)
        except:
            print (bcolors.WARNING+"Visualization directory is already created"+bcolors.ENDC)
        res = traverseAll(reconstruction)
        # go through all the reconstruction ancestral file and create a pdf visualization of them, using the group text to color group
        for file in res:
            operon = file.split('/')[-1]
            if 'mapping' in operon:
                continue
            outfile = visualization+operon
            cmd11 = './show_tree.py -i {} -g {} -o {} -r {}'.format(file, group, outfile, name)
            os.system(cmd11)    
            print (bcolors.OKCYAN+'Command 11:',cmd11+bcolors.ENDC)
        
