#!/usr/bin/env python

'''
Copyright(C) 2015 David Ream
Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
Do not remove this comment
'''
from multiprocessing import Pool
import time, os, sys, argparse
from homolog4 import *
from collections import defaultdict


def parser_code():
    parser = argparse.ArgumentParser(description = "Organize the results of a BLAST -m8 search [./blast_result/] by a \
        gene block query file [gene_block_names_and_genes.txt] and output the results into a user-specified directory \
        [./blast_parse/]")              
    parser.add_argument("--infolder", "-i", dest = "infolder", metavar = "BLAST_DB_DIR", default = './blast_result/', 
        help = "Input BLAST database folder containing every organism that you are interested in searching")  
    parser.add_argument("--outfolder", "-o", dest = "outfolder", metavar = "OUTPUT_DIR", default = './blast_parse/', 
        help = "Output BLAST results folder organized by gene block [default: ./blast_result/]")
    parser.add_argument("--gene_block_query", "-b", dest = "gene_block_query", default = 'gene_block_names_and_genes.txt', 
        help = "Input gene block query file containing the tab-delimited gene block names and genes under investigation \
        [default: gene_block_names_and_genes.txt]")
    parser.add_argument("--filter", "-f", dest = "filter", metavar = "FILTER_FILE", default = '',
        help = "Optional input file containing NCBI accession numbers of organisms to include in the parsed BLAST \
        results (each accession must be a new line) [default: NONE]")            
    parser.add_argument("--num_proc", "-n", dest = "num_proc", type = int, default = os.sysconf("SC_NPROCESSORS_CONF"), 
        help = "Number of processors this script will run on [default: every CPU that the system has]")
    parser.add_argument("--quiet", "-q", dest = "quiet", action = "store_true", default = False,
        help = "Suppresses most program text outputs [default: FALSE]")
    return parser.parse_args()


def check_options(parsed_args):
    # Input folder
    if os.path.isdir(parsed_args.infolder):
        infolder = parsed_args.infolder
    else:
        print("The input BLAST database directory %s does not exist." % parsed_args.infolder)
        sys.exit()

    # Input gene block query    
    if os.path.exists(parsed_args.gene_block_query):
        gene_block_query = parsed_args.gene_block_query
    else:
        print("The gene block query file %s does not exist." % parsed_args.gene_block_query)
        sys.exit()

    # Input filter file
    if os.path.exists(parsed_args.filter):
        filter_file = parsed_args.filter
    elif parsed_args.filter == '':
        filter_file = parsed_args.filter
    else:
        print("The filter file %s does not exist." % parsed_args.filter)
        sys.exit()

    # Output folder
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder = parsed_args.outfolder
    if outfolder[-1] != '/':
        outfolder = outfolder + '/'
    
    # CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc = os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc = 1
    else:
        num_proc = int(parsed_args.num_proc)

    # Quiet TRUE/FALSE toggle
    quiet = parsed_args.quiet

    return infolder, outfolder, gene_block_query, filter_file, num_proc, quiet


def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


# The filter file here i think should be either a vacant value (such as '') or a user defined
# value.  I do not think by default it should be given, like I have made the default behavior.    
def parallel_blast_parse_dict(in_folder, out_folder, num_proc, filter_file, gene_block_dict):
    result = {}
    gene_block_out_folder = out_folder
    if not os.path.isdir(gene_block_out_folder):
        os.makedirs(gene_block_out_folder)
    if filter_file != '':
        tmp = returnRecursiveDirFiles(in_folder)
        nc_list = [i.strip() for i in open(filter_file).readlines()]
        fname_list = [i for i in tmp if i.split('/')[-1].split('.')[0] in nc_list]
    else:
        fname_list = returnRecursiveDirFiles(in_folder)
        
    for fname in fname_list:
        for line in [i.strip() for i in open(fname).readlines()]:
            try:
                hlog = Homolog.from_blast(line)
            except:
                print("Error in function parallel_blast_parse_dict from script blast_parse.py, conversion from result \
                    to Homolog class failed.", line)
            try:
                accession = hlog.accession()
            except:
                print("There was an error in the Homolog class, found in line function parallel_blast_parse_dict from \
                    script blast_parse.py")
            predicted_gene = hlog.blast_annotation()
            
            try: # faster implementation than "if predicted_gene in gene_block_dict.keys():"
                gene_block = gene_block_dict[predicted_gene]
                if gene_block in list(result.keys()):  # check gene block is in the result dict already
                    if accession in list(result[gene_block].keys()):  # check organism has been added to the gene block
                        result[gene_block][accession].append(hlog.to_file())
                    else:  # if the organism is not part of the gene block, add it
                        result[gene_block].update({accession:[hlog.to_file()]})
                else:  # add the gene_block to the result
                    result.update({gene_block: {accession: [hlog.to_file()]}})
            except:
                pass    

    for gene_block in list(result.keys()):
        handle = open(out_folder + gene_block + '.txt', 'w')
        for accession in list(result[gene_block].keys()):
        	handle.write('\n'.join(result[gene_block][accession]) + '\n')
        handle.close()
  

'''     
@input:     provide a file name that contains all of the hits for every organism that we are interested in (hlog_list) 
            and sorts by organism, then by locus. Files have already been screened for eval and gene block membership.
@return:    dictionary for that gene block
TODO:       find best hit for the locus out of many
'''
def return_gene_block_list(fname):
    gene_block = fname.split('/')[-1].split('.')[0]
    hlog_list = [Homolog.from_file(i.strip()) for i in open(fname).readlines()]
    result_dict = {}
    for hlog in hlog_list:
        accession = hlog.accession()
        locus = hlog.locus()
        if accession not in list(result_dict.keys()):
            result_dict.update({accession: {}})
        if locus not in list(result_dict[accession].keys()):
            result_dict[accession].update({locus: [hlog]})
        else:
            result_dict[accession][locus].append(hlog)
    return result_dict
        

'''
currently not used... and might not use this: will see
def parallel_return_gene_block_list(infolder, outfolder, num_proc):
    pool = Pool(processes = num_proc)
    organism_dict_for_recovery = dict(pool.map(parallel_gene_block_fasta, genome_of_interest_list))
'''
# This function will take the organism-locus dict (per gene block file) and determine the best homolog.
def best_homolog_list(gene_block_dict, outfile):
    result = []
    
    for org in sorted(gene_block_dict.keys()):
        for locus in sorted(gene_block_dict[org].keys()):
            hlog_list = gene_block_dict[org][locus][1:]
            best_hit = gene_block_dict[org][locus][0]
            gene_count = defaultdict(int) # Identifies if a locus has >1 predicted gene, and the count ratio
            # print gene_count.keys()
            gene_count[best_hit.predicted_gene()] +=1
            for hlog in hlog_list:
                gene_count[hlog.predicted_gene()] +=1
                if best_hit.e_val() > hlog.e_val():
                    best_hit = hlog
            result.append(best_hit)
    handle = open(outfile, 'w')
    handle.write('\n'.join([i.ret_str() for i in result]))
    handle.close()


def return_gene_to_gene_block_dict(fname):
    gene_block_dict = {}
    for line in [i.strip().split('\t') for i in open(fname).readlines()]:
        gene_block = line[0]
        for gene in line[1:]:
            gene_block_dict.update({gene: gene_block})
    return gene_block_dict
 

def main():
    start = time.time()
    parsed_args = parser_code()
    infolder, outfolder, gene_block_query, filter_file, num_proc, quiet = check_options(parsed_args)
    if not quiet:
        print(infolder, outfolder, gene_block_query, filter_file, num_proc, quiet)
    
    # This code makes a dictionary mapping gene annotation to the gene block that it belong to
    gene_block_dict = return_gene_to_gene_block_dict(gene_block_query)
    parallel_blast_parse_dict(infolder, outfolder, num_proc, filter_file, gene_block_dict)

    if not quiet:
        print(time.time() - start)


if __name__ == '__main__':
    main()

