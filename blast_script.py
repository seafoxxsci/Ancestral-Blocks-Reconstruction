#!/usr/bin/env python

'''
Copyright(C) 2015 David Ream
Released under GPL version 3 licence. http://www.gnu.org/licenses/lgpl.html
Do not remove this comment

TODO: Fix --query option to read in gene block query, then take two files, protein and nucleic acid queries and run them
      (This would require that there are two seperate, and complementary, blast databases within this folder. 
      My desire is that it would be two subfolders, 'protein/' and 'rna/' which house these sets of data.

TODO: Fix to accept protein, nucleotide, and other blast queries
    #cmd="blastall -p tblastn -a %i -i %s -d %s -e %s -o %s -m 9" % (num_processors, query_file, db, eval_threshold, out_file)
    #cmd="blastall -p blastp -a %i -i %s -d %s -e %s -o %s -m 9" % (1, query_file, db, eval_threshold, out_file)
'''
import time, os, sys, argparse
from multiprocessing import Pool


def parser_code():
    parser=argparse.ArgumentParser(description='Run a BLAST -m8 search on a database directory [./db/] using a query \
        [gene_block_query.fa] that built based on the gene blocks and reference genome specified by the user. Results \
        are organized by database accession numbers in a separate output directory [./blast_result/]')
    parser.add_argument("--database_folder", "-d", dest="database_folder", metavar="BLAST_DB", default='./db/',
        help="Input BLAST database folder containing all BLAST searchable files to be used by the program [default: ./db/]")          
    parser.add_argument("--outfolder", "-o", dest="outfolder", metavar="OUTPUT_DIR", default='./blast_result/',
        help="Output BLAST results folder organized by database accession [default: ./blast_result]")
    parser.add_argument("--filter", "-f", dest="filter", metavar="FILTER_FILE", default='NONE',
        help="Optional input file specifying which accession numbers this script will process (each accession \
        must be a new line). If no file is provided, filtering is not performed [default: NONE]")                    
    parser.add_argument("--num_proc", "-n", dest="num_proc", metavar="INTEGER", default=os.sysconf("SC_NPROCESSORS_CONF"), 
        type=int, help="Number of processors this script will run on [default: every CPU that the system has]")
    parser.add_argument("--query", "-u", dest="query", default='gene_block_query.fa', metavar="QUERY_FILE",
        help="Input BLAST query file containing every gene of interest in the dataset (must be in fasta format with > \
        lines delimiting each query) [default: gene_block_query.fa]")
    parser.add_argument("--eval", "-e", dest="eval", default='1e-6', metavar="E_VALUE", type=float,
        help="Input e-value parameter for the BLAST search [default: 1e-06]")
    parser.add_argument("--quiet", "-q", dest="quiet", action="store_true", default=False,
        help="Suppresses most program text outputs [default: FALSE]") 
    return parser.parse_args()


def check_options(parsed_args):
    # Input BLAST database folder
    if os.path.isdir(parsed_args.database_folder):
        database_folder=parsed_args.database_folder
    else:
        print(("The database directory %s does not exist." % parsed_args.database_folder))
        sys.exit()

    # Input query file
    if os.path.exists(parsed_args.query):
        query_file=parsed_args.query
    else:
        print(("The query file %s does not exist." % parsed_args.query))
        sys.exit()

    # Input filter file
    if parsed_args.filter == 'NONE' or os.path.exists(parsed_args.filter):
        filter_file=parsed_args.filter
    else:
        print(("The filter file %s does not exist." % parsed_args.filter))
        sys.exit()

    # Input e-value threshold    
    e_val=parsed_args.eval

    # Output folder
    if not os.path.isdir(parsed_args.outfolder):
        os.makedirs(parsed_args.outfolder)
    outfolder=parsed_args.outfolder
    
    # CPU cores that will be used by the program
    if parsed_args.num_proc > os.sysconf("SC_NPROCESSORS_CONF"):
        num_proc=os.sysconf("SC_NPROCESSORS_CONF")
    elif parsed_args.num_proc < 1:
        num_proc=1
    else:
        num_proc=int(parsed_args.num_proc)

    # Quiet TRUE/FALSE toggle
    quiet=parsed_args.quiet
    
    return database_folder, outfolder, filter_file, num_proc, query_file, e_val, quiet
 
        
def returnRecursiveDirFiles(root_dir):
    result=[]
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname=os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result


def do_parallel_blast(arg_tuple):
    db, query_file, blast_result_folder, num_processors, eval_threshold=arg_tuple
    out_file="%s%s.txt" % (blast_result_folder, db.split('/')[-1].split('.')[0])
    cmd ="blastp -num_threads %i -query %s -db %s -evalue %s -out %s -outfmt 6 -seg yes" % (1, query_file, db, eval_threshold, out_file)
    os.system(cmd)


def parallel_blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val):
    unfiltered_db_list=[i for i in returnRecursiveDirFiles(database_folder) if i.split('/')[-1].split('.')[-1] == 'ffc']
    if filter_file == '' or filter_file == 'NONE' or len(filter_file)<5:
        db_list=unfiltered_db_list
    else:
        filter_list=[i.strip() for i in open(filter_file).readlines()]
        db_list=[i for i in unfiltered_db_list if i.split('/')[-1].split('.')[0] in filter_list]
    
    blast_arg_list=[(i, query_file, outfolder, 1, e_val) for i in db_list]
    pool=Pool(processes=num_proc)
    pool.map(do_parallel_blast, blast_arg_list)


def main():
    start=time.time()
    parsed_args=parser_code()
    database_folder, outfolder, filter_file, num_proc, query_file, e_val, quiet=check_options(parsed_args)
    if not quiet:
        print((database_folder, outfolder, filter_file, num_proc, query_file, e_val, quiet))
    parallel_blast(database_folder, outfolder, filter_file, num_proc, query_file, e_val)
    if not quiet:
        print((time.time() - start))
    
    
if __name__ == '__main__':
    main()

