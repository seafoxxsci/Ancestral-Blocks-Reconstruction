#!/usr/bin/env python

import shutil
dirlist = ["blast_parse",
	"blast_result",
	"db",
	"gene_block_query.fa",
	"new_result",
	"new_result_approx",
	"optimized_gene_block",
	"reconstruction_global",
	"result",
	"tree",
	"visualization"]
"""
Removes the results from a ROAGUE run, with the option of keeping the BLAST db
"""
def remove_results(results_root_dir,interactive=True,remove_db=False):
    for target_dir in dirlist:
        if target_dir == 'db' and not remove_db:
            continue
        elif interactive:
            rmdir = "{}/{}".format(results_root_dir,target_dir)
            rc = input("Remove {}/{}? [y/N]".format(rmdir))
            if rc.upper() == "Y":
                shutil.rmtree(rmdir)
        else:
            rmdir = "{}/{}".format(results_root_dir,target_dir)
            shutil.rmtree(rmdir)

