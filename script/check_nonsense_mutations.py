#!/usr/bin/env python3
__version__ = "0.1"
__update__ = "09/23/2024"
__project__ = "codeathon to combat antimicrobial resistance at NCBI"

import os
import sys
import argparse
import logging
import re
from Bio.Blast import NCBIXML
from Bio.Seq import Seq


def check_nonsense_mutations(seq):
    # check nonsense mutations, could be multiple places
    if '*' in seq:
        return [m.start() for m in re.finditer(re.escape('*'), seq)]
    else:
        return None
    
def main(argv):
    result_handle = open(argv.input_xml,'r')
    blast_records = NCBIXML.parse(result_handle)
    out_header=["Query", 
                "Query_start", 
                "Query_end", 
                "Query_frame", 
                "MUTATION", 
                "HIT_start", 
                "HIT_end",
                "HIT_acc",
                "HIT_desc", 
                "AA %Identity", 
                "Coverage",  
                "E-value",
                "BitScore"]
    print("\t".join(out_header))
    for blast_record in blast_records:
        query_descs= blast_record.query.split(" ")
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                query_seq = Seq(hsp.query)
                hit_seq = Seq(hsp.sbjct)
                nonsense_mutations = check_nonsense_mutations(hsp.query)
                if nonsense_mutations:
                    mutations = ",".join([  hit_seq[i] + str(i + 1 + hsp.sbjct_start) + "STOP"  for i in nonsense_mutations]) 
                    #mutations =""
                    coverage = (hsp.sbjct_end - hsp.sbjct_start + 1 )/alignment.length 
                    identity = hsp.identities/hsp.align_length 
                    if coverage > argv.cov and identity > argv.id :
                        print("\t".join([query_descs[0],
                                        str(hsp.query_start),
                                        str(hsp.query_end),
                                        str(hsp.frame[0]),
                                        mutations,
                                        str(hsp.sbjct_start), 
                                        str(hsp.sbjct_end), 
                                        alignment.hit_id,
                                        alignment.hit_def, 
                                        str(f"{hsp.identities}/{hsp.align_length} ({identity*100:.2f}%)"), 
                                        str(f"{hsp.sbjct_end - hsp.sbjct_start + 1}/{alignment.length} ({coverage * 100:.2f}%)"), 
                                        str(hsp.expect),
                                        str(hsp.bits)])
                            )

if __name__ == '__main__':
    toolname = os.path.basename(__file__)
    argv = argparse.ArgumentParser( prog=toolname,
        description = "check stop codon from blastx result xml file",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    argv.add_argument('-i', '--input', dest = 'input_xml', required = True,
        help = 'blastout.xml')
    argv.add_argument('--cov', type=float, dest = 'cov', default = .8 , required=False, help='target coverage')
    argv.add_argument('--id',  type=float, dest = 'id',  default = .8 , required=False, help='hit identity')
    argv.add_argument('--verbose', action='store_true', help='Show more information in log')
    argv.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))

    argvs = argv.parse_args()

    log_level = logging.DEBUG if argvs.verbose else logging.INFO
    logging.basicConfig(
        format='[%(asctime)s' '] %(levelname)s: %(message)s', level=log_level, datefmt='%Y-%m-%d %H:%M')
    
    sys.exit(not main(argvs))
