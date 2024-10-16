#!/usr/bin/env python3
__version__ = "0.2"
__update__ = "09/25/2024"
__project__ = "codeathon to combat antimicrobial resistance at NCBI"

import os
import sys
import argparse
import logging
import re
from Bio.Blast import NCBIXML
from Bio.Seq import Seq

def check_frameshift(query_seq):
    # run diamond blastx with options "--min-orf 1 --frameshift 15"
    return [m.start() for m in re.finditer(r"\\|/", query_seq)] if '/' in query_seq or '\\' in query_seq else [] 

def check_INDELs(query_seq, hit_seq, hit_start):
    # in-frame frameshifts 
    mutations = []
    locations = []
    if '-' in hit_seq:  # deletion
        for m in re.finditer('-+', hit_seq):
            mutations.append(hit_seq[m.start()-1] + str(m.start() + hit_start - 1) + query_seq[m.start()-1:m.end()])
            locations.append( str(m.start() + hit_start - 1))
    if '-' in query_seq: # insertion
        for m in re.finditer('-+', query_seq):
            mutations.append(hit_seq[m.start()-1:m.end()] + str(m.start() + hit_start - 1) + query_seq[m.start()-1])
            locations.append( str(m.start() + hit_start - 1))
    return mutations, locations

def check_nonsense_mutations(query_seq):
    # check nonsense mutations, could be multiple places
    if '*' in query_seq:
        return [m.start() for m in re.finditer(re.escape('*'), query_seq)]
    else:
        return []
    
def main(argv):
    result_handle = open(argv.input_xml,'r')
    blast_records = NCBIXML.parse(result_handle)
    out_header=["element_symbol",
                "contig_acc", 
                "contig_start", 
                "contig_stop", 
                "orientation", 
                "lesion",
                "lesion_type", 
                "aa_identity",
                "nt_identity",
                "ref_accession",
                "ref_desc", 
                "lesion_location",
                "ref_start", 
                "ref_stop",
                "coverage",  
                ]
    print("\t".join(out_header))
    for blast_record in blast_records:
        query_descs= blast_record.query.split(" ")
        contig_id = query_descs.pop(0)
        for alignment in blast_record.alignments:
            hit_names = alignment.hit_id.split("_")
            element = hit_names.pop(0)
            for hsp in alignment.hsps:
                query_seq = hsp.query
                hit_seq = hsp.sbjct
                coverage = (hsp.sbjct_end - hsp.sbjct_start + 1 )/alignment.length 
                identity = hsp.identities/hsp.align_length
                orientation = '+' if hsp.frame[0] > 0 else '-'
                if coverage > argv.cov and identity > argv.id : 
                    indels_mutations, indels_locations = check_INDELs(query_seq,hit_seq, hsp.sbjct_start) if argvs.indels else [],[]
                    nonsense_mutations = check_nonsense_mutations(query_seq)
                    frameshift_mutations = check_frameshift(query_seq)
                    
                    mutations = []
                    lesion_type =[]
                    lesion_location = []
                    if len(nonsense_mutations) > 0:
                        lesion_type = ['STOP']
                        lesion_location = [  str(i + hsp.sbjct_start)  for i in nonsense_mutations ] 
                        mutations = [  hit_seq[i] + str(i + hsp.sbjct_start) + "STOP"  for i in nonsense_mutations ] 
                    if len(frameshift_mutations) > 0:
                        lesion_type.extend(['FRAMESHIFT'])
                        lesion_location.extend([  str(i + hsp.sbjct_start - query_seq.count('-',0,i))  for i in frameshift_mutations ])
                        mutations.extend([  hit_seq.replace('-','')[i] + str(i + hsp.sbjct_start  - query_seq.count('-',0,i)) + "fs"  for i in frameshift_mutations ])
                    elif len(indels_mutations) > 0:
                        lesion_type.extend(['FRAMESHIFT'])
                        lesion_location.extend(indels_locations)
                        mutations.extend(indels_mutations)

                    if len(nonsense_mutations) > 0 or len(indels_mutations) > 0 or len(frameshift_mutations) > 0:
                        print("\t".join([   element,
                                            contig_id,
                                            str(hsp.query_start),
                                            str(hsp.query_end),
                                            orientation,
                                            ','.join(mutations),
                                            ','.join(lesion_type),
                                            str(f"{hsp.identities}/{hsp.align_length} ({identity*100:.2f}%)"), 
                                            '',
                                            '_'.join(hit_names),
                                            alignment.hit_def, 
                                            ','.join(lesion_location),
                                            str(hsp.sbjct_start), 
                                            str(hsp.sbjct_end), 
                                            str(f"{hsp.sbjct_end - hsp.sbjct_start + 1}/{alignment.length} ({coverage * 100:.2f}%)")
                                            ])
                                )
    result_handle.close()

if __name__ == '__main__':
    toolname = os.path.basename(__file__)
    argv = argparse.ArgumentParser( prog=toolname,
        description = r'''

    ____                 _    ____                  __        __    _ _    _             
   |  _ \  ___  __ _  __| |  / ___| ___ _ __   ___  \ \      / __ _| | | _(_)_ __   __ _ 
   | | | |/ _ \/ _` |/ _` | | |  _ / _ | '_ \ / _ \  \ \ /\ / / _` | | |/ | | '_ \ / _` |
   | |_| |  __| (_| | (_| | | |_| |  __| | | |  __/   \ V  V | (_| | |   <| | | | | (_| |
   |____/ \___|\__,_|\__,_|  \____|\___|_| |_|\___|    \_/\_/ \__,_|_|_|\_|_|_| |_|\__, |
                                                                                   |___/ 
    Check nonsensus/frameshift mutations from blastx result xml file
    ''',
        formatter_class = argparse.RawDescriptionHelpFormatter)
    argv.add_argument('-i', '--input', dest = 'input_xml', required = True,
        help = '-outfmt 5 -o blastx_out.xml')
    argv.add_argument('--cov', type=float, dest = 'cov', default = .8 , required=False, help='target coverage')
    argv.add_argument('--id',  type=float, dest = 'id',  default = .8 , required=False, help='hit identity')
    argv.add_argument('--indels', action='store_true', help='check in-frame insertion and deletion')
    argv.add_argument('--verbose', action='store_true', help='Show more information in log')
    argv.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))

    argvs = argv.parse_args()

    log_level = logging.DEBUG if argvs.verbose else logging.INFO
    logging.basicConfig(
        format='[%(asctime)s' '] %(levelname)s: %(message)s', level=log_level, datefmt='%Y-%m-%d %H:%M')
    
    sys.exit(not main(argvs))
