#!/usr/bin/env python3
__version__ = "0.3"
__update__ = "09/26/2024"
__project__ = "codeathon to combat antimicrobial resistance at NCBI"

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile

from Bio.Blast import NCBIXML

def check_diamond_db(args):
    database = os.path.expanduser(args.database)
    database = os.path.expandvars(database)
    dmnd_check = os.path.exists(database + ".dmnd") or os.path.exists(database)

    if dmnd_check:
        pass
    else:
        logging.info(f"DIAMOND database index not found for {args.database}")
        logging.info(f"Please run 'diamond makedb --in {args.database} --db {args.database}.dmnd'")
        sys.exit(1)
    
    try:
        # Run the diamond command to get version info
        result = subprocess.run(['diamond', '--version'], capture_output=True, text=True, check=True)
        
        # The version info is in the stdout of the result
        version_string = result.stdout.strip()
        
        return version_string
    except subprocess.CalledProcessError as e:
        logging.error(f"Error getting DIAMOND version: {e}")
        return f"Error getting DIAMOND version: {e}"
    except FileNotFoundError:
        logging.error("DIAMOND executable not found. Make sure it's installed and in your PATH")
        sys.exit(1)
        
def check_blast_db(args):
    database = os.path.expanduser(args.database)
    database = os.path.expandvars(database)

    psq_check = os.path.exists(database + ".psq")
    phr_check = os.path.exists(database + ".phr")
    pin_check = os.path.exists(database + ".pin")

    if all([psq_check, phr_check, pin_check]):
        pass
    else:
        logging.info(f"BLAST database indexes not found for {args.database}")
        logging.info(f"Please run 'makeblastdb -dbtype prot -in {args.database} -out {args.database}")
        sys.exit(1)

    try:
        # Run the BLAST command to get version info
        result = subprocess.run(['blastx', '-version'], capture_output=True, text=True, check=True)
        
        # The version info is in the stdout of the result
        version_string = result.stdout.strip()[0:15]
        
        return version_string
    except subprocess.CalledProcessError as e:
        logging.error(f"Error getting BLAST version: {e}")
        return f"Error getting BLAST version: {e}"
    except FileNotFoundError:
        logging.error("BLAST executable not found. Make sure it's installed and in your PATH.")
        sys.exit(1)
        

def run_blast(version: str, args, in_fasta: str, out_xml: str) -> None:
    logging.debug(version)
    logging.info("NCBI blastx executed with %s threads on %s against %s." % (args.threads, os.path.basename(in_fasta), os.path.basename(args.database)))

    blast_cmd = ["blastx",
                 "-query", in_fasta,
                 "-db", args.database,
                 "-num_threads", str(args.threads),
                 "-max_target_seqs", '100',
                 "-max_hsps", '10',
                 "-evalue", '0.001',
                 "-outfmt", '5',
                 "-out", out_xml.name ]
    logging.debug("CMD:" + " ".join(blast_cmd))
    try:
        subprocess.run(blast_cmd, check=True)
    except subprocess.CalledProcessError as error:
        error_message = error.stderr.decode('UTF-8')
        logging.error(f"Error when running BLAST:\n{error_message}")
        sys.exit(1)
    finally:
        out_xml.close()


def run_diamond(version: str, args, in_fasta: str, out_xml: str) -> None:
    logging.debug(version)
    logging.info("Diamond blastx executed with %s threads on %s against %s." % (args.threads, os.path.basename(in_fasta),os.path.basename(args.database)))
    database = args.database + ".dmnd" if os.path.exists(args.database + ".dmnd") else args.database
    diamond_cmd = ["diamond", "blastx",  
                   "--quiet",
                   "--query", in_fasta,
                   "--db", database, 
                   "--out", out_xml.name,
                   "--threads", str(args.threads),
                   "--max-target-seqs", str(100),
                   "--evalue", '0.001',
                   "--max-hsps", '10',
                   "--outfmt", '5',
                   "--min-orf", '1', 
                   "--frameshift", '15'
    ]
    logging.debug("CMD:" + " ".join(diamond_cmd))
    try:
        subprocess.run(diamond_cmd, check=True)
    except subprocess.CalledProcessError as error:
        error_message = error.stderr
        logging.error(f"Error when running DIAMOND:\n{error_message}")
        sys.exit(1)
    finally:
        out_xml.close()

def check_frameshift(query_seq):
    # run diamond blastx with options "--min-orf 1 --frameshift 15"
    return [m.start() for m in re.finditer(r"\\|/", query_seq)] if '/' in query_seq or '\\' in query_seq else [] 

def check_INDELs(query_seq, hit_seq, hit_start):
    # in-frame frameshifts 
    mutations = []
    if '-' in hit_seq:  # deletion
        for m in re.finditer('-+', hit_seq):
            mutations.append(hit_seq[m.start()-1] + str(m.start() + hit_start - 1) + query_seq[m.start()-1:m.end()])
    if '-' in query_seq: # insertion
        for m in re.finditer('-+', query_seq):
            mutations.append(hit_seq[m.start()-1:m.end()] + str(m.start() + hit_start - 1) + query_seq[m.start()-1])
    return mutations

def check_nonsense_mutations(query_seq):
    # check nonsense mutations, could be multiple places
    if '*' in query_seq:
        return [m.start() for m in re.finditer(re.escape('*'), query_seq)]
    else:
        return []
    
def main(argv):

    temp_xml = tempfile.NamedTemporaryFile(mode='w+', suffix='.xml', delete=False)
    if (argv.diamond):
        version=check_diamond_db(argv)
        run_diamond(version,argv,argv.input_fasta,temp_xml)
    else:
        version=check_blast_db(argv)
        run_blast(version,argv,argv.input_fasta,temp_xml)

    logging.info(f'Detect nonsensus/frameshift mutations from blast output')
    result_handle = open(temp_xml.name,'r')
    output_handle = open(argv.output, 'w') if (argv.output and argv.output != "-") else sys.stdout
    blast_records = NCBIXML.parse(result_handle)
    out_header=["element_symbol",
                "contig_acc", 
                "contig_start", 
                "contig_stop", 
                "orientation", 
                "lession_type", 
                "aa_identity",
                "nt_identity",
                "ref_accession",
                "ref_desc", 
                "ref_start", 
                "ref_stop",
                "coverage",  
                ]
    
    output_handle.write("\t".join(out_header)+"\n")
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
                if coverage > argv.cov and identity > argv.id : 
                    indels_mutations = check_INDELs(query_seq,hit_seq, hsp.sbjct_start) if argvs.indels else []
                    nonsense_mutations = check_nonsense_mutations(query_seq)
                    frameshift_mutations = check_frameshift(query_seq)
                    
                    mutations = []
                    if len(nonsense_mutations) > 0:
                        mutations = [  hit_seq[i] + str(i + hsp.sbjct_start) + "STOP"  for i in nonsense_mutations ] 
                    if len(frameshift_mutations) > 0:
                        mutations.extend([  hit_seq.replace('-','')[i] + str(i + hsp.sbjct_start) + "fs"  for i in frameshift_mutations ])
                    elif len(indels_mutations) > 0:
                        mutations.extend(indels_mutations)

                    if len(nonsense_mutations) > 0 or len(indels_mutations) > 0 or len(frameshift_mutations) > 0:
                        output_handle.write("\t".join([
                                            element,
                                            contig_id,
                                            str(hsp.query_start),
                                            str(hsp.query_end),
                                            str(hsp.frame[0]),
                                            ','.join(mutations),
                                            str(f"{hsp.identities}/{hsp.align_length} ({identity*100:.2f}%)"), 
                                            '',
                                            '_'.join(hit_names),
                                            alignment.hit_def, 
                                            str(hsp.sbjct_start), 
                                            str(hsp.sbjct_end), 
                                            
                                            str(f"{hsp.sbjct_end - hsp.sbjct_start + 1}/{alignment.length} ({coverage * 100:.2f}%)"), 
                                            ])
                                        +"\n")
    
    if (argv.output and argv.output != "-"):
        output_handle.close()
        logging.info(f'Output save to {argv.output}')
    
    result_handle.close()

    os.unlink(temp_xml.name)

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
  Detect nonsensus/frameshift mutations by running blastx on contig fasta against target protein database
''',
        formatter_class = argparse.RawDescriptionHelpFormatter)
    argv.add_argument('-i', '--input_fasta', dest = 'input_fasta', required = True, help = 'contig fasta file')
    argv.add_argument('-d', '--database', type=str, dest = 'database', required=True, help='protein (diamond) blast DB')
    argv.add_argument('-o', '--output', dest = 'output', required = False,help = 'output tsv file')

    argv.add_argument('--cov', type=float, dest = 'cov', default = .8 , required=False, help='target coverage')
    argv.add_argument('--id',  type=float, dest = 'id',  default = .8 , required=False, help='hit identity')

    argv.add_argument('--indels', action='store_true', help='check in-frame insertion and deletion')
    argv.add_argument('--diamond', action='store_true', help='use diamond blastx')

    argv.add_argument('--threads',  type=float, dest = 'threads',  default = 4 , required=False, help='threads')
    argv.add_argument('--verbose', action='store_true', help='Show more information in log')
    argv.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))

    argvs = argv.parse_args()

    log_level = logging.DEBUG if argvs.verbose else logging.INFO
    logging.basicConfig(
        format='[%(asctime)s' '] %(levelname)s: %(message)s', level=log_level, datefmt='%Y-%m-%d %H:%M')
    
    sys.exit(not main(argvs))
