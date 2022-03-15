#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage: hoxpipe.py                     -p FILE -d FILE [-h|--help]

    -p, --protein_file FILE         FASTA file
    -d, --protein_database FILE     FASTA file
    -h --help                       show this                    
        
"""
from docopt import docopt
import pandas as pd
import collections
from tqdm import tqdm
import pyhmmer
import pyfasta

'''
# braker is overkill and yields general gene annotations
>>> braker.pl --genome hox_genome_fastas/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.clean.fa --prot_seq hox_protein
_fastas/caenorhabditis_elegans.hox_proteins.fa --prg gth --trainFromGth --softmasking --cores 40 --gff3 --workingdir braker_celeg_hox --AUGUSTUS_CONFIG_PATH=/ceph/users/dlaetsch/.conda/envs/braker/config --GENE
MARK_PATH=/ceph/software/genemark/genemark_EX.v4.38/gm_et_linux_64/gmes_petap/ --PROTHINT_PATH=/ceph/users/dlaetsch/software/ProtHint/bin

# prothint is too fiddly to get to work, needs additional steps
>>> ~/software/ProtHint/bin/prothint.py ../hox_genome_fastas/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.clean.fa ../hox_protein
_fastas/caenorhabditis_elegans.hox_proteins.fa --workdir . --threads 20

# metaeuk is the way for de novo annotation of hoxes
# https://github.com/soedinglab/metaeuk
>>> conda install -c conda-forge -c bioconda metaeuk
>>> metaeuk easy-predict ../../hox_genome_fastas/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.clean.fa ../mab5.seed.BM.fasta test temp 
# test was BM mab5 on CELEG, yields partial protein

# pyhmmer is the way to go re HMMs
https://pyhmmer.readthedocs.io
>>> conda install -c bioconda pyhmmer=0.5
'''

'''
[Workflows]
1. single seed based pyhmmer search
2. MSA seed based pyhmmer search
3. Analysis of powerset of all MSA seed based pyhmmer searches
4. Targetted annotation and concerted single seed based pyhmmer search of all target genes
5. Targetted annotation and concerted, snow-balling (?) MSA seed based pyhmmer search of all target genes
'''

#def load_protein_fasta(protein_file):
#    fasta_dict = dict(pyfasta.Fasta(protein_file))
#    print("[+] Parsed %s protein sequences" % len(fasta_dict))
#    return fasta_dict
    
def load_protein_fasta(protein_file):
    with pyhmmer.easel.SequenceFile(protein_file) as sf:
        sequences = list(sf)
        print(sequences[0].name, sequences[0].sequence)
        return sequences

def make_msa(sequences):
    msa = pyhmmer.easel.TextMSA(name=b"msa", sequences=sequences)
    print(msa)
    msa_d = msa.digitize(alphabet)
    return msa_d

def build_hmm(msa_d):
    builder = pyhmmer.plan7.Builder(alphabet)
    print(builder)
    background = pyhmmer.plan7.Background(alphabet)
    print(background)
    hmm, _a, _b = builder.build_msa(msa_d, background)
    print('hmm', hmm)
    print('_a', _a)
    print('_b', _b)
    print('hmm.consensus', hmm.consensus)
    return hmm, background

def save_hmm(hmm, prefix='test'):
    with open("%s.hmm" % prefix, "wb") as output_file:
        hmm.write(output_file)

def parse_protein_dir(protein_dir):
    pass

def align_proteins(proteins):
    pass

#def load_genomes(genome_dir):
#    pass

def hmm_search(hmm, background, database):
    with pyhmmer.easel.SequenceFile(database, digital=True, alphabet=alphabet) as seq_file:
        sequences = list(seq_file)
        pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)
        hits = pipeline.search_hmm(query=hmm, sequences=sequences)
        return hits

if __name__ == '__main__':
    args = docopt(__doc__)
    print('args', args)
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = load_protein_fasta(args['--protein_file'])
    msa_d = make_msa(sequences)
    hmm, background = build_hmm(msa_d)
    hits = hmm_search(hmm, background, args['--protein_database'])
    print('hits', hits)
    print("\n".join(["%s\t%s\t%s\t%s\t%s" % (hit.name, hit.pre_score, hit.pvalue, hit.score, hit.sum_score) for hit in hits]))
    ali = hits[0].domains[0].alignment
    print(" "*3, ali.target_name.decode())
    print("{:3}".format(ali.hmm_from), ali.hmm_sequence[:80] + "...")
    print(" "*3, ali.identity_sequence[:80] + "...")
    print("{:3}".format(ali.target_from), ali.target_sequence[:80] + "...")
    print(" "*3, ali.hmm_name.decode())
