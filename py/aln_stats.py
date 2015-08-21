#!/usr/bin/env python

"""
Utility to extract information from a series of pairwise sequence
alignments produced by HMMER. Additionally, calculates a number of
other metrics useful for homology modelling.

Written by {0} [{1}]
"""

from __future__ import print_function, division

import os
import sys
import warnings

try:
    from Bio import BiopythonExperimentalWarning
    warnings.filterwarnings('ignore', category=BiopythonExperimentalWarning)
    from Bio import SearchIO
except ImportError as exception:
    print("[!] Could not import Biopython modules", file=sys.stderr)
    raise exception

__author__ = "Joao Rodrigues"
__email__ = "j.p.g.l.m.rodrigues@gmail.com"

#

def _check_file(apath):
    """
    Validates the existence and permissions of a given path.
    Returns the full absolute path if successful.
    """
    full_path = os.path.abspath(apath)
    if os.path.exists(full_path):
        return full_path
    else:
        raise IOError('File ({0}) could not be found or read'.format(apath))

def read_hmmer_file(f_path):
    """
    Uses Biopython's SearchIO to parse a HMMER output file.
    Returns an iterator with search hits.
    """

    f_path = _check_file(f_path)
    return SearchIO.read(f_path, 'hmmer3-text')

def get_output_stream(stream):
    """
    Returns a file-like object to output things to.
    """
    if cmd.output == sys.stdout:
        return cmd.output
    else:
        return open(stream, 'w')

def calculate_seq_id(seqA, seqB):
    """
    Calculates the pairwise sequence identity between two aligned sequences.
    Returns the sequence identity (fraction).
    """

    size_of_A = len(seqA)
    assert len(seqB) == size_of_A, \
        'Sequence length does not match: {0} vs. {1}'.format(size_of_A, len(seqB))

    identical = sum([True for (aa_A, aa_B) in zip(seqA, seqB) if aa_A == aa_B])
    return (identical/size_of_A)

#
if __name__ == '__main__':

    import argparse
    ap = argparse.ArgumentParser(description=__doc__.format(__author__, __email__))
    ap.add_argument('alignment', type=str, help='HMMER output file')
    ap.add_argument('-o', '--output', type=str, default=sys.stdout,
                    help='Output file [default: stdout]')
    cmd = ap.parse_args()

    hits_iter = read_hmmer_file(cmd.alignment)
    ostream = get_output_stream(cmd.output)

    # PDBID, Evalue, BitScore, %Id, %Coverage, Description
    header = ('PDBID', 'E-value', 'Bit Score', 'Seq. Id.', 'Seq. Cov.', 'Hit Description')
    oheader = '#{0:6s}\t{1:>10s}\t{2:>9s}\t{3:>8s}\t{4:>8s}\t{5}'
    oformat = '{0:6s}\t{1:10.1E}\t{2:9.2f}\t{3:8.2f}\t{4:8.2f}\t{5}'

    with ostream as out:
        print(oheader.format(*header), file=out)
        q_full_len = hits_iter.seq_len
        for hit in hits_iter:
            for hsp in hit.hsps:
                aln_id = calculate_seq_id(hsp.query.seq.upper(), hsp.hit.seq)
                aln_cov = 1 - ((q_full_len - hsp.hit_span)/q_full_len)
                hit_description = ' '.join(hsp.hit_description.split()[2:]) # pdb_seqres
                print(oformat.format(hit.id, hit.evalue, hit.bitscore,
                                     aln_id, aln_cov, hit_description), file=out)
