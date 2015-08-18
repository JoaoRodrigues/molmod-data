#!/usr/bin/env python

"""
Utility to match (and compare) PDB files.

WILL NOT WORK PROPERLY ON HOMO-MULTIMERS BECAUSE OF THE CHAIN
MATCHING LOOP.

Uses global sequence alignment to find equivalent positions
between the sequences. Also superimposes the structures based
on the alignments and outputs per-chain RMSDs.

Outputs several values, see the example:
         |------------------------------- Matched Chains (Ref<>Mobi)
         |
         |     |---------------------------- Full alignment seq. id
         |     |
         |     |        |------------------ Gapless seq. id (aligment without edge gaps)
         |     |        |
         |     |        |         |------- Chain RMSD (n. atoms used in RMS/total in chain)
[++] A<>A :=  95.4% |  99.0% |  1.70A
	      + --------------------------------------------------------- Sequence mismatch
	MSVPTDGAVTTSQIPASEQETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRL --- Reference Sequence (trimmed)
	MSVPTDEAVTTSQIPASEQETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRL --- Mobile Sequence (trimmed)
	***************   *   ------------------------------------------- Removed during RMS calculation

Written by {0} [{1}]
"""

from __future__ import print_function, division

from operator import itemgetter
import os
import sys
import tempfile
import warnings

try:
    from Bio import BiopythonExperimentalWarning
    warnings.filterwarnings('ignore', category=BiopythonExperimentalWarning)

    from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

except ImportError as exception:
    print("[!] Could not import Biopython modules", file=sys.stderr)
    raise exception

__author__ = "Joao Rodrigues"
__email__ = "j.p.g.l.m.rodrigues@gmail.com"

#
def parse_structure(struct):
    """
    Parses a PDB-formatted structure using the PDBParser module
    of Biopython. Filters HETATMs and Waters. Returns a SMCRA object.
    """

    sname = '.'.join(os.path.basename(struct).split('.')[:-1])
    full_path = os.path.abspath(struct)
    if not os.path.isfile(full_path):
        print('[!!] File not found or not readable'.format(struct), file=sys.stderr)
        sys.exit(1)

    s = parser.get_structure(sname, full_path)
    # Ensemble check
    if len(s) > 1:
        print('[!!] Ensembles are not supported: {0}'.format(struct), file=sys.stderr)
        print('[!!] Using first model only')

        class ModelSelector(Select):
            def accept_model(self, model):
                if model.id == 1:
                    return 1
                else:
                    return 0

        tempf = tempfile.NamedTemporaryFile()
        io.set_structure(s)
        io.save(tempf.name, ModelSelector())

        s = parser.get_structure(sname, tempf)
        tempf.close()
        
    # Double occupancy check
    for atom in list(s.get_atoms()):
        if atom.is_disordered():
            residue = atom.parent
            sel_at = atom.selected_child
            sel_at.altloc = ' '
            sel_at.disordered_flag = 0
            residue.detach_child(atom.id)
            residue.add(sel_at)

    # Remove HETATM and Solvent
    res_list = list(s.get_residues())
    n_res = len(res_list)
    _ignore = lambda r: r.id[0][0] == 'W' or r.id[0][0] == 'H'
    for res in res_list:
        if _ignore(res):
            chain = res.parent
            chain.detach_child(res.id)
    nn_res = len(list(s.get_residues()))

    # print('[+] Parsed PDB file {0} ({1}/{2} residues kept)'.format(struct, nn_res, n_res))
    return s

def _align_sequences(structA, structB, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """

    def _calculate_identity(sequenceA, sequenceB):
        """
        Returns the percentage of identical characters between two sequences.
        Assumes the sequences are aligned.
        """

        sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
        matches = [sa[i] == sb[i] for i in xrange(sl)]
        seq_id = (100 * sum(matches)) / sl

        gapless_sl = sum([1 for i in xrange(sl) if (sa[i] != '-' and sb[i] != '-')])
        gap_id = (100 * sum(matches)) / gapless_sl
        return (seq_id, gap_id)

    def _get_pdb_sequence(structure):
        """
        Retrieves the AA sequence from a PDB structure.
        """

        _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
        seq = map(_aainfo, structure.get_residues())
        return seq

    #
    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    trim_ends = kwargs.get('trim_ends', True)

    resseq_A = _get_pdb_sequence(structA)
    resseq_B = _get_pdb_sequence(structB)

    sequence_A = ''.join(map(itemgetter(1), resseq_A))
    sequence_B = ''.join(map(itemgetter(1), resseq_B))
    alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                    matrix, gap_open, gap_extend,
                                    penalize_end_gaps=(False, False) )

    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln

    # Equivalent residue numbering
    # Relative to reference
    mapping = {}
    aa_i_A, aa_i_B = 0, 0
    for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
        if aa_aln_A == '-':
            if aa_aln_B != '-':
                aa_i_B += 1
        elif aa_aln_B == '-':
            if aa_aln_A != '-':
                aa_i_A += 1
        else:
            assert resseq_A[aa_i_A][1] == aa_aln_A
            assert resseq_B[aa_i_B][1] == aa_aln_B
            mapping[resseq_A[aa_i_A][0]] = resseq_B[aa_i_B][0]
            aa_i_A += 1
            aa_i_B += 1

    # Gapless alignment
    def _trimmer(sequence):
        """Returns indices of first and last ungapped position"""

        leading = [i for (i, aa) in enumerate(sequence) if aa != '-'][0]
        trailing = [i for (i, aa) in enumerate(sequence[::-1]) if aa != '-'][0]

        trailing = len(sequence) - trailing
        return (leading, trailing)

    lead_A, trail_A = _trimmer(aligned_A)
    lead_B, trail_B = _trimmer(aligned_B)

    lead = max(lead_A, lead_B)
    trail = min(trail_A, trail_B)
    trim_aln_A = aligned_A[lead:trail]
    trim_aln_B = aligned_B[lead:trail]
    mismatch = ''.join(['+' if a!=b else ' ' for (a,b) in zip(trim_aln_A, trim_aln_B)])

    # Calculate (gapless) sequence identity
    seq_id, g_seq_id = _calculate_identity(aligned_A, aligned_B)
    return ((trim_aln_A, trim_aln_B, mismatch), seq_id, g_seq_id, mapping)

def match_chains(reference, mobile, min_id=30.0):
    """
    Matches the chains of two different structures using
    pairwise sequence alignment. Minimum sequence id is 30%
    and minimum gapless is 90% (fixed).
    """
    r_chains = list(reference.get_chains())
    m_chains = list(mobile.get_chains())

    chain_map = {}
    _mapped = set()
    for rc in r_chains:
        best_id = min_id
        for mc in m_chains:
            # Naive check to avoid remapping chains
            if mc.id in _mapped:
                continue

            aln, seq_id, gap_id, res_map = _align_sequences(rc, mc)
            if seq_id >= best_id and gap_id >= 90.0:
                _mapped.add(mc.id)
                chain_map[rc.id] = (mc.id, seq_id, gap_id, aln, res_map)
                best_id = seq_id
                if seq_id == 100.0:
                    break

        # if rc.id not in chain_map:
        #     print('[+++] {0} = No Match'.format(rc.id))
        # else:
        #     mc_id, seq_id, gap_id, _, _ = chain_map[rc.id]
        #     print('[+++] {0}<>{1} = {2:5.1f} | {3:5.1f}'.format(rc.id, mc_id, seq_id, gap_id))

    return chain_map

def match_structures(reference, mobile, mapping):
    """
    Outputs fitted structures of 'trimmed' mobile structures,
    renumbered to match the reference.
    RMSD calculated only on matching residues.
    """

    # Tag mobile chains
    for chain in list(mobile.get_chains()):
        ori = chain.id
        new = "{0}_tag".format(chain.id)
        mobile[0].child_dict[new] = chain
        del mobile[0].child_dict[ori]
        chain.id = new

    # Rename matching residues/chains
    for rc_id in mapping:
        mc_id, _, _, _, c_map = mapping[rc_id]
        mc = mobile[0]['{0}_tag'.format(mc_id)]

        inv_c_map = dict([(v, k) for k, v in c_map.items()])
        mc_res_list = list(mc.get_residues())
        for res in mc_res_list:
            if res.id[1] not in inv_c_map:
                mc.detach_child(res.id)
            else:
                res.id = (res.id[0], inv_c_map[res.id[1]], res.id[2])

        mc.id = rc_id
    # Remove leftover chains
    for chain in list(mobile.get_chains()):
        if chain.id.endswith('_tag'):
            mobile[0].detach_child(chain.id)

    # Superimposer (iterative, discards atoms)
    def _dist_per_atom(atoms_A, atoms_B):
        """
        Calculates pairwise distances of atoms.
        Returns indices of outliers (>av).
        """
        dist = [i-j for (i,j) in zip(atoms_A, atoms_B)]
        av = sum(dist) / len(dist)
        outliers = [True if d > av else False for d in dist]
        return outliers

    n_trials = 0
    r_atoms = [r['CA'] for r in reference.get_residues() if r.parent.id in mapping and r.id[1] in mapping[r.parent.id][4]]
    m_atoms = [r['CA'] for r in mobile.get_residues()]
    ori_n_atoms = len(r_atoms)

    while 1:
        n_trials += 1

        super_imposer.set_atoms(r_atoms, m_atoms)
        super_imposer.apply(mobile.get_atoms())

        n_atoms = len(r_atoms)
        outliers = _dist_per_atom(r_atoms, m_atoms)
        n_outliers = sum(outliers)
        # print(n_trials, super_imposer.rms, n_outliers, n_atoms)

        # Stop if less than 90% of starting atoms or less than 5 outliers OR after 3 trials
        if n_atoms < ori_n_atoms*0.9 or n_outliers < 5 or n_trials == 3:
            break

        # Remove outliers for better RMS
        r_atoms = [r for (r, o) in zip(r_atoms, outliers) if not o]
        m_atoms = [r for (r, o) in zip(m_atoms, outliers) if not o]

    # Build full list of outliers
    _r_atoms = set(r_atoms)
    r_atoms = [r['CA'] for r in reference.get_residues() if r.parent.id in mapping and r.id[1] in mapping[r.parent.id][4]]
    m_atoms = [r['CA'] for r in mobile.get_residues()]
    outliers = [r not in _r_atoms for r in r_atoms]

    # Mark highest deviations in the sequence
    # Calculate rms per chain (excl. outliers)
    # chain: (distances, num outliers)
    rms_info = dict([(c, [[], []]) for c in mapping.keys()])

    for r_ca, m_ca, is_outlier in zip(r_atoms, m_atoms, outliers):
        r_chain = r_ca.parent.parent.id
        if is_outlier:
            distance = 0
        else:
            distance = r_ca - m_ca
        rms_info[r_chain][0].append(distance)
        rms_info[r_chain][1].append(is_outlier)

    for r_chain in rms_info:
        chain_distances, c_outliers = rms_info[r_chain]
        n_residues = len(chain_distances) - sum(c_outliers)
        rms = sum(map(lambda x: x**2, chain_distances)) / n_residues
        rms = rms**0.5
        markers = ''.join(['*' if o else ' ' for o in c_outliers])
        rms_info[r_chain] = (rms, n_residues, markers)

    return rms_info

def print_summary(reference, mobile, mapping, rmsd_info):
    """
    """

    for ref_chain in reference.get_chains():
        r_id = ref_chain.id
        if r_id in mapping:
            mc_id, seq_id, gap_id, aln, res_map = mapping[r_id]
            c_len = len(res_map)
            chain_rms, n_res, markers = rmsd_info[r_id]
            _vals = (r_id, mc_id, seq_id, gap_id, chain_rms, n_res, c_len)
            print('[++] {0}<>{1} := {2:5.1f}% | {3:5.1f}% | {4:5.2f}A ({5}/{6})'.format(*_vals))
            print('\t{0[2]}\n\t{0[0]}\n\t{0[1]}\n\t{1}'.format(aln, markers))
        else:
            print('[++] {0}<>None := No Match'.format(r_id))

def save_mobile(reference, mobile):
    # Write matched/fit file
    io.set_structure(mobile)
    match_name = '{0}_matched_to_{1}.pdb'.format(mobile.id, reference.id)
    print('[==] Saved {0}'.format(match_name))
    io.save(match_name)

#
if __name__ == '__main__':

    import argparse
    from argparse import RawTextHelpFormatter

    ap = argparse.ArgumentParser(description=__doc__.format(__author__, __email__), formatter_class=RawTextHelpFormatter)

    ap.add_argument('pdbf_list', type=str, nargs='+', help='PDB files')
    ap.add_argument('--reference', type=str, help='Reference PDB file for comparison')
    cmd = ap.parse_args()

    # Bio.PDB classes
    parser = PDBParser(QUIET=1)
    io = PDBIO()
    super_imposer = Superimposer()

    # The Real Stuff
    # Read reference first
    refe_path = cmd.reference if cmd.reference else cmd.pdbf_list[0]
    print('[+] Matching structures to {0}'.format(refe_path))
    reference = parse_structure(refe_path)

    # Iterate over others
    for pdbf in cmd.pdbf_list:
        mobile = parse_structure(pdbf)
        print('[+] Comparing structures: {0} vs {1}'.format(reference.id, mobile.id))
        chain_mapping = match_chains(reference, mobile)
        rmsd_info = match_structures(reference, mobile, chain_mapping)
        print_summary(reference, mobile, chain_mapping, rmsd_info)
        save_mobile(reference, mobile)
