#!/usr/bin/env python

from __future__ import print_function, division
import os
import subprocess
import sys
import tempfile

try:
    from Bio.PDB import PDBParser, NeighborSearch
    from Bio.PDB import PDBIO, Select
    from Bio.PDB.Polypeptide import PPBuilder, is_aa
except ImportError as e:
    print('[!] Requires Biopython', file=sys.stderr)
    raise ImportError(e)

# Location of freesasa
FREESASA_BIN='/opt/bin/freesasa'
FREESASA_CONFIG='/opt/share/naccess.config'

# Scaling factors for relative ASA
# Calculated using extended ALA-X-ALA peptides
# Taken from NACCESS
rel_asa = {
    'total':
        {
        'ALA': 107.95,
        'CYS': 134.28,
        'ASP': 140.39,
        'GLU': 172.25,
        'PHE': 199.48,
        'GLY': 80.10,
        'HIS': 182.88,
        'ILE': 175.12,
        'LYS': 200.81,
        'LEU': 178.63,
        'MET': 194.15,
        'ASN': 143.94,
        'PRO': 136.13,
        'GLN': 178.50,
        'ARG': 238.76,
        'SER': 116.50,
        'THR': 139.27,
        'VAL': 151.44,
        'TRP': 249.36,
        'TYR': 212.76,
        },
    'bb':
        {
        'ALA': 38.54,
        'CYS': 37.53,
        'ASP': 37.70,
        'GLU': 37.51,
        'PHE': 35.37,
        'GLY': 47.77,
        'HIS': 35.80,
        'ILE': 37.16,
        'LYS': 37.51,
        'LEU': 37.51,
        'MET': 37.51,
        'ASN': 37.70,
        'PRO': 16.23,
        'GLN': 37.51,
        'ARG': 37.51,
        'SER': 38.40,
        'THR': 37.57,
        'VAL': 37.16,
        'TRP': 38.10,
        'TYR': 35.38,
        },
    'sc':
        {
        'ALA': 69.41,
        'CYS': 96.75,
        'ASP': 102.69,
        'GLU': 134.74,
        'PHE': 164.11,
        'GLY': 32.33,
        'HIS': 147.08,
        'ILE': 137.96,
        'LYS': 163.30,
        'LEU': 141.12,
        'MET': 156.64,
        'ASN': 106.24,
        'PRO': 119.90,
        'GLN': 140.99,
        'ARG': 201.25,
        'SER': 78.11,
        'THR': 101.70,
        'VAL': 114.28,
        'TRP': 211.26,
        'TYR': 177.38,
    }
}

def parse_structure(path):
    """
    Parses a PDB formatter structure using Biopython's PDB Parser
    Verifies the integrity of the structure (gaps) and its
    suitability for the calculation (is it a complex?).
    """

    print('[+] Reading structure file: {0}'.format(path))
    fname = os.path.basename(path)
    sname = '.'.join(fname.split('.')[:-1])

    try:
        s = P.get_structure(sname, path)
    except Exception as e:
        print('[!] Structure \'{0}\' could not be parsed'.format(sname), file=sys.stderr)
        raise Exception(e)

    # Double occupancy check
    for atom in list(s.get_atoms()):
        if atom.is_disordered():
            residue = atom.parent
            sel_at = atom.selected_child
            sel_at.altloc = ' '
            sel_at.disordered_flag = 0
            residue.detach_child(atom.id)
            residue.add(sel_at)

    # Remove HETATMs and solvent
    res_list = list(s.get_residues())
    n_res = len(res_list)
    _ignore = lambda r: r.id[0][0] == 'W' or r.id[0][0] == 'H'
    for res in res_list:
        if _ignore(res):
            chain = res.parent
            chain.detach_child(res.id)
        elif not is_aa(res, standard=True):
            raise ValueError('Unsupported non-standard amino acid found: {0}'.format(res.resname))

    # Detect gaps and compare with no. of chains
    pep_builder = PPBuilder()
    peptides = pep_builder.build_peptides(s)
    n_peptides = len(peptides)
    n_chains = len(set([c.id for c in s.get_chains()]))

    if n_peptides != n_chains:
        print('[!] Structure contains gaps:', file=sys.stderr)
        for i_pp, pp in enumerate(peptides):
            print('\t{1.parent.id} {1.resname}{1.id[1]} < Fragment {0} > {2.parent.id} {2.resname}{2.id[1]}'.format(i_pp, pp[0], pp[-1]), file=sys.stderr)
        #raise Exception('Calculation cannot proceed')

    return (s, n_chains, n_res)

def parse_freesasa_output(fpath):
    """
    Returns per-residue relative accessibility of side-chain and main-chain
    atoms as calculated by freesasa.
    """

    rsa_data = {}

    _rsa = rel_asa
    _bb = set(('CA', 'C', 'N', 'O'))

    s = P.get_structure('bogus', fpath.name)
    for res in s.get_residues():
        res_id = (res.parent.id, res.resname, res.id[1])
        asa_bb, asa_sc, total_asa = 0, 0, 0
        for atom in res:
            asa = atom.bfactor
            if atom.name in _bb:
                asa_bb += asa
            else:
                asa_sc += asa
            total_asa += asa

        rsa_data[res_id] = (100*total_asa / _rsa['total'][res.resname],
                            100*asa_bb / _rsa['bb'][res.resname],
                            100*asa_sc / _rsa['sc'][res.resname])

    return rsa_data

def execute_freesasa(structure, pdb_selection=None):
    """
    Runs the freesasa executable on a PDB file.

    You can get the executable from:
        https://github.com/mittinatten/freesasa

    The binding affinity models are calibrated with the parameter
    set for vdW radii used in NACCESS:
        http://www.ncbi.nlm.nih.gov/pubmed/994183
    """

    freesasa, param_f= FREESASA_BIN, FREESASA_CONFIG
    if not os.path.isfile(freesasa):
        raise IOError('[!] freesasa binary not found at `{0}`'.format(freesasa))
    if not os.path.isfile(param_f):
        raise IOError('[!] Atomic radii file not found at `{0}`'.format(param_f))

    # Rewrite PDB using Biopython to have a proper format
    # freesasa is very picky with line width (80 characters or fails!)
    # Select chains if necessary
    class ChainSelector(Select):
        def accept_chain(self, chain):
            if pdb_selection and chain.id in pdb_selection:
                return 1
            elif not pdb_selection:
                return 1
            else:
                return 0

    _pdbf = tempfile.NamedTemporaryFile()
    io.set_structure(structure)
    io.save(_pdbf.name, ChainSelector())

    # Run freesasa
    # Save atomic asa output to another temp file
    _outf = tempfile.NamedTemporaryFile()
    cmd = '{0} -B {1} -L -d 0.05 -c {2} {3}'.format(freesasa, _outf.name, param_f, _pdbf.name)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    if p.returncode:
        print('[!] freesasa did not run successfully', file=sys.stderr)
        raise Exception(stderr)

    # Rewind & Parse results file
    # Save
    _outf.seek(0)
    rsa = parse_freesasa_output(_outf)

    _pdbf.close()
    _outf.close()

    return rsa


if __name__ == '__main__':
    P = PDBParser(QUIET=1)
    io = PDBIO()

    # Parse structure
    pdb_path = sys.argv[1]
    structure, n_chains, n_res = parse_structure(pdb_path)
    print('[+] Parsed PDB file {0} ({1} chains, {2} residues)'.format(structure.id, n_chains, n_res))
    
    cmplx_sasa = execute_freesasa(structure)
    print('# Residue\tbbRSA\tscRSA')
    for res in sorted(cmplx_sasa, key=lambda x: x[2]):
        print('{0[1]} {0[2]:<4d}\t{1[1]:>6.2f}\t{1[2]:>6.2f}'.format(res, cmplx_sasa[res]))
