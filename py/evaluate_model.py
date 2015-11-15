#!/usr/bin/env python

"""
Produces a per-residue energy evaluation of the model(s).
"""

from __future__ import print_function
import argparse
import os
import sys

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('pdb_f', nargs='+',
                help='PDB file names')

cmd = ap.parse_args()

from modeller import *
from modeller.scripts import complete_pdb

env = environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

for pdb in cmd.pdb_f:
    # read model file
    mdl = complete_pdb(env, pdb)

    # Assess with DOPE:
    s = selection(mdl)   # all atom selection
    s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=pdb[:-4]+'.dope_profile',
                  normalize_profile=True, smoothing_window=15)
