#!/usr/bin/env python

"""
Python script to run homology modelling calculations
using MODELLER without much hassle.

JR @ Utrecht, 2015
"""

from __future__ import print_function

import argparse
import os
import sys

try:
    from modeller import *                       
    from modeller.automodel import *             
    from modeller.parallel import *              
except ImportError as e:
    print('[!] Could not import MODELLER. Check your installation and environment', file=sys.stderr)
    print(e, file=sys.stderr)
    sys.exit(1)

# Argument Handling
ap = argparse.ArgumentParser(description=__doc__)

io = ap.add_argument_group('Input options')
io.add_argument('-a', '--alignment', type=str, required=True,
                help='Alignment file in PIR format')
io.add_argument('-t', '--template', type=str, nargs='+', required=True,
                help='Templates to use in PDB format')

ro = ap.add_argument_group('Runtime Options')
ro.add_argument('--serial', action='store_true',
                help='Performs the calculation using a single process')
ro.add_argument('--num_jobs', type=int,
                help='Number of processes to spawn during modelling.')

mo = ap.add_argument_group('Modelling Options')
mo.add_argument('--num_models', type=int, default=10,
                help='Number of models to create.')

mo.add_argument('--loopmodel', action='store_true',
                help='Uses the loop modelling routines of MODELLER')
mo.add_argument('--num_loop_models', type=int, default=10,
                help='Number of loop models to create per model')

mo.add_argument('--use_dope', action='store_true',
                help='Uses the DOPE function to score the models')

mo.add_argument('--cluster', action='store_true',
                help='Clusters the generated models by RMSD')
mo.add_argument('--cluster_cutoff', type=float, default=1.5,
                help='Clustering cutoff in Angstrom [default=1.5]')

cmd = ap.parse_args()

### Check input arguments
templates = []
for template_f in cmd.template:
    if os.path.isfile(template_f):
        sname = '.'.join(template_f.split('.')[:-1])
        templates.append(sname)
    else:
        print('[!] Could not open template file: {0}'.format(template_f), file=sys.stderr)
        sys.exit(1)
templates = tuple(templates)

if os.path.isfile(cmd.alignment):
    ali_templates = set()
    with open(cmd.alignment) as handle:
        for lineno, line in enumerate(handle, start=1):
            if line.startswith('>'):
                entry_code = line.split(';')[1].strip()
            elif line.startswith('sequence'):
                entry_code2 = line.split(':')[1].strip()
                if entry_code != entry_code2:
                    print('[!] Mismatch in your alignment at line {0}: {1!r} =/= {2!r})'.format(lineno, entry_code, entry_code2), file=sys.stderr)
                    sys.exit(1)
                seq_name = entry_code2
            elif line.startswith('structure'):
                entry_code2 = line.split(':')[1].strip()
                if entry_code != entry_code2:
                    print('[!] Mismatch in your alignment at line {0}: {1!r} =/= {2!r})'.format(lineno, entry_code, entry_code2), file=sys.stderr)
                    sys.exit(1)
                if entry_code2 not in templates:
                    print('[!] Alignment template missing from your input: {0}'.format(entry_code2), file=sys.stderr)
                    sys.exit(1)
                ali_templates.add(entry_code2)
    common = set(templates) - ali_templates
    if common:
        common = ','.join(sorted(common))
        print('[!] Your alignment contains more templates than specified: {0}'.format(common), file=sys.stderr)
        sys.exit(1)
else:
    print('[!] Could not open alignment file: {0}'.format(cmd.alignment), file=sys.stderr)
    sys.exit(1)

### Run Modelling
## Setup MODELLER environment
j = job()
j.append(local_slave())
if not cmd.serial and cmd.num_jobs:
    for req_job in range(cmd.num_jobs - 1):
        j.append(local_slave())

log.minimal()                                
env = environ()                              
env.io.atom_files_directory = [os.curdir]
env.io.hetatm = True

## Model Building
if cmd.use_dope:
    assess_methods = (assess.DOPE,)
else:
    assess_methods = ()

if cmd.loopmodel:
    a = loopmodel(env,
                  alnfile  = cmd.alignment,	
                  knowns   = templates,		
                  sequence = seq_name,		
                  assess_methods=assess_methods,	
                  loop_assess_methods=assess_methods,
                  )
else:
    a = automodel(env,
                  alnfile  = cmd.alignment,	
                  knowns   = templates,		
                  sequence = seq_name,		
                  assess_methods=assess_methods,	
                  )

## Number of Models to generate
# Backbone models
a.starting_model= 1				
a.ending_model  = cmd.num_models				
# Loop models (refined, if required)
if cmd.loopmodel:
    a.loop.starting_model = 1           
    a.loop.ending_model   = cmd.num_loop_models          
    a.loop.md_level       = refine.slow 

## Optimization Settings
# Conjugate Gradient
a.library_schedule = autosched.slow
a.max_var_iterations = 300			
a.max_molpdf = 1e6				
#MD
a.md_level = refine.slow			

## Finalize settings and launch modelling
a.use_parallel_job(j)                        
a.make()                                     

## Cluster the resulting models
if cmd.cluster:
    a.cluster(cluster_cut=cmd.cluster_cutoff)
