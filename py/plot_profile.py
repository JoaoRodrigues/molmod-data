#!/usr/bin/env python

"""
Produces line plots from MODELLER .profile files.
"""

from __future__ import print_function
import argparse
import os
import sys

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('profile_f', nargs='+', 
                help='MODELLER profile file names')
ap.add_argument('-o', '--output', type=str,
                help='Name to give the resulting plot file')

cmd = ap.parse_args()

import matplotlib
if cmd.output:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

##
xy_data = []

for profile in cmd.profile_f:
    xy_data.append([])
    with open(profile) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            xy_data[-1].append((int(fields[0]), float(fields[41])))

##
f = plt.figure()
ax = plt.gca()

##
res_list = sorted([res[0] for profile in xy_data for res in profile])
min_res = res_list[0]
max_res = res_list[-1]

for profile_data in xy_data:
    res_list, dope_ene = zip(*profile_data)
    ax.plot(res_list, dope_ene)

ax.plot(range(min_res, max_res+1), [-0.03]*(max_res-min_res+1), c='black', linewidth=2.0)

ax.set_xlim(min_res, max_res)
ax.set_xlabel('Residue Number')
ax.set_ylabel('Per-Residue DOPE Energy')
ax.grid('on')

plt.legend(cmd.profile_f, bbox_to_anchor=(0.5, 1.05), ncol=3, loc='upper center')

if cmd.output:
    plt.savefig(cmd.output)
else:
    plt.show()
##

