#! /usr/bin/env python3

import sys
import logging
import time

from spatial.divide import *
from spatial.utils import *
from spatial.plot import *
from spatial.opt_dict import *
from spatial.integrate import * 
#from spatial.plot_new3D import *

# Set options from command line
options = []
if sys.argv[1]!="--help":
    input_name = sys.argv[1]
else:
    options.append(sys.argv[1])

if len(sys.argv)>2:
    options += sys.argv[2:]

opt_dict = options_dict(options)

np.random.seed(opt_dict["seed"])

# Configure warnings
logging.basicConfig(level=logging.WARN, filemode="w", format="  %(levelname)s - %(message)s")

# Print output
print(" # ----------------- NCICLUSTER ------------------------")
print(" # -----------------------------------------------------")
print(" Start -- {} \n".format(time.ctime()))


# Read input file
files = []
with open(input_name, "r") as f:
    for line in f:
        files.append(line[:-1])

# Perform clustering for every file in input
for filename in files:
    print("+ Reading cube files: {}-dens.cube".format(filename))
    print("                      {}-grad.cube".format(filename))
    print(" ")
    header, X, atcoords = pos_dens_grad_matrix(filename)
    incr = np.array([float(header[3].split()[1]), float(header[4].split()[2]), float(header[5].split()[3])])
    X_iso, labels = divide_nci_regions(
        X, opt_dict["n"], isovalue=opt_dict["isovalue"], size_sample=opt_dict["size"], only_pos=opt_dict["onlypos"], method=opt_dict["method"], discard_tails_thr=None #1e-5
    )

    if len(np.unique(labels)) == 1:
        raise ValueError("Only one cluster has been found, you might want to manually set the number of clusters.")
    
    plot_2d(X_iso, labels, filename, X=X, verbose=opt_dict["verbose"], c3=opt_dict["c3"])
    plot_3d(X_iso, labels, filename, verbose=opt_dict["verbose"])
    plot_heatmap_distances(X_iso, labels, filename, verbose=opt_dict["verbose"])
    for cl in set(labels):
        writecube(filename, cl, X_iso, labels, header, verbose=opt_dict["verbose"])
        writedat(filename, cl, X_iso, labels, verbose=opt_dict["verbose"])
        if opt_dict["doint"]:
            print("  Integral of density in cluster {} : {:.8f}".format(cl, integrate_density_cl_cube(cl, X_iso, opt_dict["isovalue"], labels, incr)))
    writevmd(filename, labels, opt_dict["isovalue"], verbose=opt_dict["verbose"], c3=opt_dict["c3"])
    print("\n")
    
    if opt_dict["doint"]:
        int_neg, int_tiny, int_pos = integrate_ranges(X_iso, incr, opt_dict["outer"], opt_dict["inner"])
        print("  Integral of density in ranges :")
        print("  -{} to {} : {:.8f}".format(opt_dict["outer"], opt_dict["inner"], int_neg))
        print("  -{} to  {} : {:.8f}".format(opt_dict["inner"], opt_dict["inner"], int_tiny))
        print("   {} to  {} : {:.8f}".format(opt_dict["inner"], opt_dict["outer"], int_pos))

print(" # -----------------------------------------------------")
print(" End -- {} \n".format(time.ctime()))
