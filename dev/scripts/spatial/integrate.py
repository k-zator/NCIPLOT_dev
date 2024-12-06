#! /usr/bin/env python3

import numpy as np

def integrate_density_cl_cube(cl, X_iso, iso, labels, incr):
    cl_idx = np.where(labels == cl)[0]
    density = np.absolute(0.01*X_iso[cl_idx, 3])
    density_all = np.absolute(0.01*X_iso[:, 3])    
    dvol = incr[0]*incr[1]*incr[2] # and those come from a file and are equal so the incremenrs are assumed to be constant
    # well, should be a fraction of total to cluster_0 number of points - but still needs a reference full intergral value
    # which, btw, is present in the nci_ouput file so should be able to be carried over anyway
    # is the dvol is constant assumption true with the CG2FG function
    return np.sum(density*dvol) 
    

def integrate_ranges(X_iso, incr, outer, inner):
    dvol = incr[0]*incr[1]*incr[2]
    density_sgn = 0.01*X_iso[:, 3]    

    idx_dens_neg = np.where(np.logical_and(density_sgn>-outer, density_sgn<-inner))[0]
    idx_dens_tiny = np.where(np.logical_and(density_sgn>=-inner, density_sgn<=inner))[0]
    idx_dens_pos = np.where(np.logical_and(density_sgn>inner, density_sgn<outer))[0]

    return np.sum(np.absolute(density_sgn)[idx_dens_neg]*dvol), np.sum(np.absolute(density_sgn)[idx_dens_tiny]*dvol), np.sum(np.absolute(density_sgn)[idx_dens_pos]*dvol),

# def integrate_density_cl_cube(cl, X, X_iso, labels, incr):
#     cl_idx = np.where(labels == cl)[0]
#     density = np.absolute(X_iso[cl_idx, 3])
    
#     dvol = incr[0]*incr[1]*incr[2]
#     idx_dens = np.where(density<1000)[0]
#     return np.sum(density[idx_dens]*dvol)
