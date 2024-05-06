#! /usr/bin/env python3

import numpy as np

def integrate_density_cl_cube(cl, X_iso, iso, labels, incr):
    cl_idx = np.where(labels == cl)[0]
    density = np.absolute(0.01*X_iso[cl_idx, 3])
    density_all = np.absolute(0.01*X_iso[:, 3])    
    dvol = incr[0]*incr[1]*incr[2]
    return np.sum(density*dvol)#/np.sum(density_all*dvol)

def integrate_ranges(X_iso, incr):
    dvol = incr[0]*incr[1]*incr[2]
    density_sgn = 0.01*X_iso[:, 3]    

    idx_dens_neg = np.where(np.logical_and(density_sgn>-0.07, density_sgn<-0.01))[0]
    idx_dens_tiny = np.where(np.logical_and(density_sgn>=-0.01, density_sgn<=0.01))[0]
    idx_dens_pos = np.where(np.logical_and(density_sgn>0.01, density_sgn<0.07))[0]

    return np.sum(np.absolute(density_sgn)[idx_dens_neg]*dvol), np.sum(np.absolute(density_sgn)[idx_dens_tiny]*dvol), np.sum(np.absolute(density_sgn)[idx_dens_pos]*dvol),

# def integrate_density_cl_cube(cl, X, X_iso, labels, incr):
#     cl_idx = np.where(labels == cl)[0]
#     density = np.absolute(X_iso[cl_idx, 3])
    
#     dvol = incr[0]*incr[1]*incr[2]
#     idx_dens = np.where(density<1000)[0]
#     return np.sum(density[idx_dens]*dvol)
