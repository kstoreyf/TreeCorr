#!/usr/bin/env python
import pandas as pd
from astropy.cosmology import LambdaCDM
import numpy as np

import treecorr



nd = 10
data1fn = '../../lss/mangler/samples/a0.6452_0001.v5_ngc_ifield_ndata{}.rdzw'.format(nd)
rand1fn = '../../lss/mangler/samples/a0.6452_rand20x.dr12d_cmass_ngc_ifield_ndata{}.rdz'.format(nd)
data2fn = data1fn
rand2fn = rand1fn

data1 = pd.read_csv(data1fn)
rand1 = pd.read_csv(rand1fn)

cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

def get_comoving_dist(z):
    comov = cosmo.comoving_distance(z)
    return comov.value*cosmo.h

ra = data1['ra'].values
dec = data1['dec'].values
dist = data1['z'].apply(get_comoving_dist).values
print dist

ra_rand = rand1['ra'].values
dec_rand = rand1['dec'].values
dist_rand = rand1['z'].apply(get_comoving_dist).values

ndata = len(ra)
nrand = len(ra_rand)


idx = np.arange(ndata, dtype=long)
idx_rand = np.arange(nrand, dtype=long)


print 'Processing'
cat_data = treecorr.Catalog(ra=ra, dec=dec, r=dist, idx=idx, ra_units='deg', dec_units='deg')
cat_rand = treecorr.Catalog(ra=ra_rand, dec=dec_rand, r=dist_rand, idx=idx_rand, ra_units='deg', dec_units='deg')

dd = treecorr.NNCorrelation(min_sep=10**0, max_sep=10**2, bin_size=1, res_size=ndata**2)
dd.process(cat_data)
print dd.npairs
print dd.idxpairs
print len(dd.idxpairs)
print
dr = treecorr.NNCorrelation(min_sep=10**0, max_sep=10**2, bin_size=1, res_size=ndata*nrand)
dr.process(cat_data, cat_rand)
print dr.npairs
print dr.idxpairs
print len(dr.idxpairs)
print
rr = treecorr.NNCorrelation(min_sep=10**0, max_sep=10**2, bin_size=1, res_size=nrand**2)
rr.process(cat_rand)
print rr.npairs
print rr.idxpairs
print len(rr.idxpairs)

xi, varxi = dd.calculateXi(rr, dr)
print xi
