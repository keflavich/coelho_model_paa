import pylab as pl
from astropy import table
from astropy import visualization
pl.ion()


tblfn = 'coelho14_model_paa.fits'
tbl = table.Table.read(tblfn)

cont_est = (tbl['paac_l'] + tbl['paac_h'])/4

teff = tbl['teff']
normcolor = (teff - 3000)/10000
norm = visualization.simple_norm(teff)

pl.figure(1)
pl.clf()
# arbitrary number added to make zero correspond to A0
scat = pl.scatter(-(tbl['paach-paacl'] - 0.47),
                  tbl['cont_m_paa'] / cont_est,
                  c=teff,
                  norm=norm,
                  cmap='jet_r',
                  alpha=0.5,
                 )
pl.colorbar(mappable=scat)
pl.xlabel('Continuum Color ($M_{low} - M_{high}$)')
pl.ylabel('Fractional Line Absorption Depth')
pl.savefig('../paper/figures/paa_cont_colordepth.pdf', bbox_inches='tight')


pl.figure(2)
pl.clf()
# arbitrary numbers added to make zero correspond to A0
scat = pl.scatter(tbl['H-K'] + 1.15,
                  tbl['mag_paa'] - tbl['mK'] + 0.55,
                  c=teff,
                  norm=norm,
                  cmap='jet_r',
                  alpha=0.5,
                 )
pl.colorbar(mappable=scat)
pl.xlabel('H-K')
pl.ylabel('M$_{Pa\\alpha}$-K')
pl.savefig('../paper/figures/paa_hk_colorcolor.pdf', bbox_inches='tight')
