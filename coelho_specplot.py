import numpy as np
from astropy.io import fits
import glob
from spectral_cube import lower_dimensional_structures
from astropy import units as u
import pylab as pl
from astropy.utils.console import ProgressBar

fig = pl.figure(1)
fig.clf()
ax = pl.gca()

for fn in ProgressBar(glob.glob("s_coelho14_sed/*fits")):
    fh = fits.open(fn)
    header = fh[0].header
    sp = lower_dimensional_structures.OneDSpectrum.from_hdu(fh)
    #sp = specutils.Spectrum1D(data=fh[0].data, wcs=wcs.WCS(header), meta={'header': header})
    x = 10**sp.spectral_axis * u.AA

    sel = (x > 18600 * u.AA) & (x < 18900 * u.AA)

    teff = header['TEFF']
    normcolor = (teff - 3000)/10000
    color = pl.cm.jet(normcolor)

    ax.plot(x[sel], sp[sel]/np.mean(sp[sel]), linewidth=0.1, alpha=0.1, color=color)

ax.set_xlabel("Wavelength [$\\circ{A}$]")
ax.set_ylabel("Normalized Spectral Luminosity")
