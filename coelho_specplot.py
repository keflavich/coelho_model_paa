import numpy as np
from astropy.io import fits
import glob
from spectral_cube import lower_dimensional_structures
from astropy import units as u
import pylab as pl
from astropy.utils.console import ProgressBar
from astropy import visualization

teffs = []
data = []

for fn in ProgressBar(glob.glob("s_coelho14_sed/*fits")):
    fh = fits.open(fn)
    header = fh[0].header
    sp = lower_dimensional_structures.OneDSpectrum.from_hdu(fh)
    #sp = specutils.Spectrum1D(data=fh[0].data, wcs=wcs.WCS(header), meta={'header': header})
    x = 10**sp.spectral_axis * u.AA

    sel = (x > 18400 * u.AA) & (x < 19100 * u.AA)

    teff = header['TEFF']
    normcolor = (teff - 3000)/10000
    color = pl.cm.jet(normcolor)

    teffs.append(teff)
    data.append(sp[sel])
    xsel = x[sel]

data = np.array(data)
ndata = data / data[:,-1:]
newx = np.linspace(xsel.min(), xsel.max(), 200)
from scipy.interpolate import interp1d
ndata = interp1d(xsel, ndata, kind='cubic')(newx)

norm = visualization.simple_norm(teffs)

segments = np.array([list(zip(newx.value,d)) for d in ndata])


fig = pl.figure(1)
fig.clf()
ax = pl.gca()

lines = pl.matplotlib.collections.LineCollection(segments=segments,
                                                 cmap='jet_r',
                                                 alpha=0.05,
                                                 norm=norm)
lines.set_array(np.array(teffs))
ax.add_collection(lines)

transmission_curve_lower = (0.95 + np.random.randn(newx.size)/1000) * ((newx > 18610*u.AA) & (newx < 18710*u.AA))
transmission_curve_paa = (0.95 + np.random.randn(newx.size)/1000) * ((newx > 18731*u.AA) & (newx < 18781*u.AA))
transmission_curve_upper = (0.95 + np.random.randn(newx.size)/1000) * ((newx > 18800*u.AA) & (newx < 18900*u.AA))

ax.plot(newx, np.array([transmission_curve_lower, transmission_curve_paa, transmission_curve_upper]).T, color='k')

#ax.plot(xsel, ndata, linewidth=0.1, alpha=0.1, color=color)

ax.set_xlim(18500, 19000)
ax.set_ylim(0.9, 1.1)
ax.set_xlabel("Wavelength [$\\AA$]")
ax.set_ylabel("Normalized Spectral Luminosity")
cb = pl.colorbar(mappable=lines)
cb.set_alpha(1)
cb.draw_all()
cb.set_label('Effective Temperature [K]')

pl.savefig("model_stellar_spectra.png", bbox_inches='tight')
pl.savefig("model_stellar_spectra.pdf", bbox_inches='tight')
