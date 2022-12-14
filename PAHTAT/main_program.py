from astropy.table import Table
import PAHTAT as pahtat

# Load MIRI stitched spectrum you want to fit
vv114_stitched_spectrum = Table.read('VV114_MIRI_MRS_stitched_spectrum.txt', format='ascii')

# Initialize the class to use pahtat, be sure to add the redshift of the object
stitched_spectrum = pahtat.Spectrum(vv114_stitched_spectrum['values'],
                                    vv114_stitched_spectrum['wavelength'],
                                    resolution=1000,
                                    redshift=0.020067)

# Fit the spectrum and plot results
stitched_spectrum.plot_results()

# Save pahtat results as well as fitted emission lines
stitched_spectrum.save_results(filename='VV114_stitched_MIRI_pahtat_results', gas_lines=True)
