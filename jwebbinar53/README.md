# JWebbinar 53 workshop notebooks

Start with the three numbered notebooks in this folder. Their code cells are ready for participants to run. Completed notebooks with saved figures and matching PDFs are in [`solutions/`](solutions/).

| Notebook | Topic |
| --- | --- |
| `01_JWebbinar_OCDb_to_ModelReadyData.ipynb` | Turn optical constants into model ready data |
| `02_JWebbinar_WaterClouds.ipynb` | Model a 250 K brown dwarf with and without water clouds, then compare JWST spectra |
| `03_JWebbinar_OCDb_Hazy_SubNeptune.ipynb` | Use OCdb data in a hazy sub Neptune model |

Install [`requirements.txt`](requirements.txt) in a Python 3.11 or 3.12 virtual environment, then start JupyterLab from this folder and run the notebooks in order. The notebooks also work when opened from `solutions/`.

All small inputs used by the exercises are in [`data/`](data/): `ocdb/` holds downloaded optical constants, `observations/` holds spectra, `processed/` holds prepared optical data, and `stellar/` holds the two extracted stellar grid files used by notebook 03. Notebook 02 uses the prepared `H2O.mieff` file in `picaso_reference/virga/`, so participants do not need a separate water-ice OCdb download.

[`picaso_reference/`](picaso_reference/) contains the small PICASO inputs. The large PICASO gas opacity database is kept locally but excluded from Git. On first run, notebooks 02 and 03 download it into `picaso_reference/opacities/` if it is missing. Allow time and several gigabytes of disk space for this first download. The saved notebook 02 solution used the local gas database for both clear and cloudy spectra.

The [`archive/`](archive/) folder contains earlier versions and collaborator material and is excluded from Git. It is not needed for the workshop.
