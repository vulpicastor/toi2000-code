# Data files for TOI-2000

The main data tables are enumerated below.

- `toi_2000_table_02.csv`: TESS light curve of TOI-2000 from year 1 at 30-minute cadence.
- `toi_2000_table_03.csv`: TESS light curve of TOI-2000 from year 3 at 20-second cadence.
- `toi_2000_table_04_B.csv`: LCOGT (Siding Spring Observatory) light curve of TOI-2000, in *B* band.
- `toi_2000_table_04_zs.csv`: LCOGT (Siding Spring Observatory) light curve of TOI-2000, in *z*_s band.
- `toi_2000_table_05_B.csv`: PEST light curve of TOI-2000, in *B* band.
- `toi_2000_table_05_Ic.csv`: PEST light curve of TOI-2000, in *I*_c band.
- `toi_2000_table_06.csv`: ASTEP light curve of TOI-2000.
- `toi_2000_table_08.csv`: Radial velocity measurements of TOI-2000.
- `toi_2000_table_A1.csv`: Multiplanet systems with transiting hot or warm gas
  giants (*P* < 100 d).
- `toi_2000_table_B1.rdb_`: Diagnostics and stellar activity indicators for the
  HARPS spectra of TOI-2000.
- `Gianthosts.csv`: Properties of host stars of multiplanet systems with
  transiting hot or warm gas giants.
- `PS_2022.08.19_15.48.34.csv`: Snapshot of Planetary Systems table from the
  NASA Exoplanet Archive.


## Isochrones

The main notebook `FittingToi2000Xo.ipynb` depends on
[MIST](http://waps.cfa.harvard.edu/MIST/model_grids.html)
stellar evolutionary tracks and bolometric correction grids as repackaged by
the `isochrones`
[Python package](https://github.com/timothydmorton/isochrones)
([Morton et al. 2015](http://ascl.net/1503.010)).
The user is advised to install the `isochrones` package and trigger download
of the full MIST model grids. For the Jupyter notebook in this repository,
the following files must be copied or linked from isochrone's data directory
(by default `~/.isochrones`) to a newly created directory under this repo
`data/isochrones`, and placed in exactly the relative location. The files are:

- `mist/tracks/full_grid_v1.2_vvcrit0.0.npz`
  (If this does not work, `mist/tracks/full_grid_v1.2_vvcrit0.4.npz` can be
  substituted, provided the corresponding filename in `FittingToi2000Xo.ipynb`
  is replaced.)
- `BC/mist/UBVRIplus.h5`
- `BC/mist/WISE.h5`
