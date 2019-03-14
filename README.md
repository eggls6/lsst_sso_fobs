# lsst_sso_fobs
Simulate Solar System Object (SSO) observations with LSST

Intended as independent validation to sims_movingObjects (https://github.com/lsst/sims_movingObjects)

needs: 1)LSST observtion campaign database such as kraken_2026.db (exported SQL to ASCII at the moment)
         http://astro-lsst-01.astro.washington.edu:8080/
       
       2) File containing SSO orbits (OpenOrb cometary elements)
       
       3)(optional) JPL planetary ephemeris and spice kernels (e.g. DE431 etc.) for Earth position and orientation
          https://naif.jpl.nasa.gov/naif/toolkit.html
          
 
use: lsst_mba.exe [LSST obs database] [SSO orbit file] [Output file]




