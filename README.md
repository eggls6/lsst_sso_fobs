# openobs
Open source Solar System Object (SSO) observations simulation software

Orininally developed as independent validation software for the LSST project's sims_movingObjects (https://github.com/lsst/sims_movingObjects) package, openobs has evolved into a stand-alone open source SSO observation simulator. 

needs: 1) Simulated observation campaign database (for LSST e.g. kraken_2026.db (exported SQL to ASCII at the moment)
         (http://astro-lsst-01.astro.washington.edu:8080/)
       2) File containing SSO orbits (OpenOrb cometary elements)
       3)(optional) JPL planetary ephemeris and spice kernels (e.g. DE431 etc.) for Earth position and orientation
         (https://naif.jpl.nasa.gov/naif/toolkit.html)
          

## How to?

lsst_mba [LSST obs database] [SSO orbit file] [Output file]





