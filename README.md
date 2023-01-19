MetReader
==========

MetReader is a library written in fortran 90 that provides an interface to
numerical weather prediction (NWP) data, or other forms of meteorological
data, such as radiosonde or other 1-d data.  

This library was originally written as a component of the USGS volcanic ash
transport and dispersion model, Ash3d.  However, since it is useful for
other programs other than Ash3d, this interface to NWP files is provided
as a separate repository that can either be compiled as a library or simply
compiled directly with other source code.

NWP data are generally made available by agencies (NCEP, NOAA, NASA, etc.)
in a variety of formats (NetCDF, GRIB, ASCII); each product having
its own data structure, naming convention, units, etc.  This library 
isolates the calling program from the peculiarities of interfacing with
a particular NWP product.  Data can be returned to the calling program on
the native grid of the NWP product, or on any grid needed by the calling
program.  Projection and interpolation of NWP data to the required grid, along
with any rotation of velocity vectors to grid-relative, is calculated internally by 
MetReader.

For details on usage, please see the User's Guide and look through the example
programs.

This library requires two additional libraries made available on GitHub and USGS GitLab:

- [HoursSince](https://github.com/DOI-USGS/volcano-ash3d-hourssince)
- [projection](https://github.com/DOI-USGS/volcano-ash3d-projection)

Additionally, the default makefile will build MetReader with both NetCDF and GRIB2
enabled.  If either of these libraries are unavailable on your system, you can
deactivate those options by setting the corresponding flags to 'F' in the makefile.

To compile as a library, simple type:

  `make all`

This will build the requested components of the library.  If GRIB2 is enabled, 
it is recommended to also build the GRIB2 indexer:

  `make gen_GRIB2_index`

This is a tool that generates an index file of the GRIB records which speeds
access time to individual records substantially.

To install the library, module files and tools, edit the `INSTALLDIR` variable of
the makefile (the default is `/opt/USGS`) and type:

  `make install`

This will also install scripts that can be used to download the 0.5 degree GFS
forecast files and the NCEP 2.5-degree Reannalysis files.

You will need to have write permission in `${INSTALLDIR}` or install as root.


Authors
-------

Hans F. Schwaiger <hschwaiger@usgs.gov>  
Larry G. Mastin <lgmastin@usgs.gov>  
Roger P. Denlinger <rdenlinger@usgs.gov>
