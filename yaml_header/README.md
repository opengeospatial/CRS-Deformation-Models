deformation_model_to_yaml.py
============================

Work in progress on non-generallised JSON/GeoTIFF to GGXF text format 
converter converting a JSON/GeoTIFF deformation model format file.  

The current version has been used to build a subset of the LINZ deformation model,
in nz_linz_nzgd2000-20180701.yaml.  The grids used by this are in nzdm/

ggxf_yaml_to_netcdf4.py
=======================

Very early proof of concept of conversion of GGXF to a NetCDF4 file.  Usage:

```
python3 ggxf_yaml_to_netcdf4.py -g nzdm nz_linz_nzgd2000-20180701.yaml nz_linz_nzgd2000-20180701.nc
```

Notes on this implementation:

* This only has a minimal subset of GGXF metadata attributes to demonstrate possible encoding
* This is built within the constraints of the NetCDF4 python library (1.5.7, libnetcdf 18.0.0, libhdf5 200.0.0 in python site-packages)
* For a full implementation it would probably make sense to generate the code for saving the NetCDF file from a GGXF schema, or maybe use a schema directly in the code to write the NetCDF file.


Questions for implementation:

* Attributes vs. variables. Complex types and array attributes (eg such as GGXF group parameters and grid affineCoefs) are stored as NetCDF Variables, whereas scalar values such as remark are stored as NetCDF attributes.  Is the dichotomy the appropriate way to store these
* Complex types vs groups.  Some complex types could be stored as subgroups with a set of attributes.  For example some structured metadata such as citation could be stored as a subgroup with attributes (or possibly a nested set of subgroups with attributes).  This might be simple and less limiting than compound types.  In particular compound types (at least with the Python implementation) only support fixed length strings.
* At the September OGC 2021 meeting it was suggested we look at [netcdf-ld](https://github.com/opengeospatial/netcdf-ld).  I have looked at this, though not in detail.  To be honest I am struggling to figure out how to use it or fit it in to this work.  However it may provide a model for encoding metadata.  