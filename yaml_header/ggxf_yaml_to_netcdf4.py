#!/usr/bin/python3
import yaml
import argparse
import netCDF4
import numpy as np
from osgeo import gdal
import os
import os.path
import logging
from collections import namedtuple
import re

logging.basicConfig()
logging.getLogger().setLevel(logging.DEBUG)


class Error(RuntimeError):
    pass


class SchemaError(Error):
    pass


class GGXF:

    RecordType=namedtuple('RecordType','dtype netcdfType')

    def __init__(self, yaml_file: str = None, options: dict = None):
        self._data = None
        self._filename = None
        self._logger = logging.getLogger("GGXF")
        self._options = options or {}
        if yaml_file:
            self.loadYaml(yaml_file)

    def loadYaml(self, yaml_file: str) -> None:
        if not os.path.isfile(yaml_file):
            raise Error(f"{yaml_file} does not exist or is not a file")
        try:
            self._logger.debug(f"Loading YAML {yaml_file}")
            self._data = yaml.load(open(yaml_file).read(), Loader=yaml.SafeLoader)
            self._filename = yaml_file
        except:
            raise Error(f"Cannot load YAML file {yaml_file}")
        self.validate()
        self.loadGrids()

    def loadGrids(self) -> None:
        startdir = os.getcwd()
        griddir = self._options.get("grid_directory", startdir)
        self._logger.debug(f"Changing directory to {griddir} to load grids")
        os.chdir(griddir)
        try:
            for group in self._data["groups"]:
                group['_root']=self._data
                for grid in group["grids"]:
                    grid['_group']=group
                    if "gridDataSource" in grid:
                        source = grid.pop("gridDataSource")
                        self._logger.debug(f"Loading {source}")
                        dataset = gdal.Open(source)
                        gridData = dataset.ReadAsArray()
                        gridData = gridData.transpose()
                        grid["gridData"] = gridData
                        self._logger.debug(f"Loaded grid with dimension {gridData.shape}")
                        # Get grid metadata to validate against YAML...
        finally:
            os.chdir(startdir)

    def validate(self) -> None:
        """
        Not implemented: Currently have a placeholder implementation pending a genuine schema implementation
        """
        data = self._data
        if "groups" not in data or not isinstance(data["groups"], list):
            raise SchemaError("groups not defined or not a list")
        for igroup, group in enumerate(data["groups"]):
            if "groupName" not in group:
                raise SchemaError(f"groupName not defined in groups[{igroup}]")
            if "grids" not in group or not isinstance(group["grids"], list):
                raise SchemaError(
                    f"grids not defined or not a list in group {group['groupName']}"
                )

    def save(self, netcdf4_file: str):
        self._logger.debug(f"Saving NetCDF4 grid as {netcdf4_file}")
        if os.path.isfile(netcdf4_file):
            os.remove(netcdf4_file)
        root=netCDF4.Dataset(netcdf4_file,"w",format="NETCDF4")
        self._saveNetCdf4( root )

    def _saveNetCdf4( self, root ):
        data=self._data
        nctypes={}

        # Common dimensions
        root.createDimension("affine",6)
        # Tried using unlimited dimension for parameter list.  Cause ncdump 4.7.3 to crash.
        # root.createDimension("parameter",None)

        # Common compound types
        paramtype=np.dtype([("parameterName",'S32'),("unit",'S16'),("unitSiRatio",np.float64)])
        ncparamtype=root.createCompoundType(paramtype,'ggxfParameterType')
        nctypes['ggxfParameterType']=GGXF.RecordType(paramtype,ncparamtype)

        # Base object metdata
        root.setncatts({
            'version':data['version'],
            'content':data['content'],
            'operationAccuracy':float(data['operationAccuracy'])
            })

        # Store each of the groups
        for group in data['groups']:
            self._saveGroupNetCdf4(root,group,nctypes)

    def _saveGroupNetCdf4( self, root:netCDF4.Dataset, group: dict, nctypes: dict ):
        name=group.get('groupName')
        cdfgroup=root.createGroup(name)

        # Store group attributes
        cdfgroup.setncatts({
            'remark': group.get('remark','')
        })
        # Store the parameters saved as array of compound type
        parameters=group.get('parameters',[])
        cdfgroup.createDimension('nParam',len(parameters))
        paramdata=[(p['parameterName'],p['unit'],p['unitSiRatio']) for p in parameters]
        paramtype=nctypes['ggxfParameterType']
        paramdata=np.array(paramdata,dtype=paramtype.dtype)
        paramvar=cdfgroup.createVariable('parameters',paramtype.netcdfType,'nParam')
        paramvar[:]=paramdata

        # Store each of the grids
        for grid in group['grids']:
            self._saveGridNetCdf4(cdfgroup,grid,nctypes)

    def _saveGridNetCdf4( self, cdfgroup: netCDF4.Group, grid:dict, nctypes: dict  ):
        name=grid.get('gridName')
        cdfgrid=cdfgroup.createGroup(name)

        # Store group attributes
        cdfgrid.setncatts({
            'remark': grid.get('remark','')
        })

        # Store the grid metadata
        affinevar=cdfgrid.createVariable('affineCoeffs',np.float64,'affine')
        self._logger.debug(f"Filling affineCoeffs {affinevar.shape} {affinevar.datatype} from {grid['affineCoeffs']}")
        affinevar[:]=np.array(grid['affineCoeffs'])
        cdfgrid.createDimension('nCol',grid['iNodeMaximum']+1)
        cdfgrid.createDimension('nRow',grid['jNodeMaximum']+1)

        # Store the grid data
        datavar=cdfgrid.createVariable('data',np.float64,['nCol','nRow','nParam'])
        datavar[:,:,:]=grid['gridData']

def main():
    parser = argparse.ArgumentParser(
        "Convert a YAML GGXF specification to a NetCDF4 file"
    )
    parser.add_argument("ggxf_yaml_file", help="Name of YAML GGXF file to load")
    parser.add_argument("netcdf4_file", help="NetCDF4 file to create")
    parser.add_argument(
        "-g", "--grid-directory", help="Directory to search for grid files"
    )
    args = parser.parse_args()

    yaml_file = args.ggxf_yaml_file
    netcdf4_file = args.netcdf4_file
    print(args)

    options = {}
    if args.grid_directory:
        options["grid_directory"] = args.grid_directory

    ggxf = GGXF(yaml_file, options=options)
    ggxf.save(netcdf4_file)


if __name__ == "__main__":
    main()
