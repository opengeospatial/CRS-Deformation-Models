#!/usr/bin/python3

import json
import os.path
import sys
import argparse
import yaml
import re
import subprocess
from collections import OrderedDict

def represent_ordereddict(dumper, data):
    value = []

    for item_key, item_value in data.items():
        node_key = dumper.represent_data(item_key)
        node_value = dumper.represent_data(item_value)

        value.append((node_key, node_value))

    return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)

yaml.add_representer(OrderedDict, represent_ordereddict)

displacementParams={
    'none': [],
    'horizontal': ['displacementEast metre length 1.0','displacementNorth metre length 1.0'],
    'vertical': ['displacementUp metre length 1.0'],
    '3d': ['displacementEast metre length 1.0','displacementNorth metre length 1.0','displacementUp metre length 1.0']
   }

uncertaintyParams={
    'none': []
    }

cleanstr=lambda s: re.sub(r'[\r\n]+','\n',s.strip())

epsgwkt={}
crsnamedict={}

crsdef=lambda s: epsgwkt.get(s,s)
crsname=lambda s: crsnamedict.get(s,s)
epsgurl=lambda s: re.sub(r'^EPSG\:','http://www.opengis.net/def/crs/epsg/0/',s)

def main():
    parser=argparse.ArgumentParser('Convert JSON/GeoTIFF model to draft GGXF YAML header')
    parser.add_argument('deformation_json_file',help='Name of deformation model JSON file')
    parser.add_argument('ggxf_yaml_header',nargs='?',help='Name of output YAML GGXF text format file')

    args=parser.parse_args()
    jsonfile=args.deformation_json_file
    yamlfile=args.ggxf_yaml_header
    if yamlfile is None:
        yamlfile=os.path.splitext(jsonfile)[0]+'.gxb'
    smodel=loadJsonGeoTiff(jsonfile)
    gmodel=ggxfModel(smodel)
    dumpGGXFYaml(gmodel,yamlfile)


def loadGTiffGridData(source,sourceref=None,tiffdir=None):
    print(f'Loading {source}')
    result=subprocess.run(['gdalinfo','-json',source],capture_output=True,cwd=tiffdir or None)
    if result.returncode != 0:
        raise RuntimeError(f'Failed to load GeoTiff {result.returncode}: {result.stderr}')
    gdata=OrderedDict()
    subgrids=[]
    crsdef=''
    try:
        griddata=json.loads(result.stdout,object_pairs_hook=OrderedDict)
        md=griddata['metadata']
        gmd=md['']
        gdata['gridName']=gmd['grid_name']
        nchild=int(gmd.get('number_of_nested_grids',0))
        if nchild:
            gdata['numberOfChildGrids']=nchild
        parent=gmd.get('parent_grid_name')
        if parent:
            gdata['parentGridName']=parent
        affine=list(griddata['geoTransform'])
        affine[0] += affine[1]/2.0
        affine[3] += affine[5]/2.0
        coeffs=affine[3:]
        coeffs.extend(affine[:3])
        gdata['affineCoeffs']=coeffs
        size=griddata['size']
        gdata['iNodeMaximum'] = size[0]
        gdata['jNodeMaximum'] = size[1]
        remark=gmd.get('TIFFTAG_IMAGEDESCRIPTION')
        if remark:
            gdata['remark']=remark
        gdata['gridDataSource']=sourceref or source
        for k,v in md.get('SUBDATASETS',{}).items():
            if m := re.match(r'^SUBDATASET_(\d+)_NAME',k):
                if m.group(1) != '1':
                    subgrids.append(v)
        crsdef=griddata['coordinateSystem']['wkt']
    except Exception as ex:
        raise RuntimeError(f'Failed to load {source}: {ex}')
    return gdata, subgrids, crsdef

def loadJsonGeoTiff(jsonfile,object_pairs_hook=OrderedDict):
    model=json.loads(open(jsonfile).read(),object_pairs_hook=OrderedDict)
    tiffdir=os.path.dirname(jsonfile)
    grids=[]
    for c in model['components']:
        sm=c['spatial_model']
        source=sm['filename']
        gdalref='GTIFF_DIR:1:'+source
        gdata,subgrids,crsdef =loadGTiffGridData(source,tiffdir=tiffdir,sourceref=gdalref)
        grids=[gdata]
        for subgridref in subgrids:
            gdata,subgrids,crsdef=loadGTiffGridData(subgridref,tiffdir=tiffdir)
            grids.append(gdata)
        sm['grids']=grids
    model['gridcrs']=crsdef
    return model

def ggxfTimeFunction( tf, extent ):
    tftype=tf['type']
    params=tf['parameters']
    gtf=OrderedDict()
    gtf['minEpoch']=extent['first']
    gtf['maxEpoch']=extent['last']
    functions=[]
    gtf['baseFunctions']=functions
    if tftype == 'velocity':
        bf=OrderedDict()
        bf['type']='polynomial'
        bf['referenceEpoch']=params['reference_epoch']
        bf['velocityMultiplier']=1.0
        functions.append(bf)
    elif tftype == 'step':
        epoch=params['step_epoch']
        gtf['minEpoch']=epoch
        bf=OrderedDict()
        bf['type']='ramp'
        bf['startEpoch']=epoch
        bf['startMultiplier']=0.0
        bf['endEpoch']=epoch
        bf['endMultiplier']=1.0
        functions.append(bf)
    elif tftype == 'reverse_step':
        epoch=params['step_epoch']
        gtf['maxEpoch']=epoch
        bf=OrderedDict()
        bf['type']='ramp'
        bf['startEpoch']=epoch
        bf['startMultiplier']=1.0
        bf['endEpoch']=epoch
        bf['endMultiplier']=0.0    
    elif tftype == 'piecewise':
        model=params['model']
        before=params['before_first']
        after=params['after_last']
        if before == 'zero':
            gtf['minEpoch']=model[0]['epoch']
        elif before != 'constant':
            raise RuntimeError(f'Cannot handle piecewise before_first={before}')
        if after == 'zero':
            gtf['maxEpoch']=model[-1]['epoch']
        elif after != 'constant':
            raise RuntimeError(f'Cannot handle piecewise after_first={after}')
        lastsf=0.0            
        for ms,me in zip(model[:-1],model[1:]):
            bf=OrderedDict()
            bf['type']='ramp'
            bf['startEpoch']=ms['epoch']
            bf['startMultiplier']=ms['scale_factor']-lastsf
            bf['endEpoch']=me['epoch']
            bf['endMultiplier']=me['scale_factor']-lastsf
            functions.append(bf)
            lastsf=me['scale_factor']
    return gtf

def ggxfModel(model):
    gmodel=OrderedDict()
    gmodel['ggxfVersion']='1.0'
    gmodel['content']='deformation model'
    gmodel['version']=model['version']
    gmodel['remark']=cleanstr(model['description'])
    gmodel['organisationName']=model['authority']['name']
    address=model['authority']['address']
    city='Unknown'
    postcode='Unknown'
    # Crude implementation!
    m=re.match(r'(.*)\r?\n([\w\s]*?)\s*(\d+)\s*$',address)
    if m:
        address=m.group(1)
        city=m.group(2)
        postcode=m.group(3)
    gmodel['addressDeliveryPoint']=address
    gmodel['addressCity']=city
    gmodel['postalCode']=postcode
    gmodel['electronicMailAddress']=model['authority']['email']
    about=[link for link in model['links'] if link['rel'] == 'about'][0]
    gmodel['onlineResourceLinkage']=about['href']
    gmodel['publicationDate']=model['publication_date']

    extent=model['extent']['parameters']['bbox']
    gmodel['contentApplicabilityExtent']=OrderedDict((
        ('southBoundLatitude',extent[1]),
        ('westBoundLongitude',extent[0]),
        ('northBoundLatitude',extent[3]),
        ('eastBoundLongitude',extent[2])

    ))
    gmodel['contentApplicabilityExtentDescription']='New Zealand EEZ'
    gmodel['startDate']=model['time_extent']['first']
    gmodel['endDate']=model['time_extent']['last']
    gmodel['interpolationCrsWkt']=crsdef(model.get('gridcrs',model['definition_crs']))
    gmodel['sourceCrsName']=crsname(model['source_crs'])
    gmodel['sourceCrsHref']=epsgurl(model['source_crs'])
    gmodel['targetCrsName']=crsname(model['target_crs'])
    gmodel['targetCrsHref']=epsgurl(model['target_crs'])
    gmodel['operationalAccuracy']='0.01'
    gmodel['horizontalErrorType']='circular 95% confidence'
    gmodel['verticalErrorType']='95% confidence limit'
    groups=[]
    gmodel['groups']=groups
    for c in model['components']:
        sm=c['spatial_model']
        group=OrderedDict()
        gname=os.path.basename(sm['filename'])
        gname=os.path.splitext(gname)[0]
        group['groupName']=gname
        if remark := c.get('description'):
            group['remark']=remark
        params=displacementParams[c['displacement_type']][:]
        params.extend(uncertaintyParams[c['uncertainty_type']])
        group['parameters']=params
        group['timeFunction']=ggxfTimeFunction(c['time_function'],model['time_extent'])
        groups.append(group)
        group['grids']=sm['grids']
    return gmodel

def blockLongLines(yamldef):
    def createBlock(match):
        indent=match.group(1)
        prefix=match.group(2)
        value=yaml.load(match.group(3),Loader=yaml.Loader)
        value=value.replace('\r','')
        lines=value.split('\n')
        lines.insert(0,indent+prefix+'|')
        return f'\n{indent}  '.join(lines)
    def createBlock2(match):
        indent=match.group(1)
        prefix=match.group(2)
        value=yaml.load(match.group(3),Loader=yaml.Loader)
        if '\n' not in value:
            return match.group(0)
        value=value.replace('\r','')
        lines=[l.strip() for l in value.split('\n')]
        lines.insert(0,indent+prefix+'|')
        return f'\n{indent}  '.join(lines)        
    blockre=r'^(\s*)(\w+\:\s+)(\".*\\n.*\")\s*$'
    yamldef=re.sub(blockre,createBlock,yamldef,flags=re.M)
    blockre=r'^(\s*)(\w+\:\s+)(\'(?:[^\']+|\\\')*\')\s*$'
    yamldef=re.sub(blockre,createBlock2,yamldef,flags=re.M | re.S)
    return yamldef

def blockAffine(yamldef):
    def affsub(m):
        prefix=m.group(1)+m.group(2)
        coeffs=[v for v in m.group(3).split() if v != '-']
        return prefix+' ['+','.join(coeffs)+']'
    affre=r'^(\s*)(affineCoeffs\:)\s*$((?:\n\1\s*\-\s+\S+\s*$){6})'
    yamldef=re.sub(affre,affsub,yamldef,flags=re.M | re.S)
    return yamldef

def dumpGGXFYaml(gmodel,yamlfile):
    gmodel['filename']=os.path.basename(yamlfile)
    yamldef=yaml.dump(gmodel,width=2048)
    yamldef=blockLongLines(yamldef)
    yamldef=blockAffine(yamldef)
    check=yaml.load(yamldef,Loader=yaml.Loader)
    open(yamlfile,'w').write(yamldef)

if __name__== '__main__':
    main()
