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

useramp=False
usebasefunctions=False

cleanstr=lambda s: re.sub(r'[\r\n]+','\n',s.strip())

epsgwkt={}
crsnamedict={}

crsdef=lambda s: epsgwkt.get(s,s)
crsname=lambda s: crsnamedict.get(s,s)
epsgurl=lambda s: re.sub(r'^EPSG\:','http://www.opengis.net/def/crs/epsg/0/',s)

def main():
    global useramp
    parser=argparse.ArgumentParser('Convert JSON/GeoTIFF model to draft GGXF YAML header')
    parser.add_argument('deformation_json_file',help='Name of deformation model JSON file')
    parser.add_argument('ggxf_yaml_header',nargs='?',help='Name of output YAML GGXF text format file')
    parser.add_argument('-g','--group',nargs='*',help='Groups to include in subset file')
    parser.add_argument('-d','--depth',type=int,help="Maximum grid nesting depth in example file")
    parser.add_argument('-w','--width',type=int,help="Maximum grid nesting width in example file")
    parser.add_argument('-r','--use-ramp',action='store_true',help="Use ramp functions instead of piecewise")
    parser.add_argument('-b','--use-base-time-function',action='store_true',help="Use ramp functions instead of piecewise")

    args=parser.parse_args()
    jsonfile=args.deformation_json_file
    yamlfile=args.ggxf_yaml_header
    if yamlfile is None:
        yamlfile=os.path.splitext(jsonfile)[0]+'.gxb'
    ggxfTimeFunction.useramp=args.use_ramp
    ggxfTimeFunction.usebasefunc=args.use_base_time_function
    smodel=loadJsonGeoTiff(jsonfile)
    gmodel=ggxfModel(smodel,usegroups=args.group,maxdepth=args.depth,maxwidth=args.width)
    dumpGGXFYaml(gmodel,yamlfile)


def pruneGrids( grids, maxdepth, maxwidth ):
    idx={}
    rootgrids=[]
    for g in grids:
        g['children']=[]
        idx[g['gridName']]=g
    for g in grids:
        if p := g.get('parentGridName'):
            idx[p]['children'].append(g)
        else:
            rootgrids.append(g)
    if maxwidth:
        rootgrids=rootgrids[:maxwidth]
    usedgrids=[]
    def trimGrid( g, depth, maxdepth, maxwidth):
        usedgrids.append(g)
        children=g.pop('children')
        if depth == maxdepth:
            children=[]
        elif maxwidth:
            children=children[:maxwidth]
        for c in children:
            trimGrid(c,depth+1,maxdepth,maxwidth)
    for g in rootgrids:
        trimGrid(g,1,maxdepth,maxwidth)
    return usedgrids

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
        gdata['iNodeMaximum'] = size[0]-1
        gdata['jNodeMaximum'] = size[1]-1
        # remark=gmd.get('TIFFTAG_IMAGEDESCRIPTION')
        # if remark:
        #     gdata['remark']=remark
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

def ggxfTimeFunction( tf ):
    global useramp
    tftype=tf['type']
    params=tf['parameters']
    gtf=OrderedDict()
    functions=[]
    if tftype == 'velocity':
        bf=OrderedDict()
        bf['type']='velocity'
        bf['referenceEpoch']=params['reference_epoch']
        functions.append(bf)
    elif tftype == 'step':
        bf=OrderedDict()
        bf['type']='step'
        bf['referenceEpoch']=params['step_epoch']        
        functions.append(bf)
    elif tftype == 'reverse_step':
        bf=OrderedDict()
        bf['type']='reverse_step'
        bf['referenceEpoch']=params['step_epoch']        
        functions.append(bf)  
    elif tftype == 'piecewise':
        model=params['model']
        before=params['before_first']
        after=params['after_last']
        epoch0=model[0]['epoch']
        epochn=model[-1]['epoch']
        if before == 'zero':
            if model[0]['scale_factor'] != 0.0:
                model.insert(0,{'epoch': epoch0,'scale_factor': 0.0})
        elif before != 'constant':
            raise RuntimeError(f'Cannot handle piecewise before_first={before}')
        if after == 'zero':
            if model[-1]['scale_factor'] != 0.0:
                model.append({'epoch': epochn,'scale_factor': 0.0})
        elif after != 'constant':
            raise RuntimeError(f'Cannot handle piecewise after_first={after}')
        lastsf=0.0 
        if ggxfTimeFunction.useramp:           
            for ms,me in zip(model[:-1],model[1:]):
                bf=OrderedDict()
                bf['type']='ramp'
                bf['startEpoch']=ms['epoch']
                bf['startMultiplier']=ms['scale_factor']-lastsf
                bf['endEpoch']=me['epoch']
                bf['endMultiplier']=me['scale_factor']-lastsf
                functions.append(bf)
                lastsf=me['scale_factor']
        else:
            bf=OrderedDict()
            bf['type']='piecewise'
            bf['epochMultipliers']=[
                OrderedDict([('epoch',mp['epoch']),('multiplier',mp['scale_factor'])])
                for mp in model]
            functions.append(bf)
    if ggxfTimeFunction.usebasefunc:
        gtf['baseFunctions']=functions
    else:
        gtf=functions
    return gtf

ggxfTimeFunction.useramp=False    
ggxfTimeFunction.usebasefunc=False    

def ggxfModel(model,usegroups=None,maxwidth=None,maxdepth=None):
    gmodel=OrderedDict()
    gmodel['ggxfVersion']='1.0'
    gmodel['filename']='unknown'
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
    gmodel['operationAccuracy']=0.01
    gmodel['horizontalErrorType']='circular 95% confidence'
    gmodel['verticalErrorType']='95% confidence limit'
    groups=[]
    gmodel['groups']=groups
    for c in model['components']:
        sm=c['spatial_model']
        gname=os.path.basename(sm['filename'])
        gname=os.path.splitext(gname)[0]
        if usegroups and gname not in usegroups:
            print(f"Skipping component {gname}")
            continue
        group=OrderedDict()
        group['groupName']=gname
        if remark := c.get('description'):
            group['remark']=remark
        params=displacementParams[c['displacement_type']][:]
        params.extend(uncertaintyParams[c['uncertainty_type']])
        group['parameters']=params
        group['timeFunction']=ggxfTimeFunction(c['time_function'])
        groups.append(group)
        group['grids']=pruneGrids(sm['grids'],maxdepth,maxwidth)
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
