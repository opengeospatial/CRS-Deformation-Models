#!/usr/bin/python3

import json
import os
import os.path
import sys
import argparse
import yaml
import re
import subprocess
from collections import OrderedDict
from urllib import request

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
    'horizontal': ['displacementEast','displacementNorth'],
    'vertical': ['displacementUp'],
    '3d': ['displacementEast','displacementNorth','displacementUp']
   }

uncertaintyParams={
    'none': []
    }

cleanstr=lambda s: re.sub(r'[\r\n]+','\n',s.strip())

epsgapiurl='https://apps.epsg.org/api/v1/CoordRefSystem/{crdsysid}/export/?format=wkt'

def main():
    parser=argparse.ArgumentParser('Convert JSON/GeoTIFF model to draft GGXF YAML header')
    parser.add_argument('deformation_json_file',help='Name of deformation model JSON file')
    parser.add_argument('ggxf_yaml_header',nargs='?',help='Name of output YAML GGXF text format file')
    parser.add_argument('-g','--group',action='append',help='Groups to include in subset file')
    parser.add_argument('-d','--depth',type=int,help="Maximum grid nesting depth in example file")
    parser.add_argument('-w','--width',type=int,help="Maximum grid nesting width in example file")
    #parser.add_argument('-r','--use-ramp',action='store_true',help="Use ramp functions instead of piecewise")
    #parser.add_argument('-b','--use-base-time-function',action='store_true',help="Use ramp functions instead of piecewise")

    args=parser.parse_args()
    jsonfile=args.deformation_json_file
    yamlfile=args.ggxf_yaml_header
    if yamlfile is None:
        yamlfile=os.path.splitext(jsonfile)[0]+'.gxb'
    # ggxfTimeFunction.useramp=args.use_ramp
    # ggxfTimeFunction.usebasefunc=args.use_base_time_function
    smodel=loadJsonGeoTiff(jsonfile)
    gmodel=ggxfModel(smodel,usegroups=args.group,maxdepth=args.depth,maxwidth=args.width)
    dumpGGXFYaml(gmodel,yamlfile)

def getepsg( crdsysid ):
    if type(crdsysid) == str:
        if crdsysid.lower().startswith('epsg:'):
            crdsysid=crdsysid[5:]
    epsgcache=os.path.expanduser(f'~/.epsg.cache.wkt.{crdsysid}')
    url=epsgapiurl.replace('{crdsysid}',f'{crdsysid}')
    if os.path.isfile(epsgcache):
        wkt=open(epsgcache).read()
    else:
        with request.urlopen(url) as response:
            status=response.getcode()
            if status != 200:
                raise RuntimeError(f"Invalid EPSG code {crdsysid}")
            wkt=response.read().decode('utf8')
            try:
                with open(epsgcache,'w') as ch:
                    ch.write(wkt)
            except:
                try:
                    os.path.remove(epsgcache)
                except:
                    pass
                pass
    name=re.match(r'\s*\w+CRS\[\"([^\"]*)\"',wkt).group(1)
    return OrderedDict([('crsName',name),('crsWkt',wkt),('crsUrl',url)])

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
    rank=1
    try:
        griddata=json.loads(result.stdout,object_pairs_hook=OrderedDict)
        rank = rank+1
        md=griddata['metadata']
        gmd=md['']
        gdata['gridName']=gmd['grid_name']
        nchild=int(gmd.get('number_of_nested_grids',0))
        parent=gmd.get('parent_grid_name')
        if parent:
            gdata['parentGridName']=parent
        gdata['hierarchyRank']=rank
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
        bf['timeFunctionType']='velocityTimeFunction'
        bf['functionReferenceDate']=params['reference_epoch']
        functions.append(bf)
    elif tftype == 'step':
        bf=OrderedDict()
        bf['timeFunctionType']='stepTimeFunction'
        bf['stepDate']=params['step_epoch']        
        functions.append(bf)
    elif tftype == 'reverse_step':
        bf=OrderedDict()
        bf['timeFunctionType']='reverseStepTimeFunction'
        bf['stepDate']=params['step_epoch']        
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
                bf['timeFunctionType']='piecewiseTimeFunction'
                bf['startDate']=ms['epoch']
                bf['startScaleFactor']=ms['scale_factor']-lastsf
                bf['endDate']=me['epoch']
                bf['endScaleFactor']=me['scale_factor']-lastsf
                functions.append(bf)
                lastsf=me['scale_factor']
        else:
            raise RuntimeError('piecewise function not supported by current UML')
            bf=OrderedDict()
            bf['timeFunctionType']='piecewiseTimeFunction'
            bf['epochMultipliers']=[
                OrderedDict([('epoch',mp['epoch']),('multiplier',mp['scale_factor'])])
                for mp in model]
            functions.append(bf)
    if ggxfTimeFunction.usebasefunc:
        raise RuntimeError('baseFunctions keyword not supported by current UML')
        gtf['baseFunctions']=functions
    else:
        gtf=functions
    return gtf

ggxfTimeFunction.useramp=True    
ggxfTimeFunction.usebasefunc=False    

def ggxfModel(model,usegroups=None,maxwidth=None,maxdepth=None):
    gmodel=OrderedDict()
    gmodel['ggxfVersion']='1.0'
    gmodel['fileName']='unknown'
    gmodel['version']=model['version']
    gmodel['content']='deformation model'
    gmodel['remark']=cleanstr(model['description'])
    authority=OrderedDict()
    authority['organisationName']=model['authority']['name']
    address=model['authority']['address']
    address=address.replace("\r\n","\n")
    city='Unknown'
    postcode='Unknown'
    # Crude implementation!
    address,city=address.rsplit("\n",maxsplit=1)
    if match := re.match(r"^(.*?)\s*(\d+)\s*$",city):
        city=match.group(1)
        postcode=match.group(2)
    authority['addressDeliveryPoint']=address
    authority['addressCity']=city
    authority['postalCode']=postcode
    authority['electronicMailAddress']=model['authority']['email']
    about=[link for link in model['links'] if link['rel'] == 'about']
    if about:
        authority['onlineResourceLinkage']=about[0]['href']
    gmodel['authority']=authority
    
    gmodel['publicationDate']=model['publication_date']

    extent=model['extent']['parameters']['bbox']
    gmodel['contentApplicabilityExtent']=OrderedDict((
        ('southBoundLatitude',extent[1]),
        ('westBoundLongitude',extent[0]),
        ('northBoundLatitude',extent[3]),
        ('eastBoundLongitude',extent[2]),
        ('extentDescription','New Zealand EEZ')

    ))
    gmodel['startDate']=model['time_extent']['first']
    gmodel['endDate']=model['time_extent']['last']
    gmodel['sourceCrs']=getepsg(model['source_crs'])
    gmodel['targetCrs']=getepsg(model['target_crs'])
    gmodel['interpolationCrs']=getepsg(model['definition_crs']) # gridcrs is WKT in current file - ignoring
    gmodel['operationAccuracy']=0.01
    gmodel['deformationUncertaintyType']='2CE/2SE'
    gmodel['deformationApplicationMethod']='addition'
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
        params=[{            
            'parameterID': p,
            'unitID': 'metre',
            'unitType': 'length',
            'unitSiRatio': 1.0
            } for p in displacementParams[c['displacement_type']]
            ]
        group['parameters']=params
        group['deformationTimeFunction']=ggxfTimeFunction(c['time_function'])
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
    gmodel['fileName']=os.path.basename(yamlfile)
    yamldef=yaml.dump(gmodel,width=2048)
    yamldef=blockLongLines(yamldef)
    yamldef=blockAffine(yamldef)
    check=yaml.load(yamldef,Loader=yaml.Loader)
    open(yamlfile,'w').write(yamldef)

if __name__== '__main__':
    main()
