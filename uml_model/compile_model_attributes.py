#!/usr/bin/python3
from collections import OrderedDict, namedtuple
import re
import argparse
import csv

class Class():

    class Attribute():

        def __init__(self,name=None,aclass=None,dtype=None,isarray=False,optional=False,copy=None):
            copyatt=lambda att: None if copy is None else getattr(copy,att)
            self.name=name or copyatt('name')
            self.aclass=aclass or copyatt('aclass')
            self.dtype=dtype if dtype is not None else copyatt('dtype')
            self.isarray=isarray or copyatt('isarray')
            self.optional=optional or copyatt('optional')

        def copy( self, **params ):
            return Class.Attribute(copy=self,**params)

    baseTypes=['characterString','real','integer','date','dateTime','']
 
    def __init__(self,name):
        self.name=name
        self.rootclass=None
        self.baseclass=[]
        self.subclass=[]
        self.attributes=[]
        self.root=False
        self.used=False

    def allAttributes( self ):
        attlist=[]
        for baseclass in self.baseclass:
            attlist.extend(baseclass.allAttributes())
        attlist.extend(self.attributes)
        return attlist

    def isDeformation(self):
        if 'Deformation' in self.name:
            return True
        for baseclass in self.baseclass:
            if baseclass.isDeformation():
                return True
        return False

class ClassModel:

    def __init__(self, plantumlfile):
        self.filename=plantumlfile
        self.classIndex=OrderedDict()
        self._parse()
        self._checkAttributeTypes()

    def _parse( self ):
        with open(self.filename) as pfh:
            source=pfh.read()
        classre=re.compile(r'class\s+(\w+)\s*\{([^\}]*)\}',re.S)
        attre=re.compile(r'\s*(\w+)(?:\s\((optional)\))?\s*\:\s+(\?|\w+\??)(\[.*\])?\s*$')
        classindex=self.classIndex
        isroot=True
        for match in classre.finditer(source):
            classname=match.group(1)
            classobj=Class(classname)
            if isroot:
                self.rootclass=classobj
                classobj.root=True
                isroot=False
            classindex[classname]=classobj
            attlist=match.group(2)
            for attdef in attlist.split('\n'):
                amatch=attre.match(attdef)
                if not amatch:
                    if not re.match(r'^\s*$',attdef):
                        print(f"Fail att {attdef}")
                    continue
                name=amatch.group(1)
                optional=amatch.group(2) == 'optional'
                dtype=amatch.group(3).replace('?','')
                isarray=bool(amatch.group(4))
                classobj.attributes.append(Class.Attribute(name,aclass=classobj,dtype=dtype,isarray=isarray,optional=optional))

        inheritre=re.compile(r'^\s*(\w+)\s+\<\|\-\-\s+(\w+)\s*$',re.M)       
        for imatch in inheritre.finditer(source):
            basename=imatch.group(1)
            classname=imatch.group(2)
            classindex[classname].baseclass.append(classindex[basename])
            classindex[basename].subclass.append(classindex[classname])

    def _checkAttributeTypes(self):
        for aclass in self.classIndex.values():
            for att in aclass.attributes:
                if att.dtype in self.classIndex:
                    self.classIndex[att.dtype].used=True
                elif att.dtype not in Class.baseTypes:
                    print(f"Invalid attribute type {att.dtype} for {aclass.name}.{att.name}")


    def compileAttributes(self):
        attributes={}
        for aclass in self.classIndex.values():
            for att in aclass.allAttributes():
                if att.name not in attributes:
                    attributes[att.name]=[]
                attributes[att.name].append(att.copy(aclass=aclass))
        return attributes


attributeDefinitionFields='name datatype domain context definition remark'.split()
AttributeDefinition=namedtuple('AttributeDefinition',attributeDefinitionFields)

def mergeAttributes( attlist, attdef, domain ):
    merged={}
    merged['name'] = attdef.get('name',attlist[0].name)
    merged['definition'] = attdef.get('definition','')
    merged['remark'] = attdef.get('remark','')
    merged['domain'] = domain
    datatypes=set()
    context=set()
    for a in attlist:
        # Assume is ABC if has subclass
        if not a.aclass.subclass:
            typename=a.dtype+('[]' if a.isarray else '')
            datatypes.add(typename)
            context.add(a.aclass.name)

    merged['datatype']=', '.join(sorted(datatypes))
    merged['context']=', '.join(sorted(context))
    return AttributeDefinition(**merged)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('plantuml_source_file')
    parser.add_argument('root_output_name')
    parser.add_argument('-d','--ggxf_keyword_description_csv',action='append')
    args=parser.parse_args()
    model=ClassModel(args.plantuml_source_file)
    attributes=model.compileAttributes()
    attspec={}
    attspecFields=['name','definition','remark']
    for specfile in args.ggxf_keyword_description_csv:
        with open(specfile) as gth:
            csvh=csv.DictReader(gth)
            for row in csvh:
                attname=row.get('name')
                if attname:
                    attspec[attname]={f:row.get(f,'') for f in attspecFields}

    deformation_attributes=[]
    unused_attributes=[attspec[a] for a in attspec if a not in attributes]

    for attname in sorted(attributes):
        attlist=attributes[attname]
        clist=[att.aclass.name for att in attlist]
        domain='deformation'
        for att in attlist:
            if not att.aclass.isDeformation():
                domain='ggxf'
        merged=mergeAttributes(attlist,attspec.get(attname,{}),domain)
        deformation_attributes.append(merged)
    
    with open(args.root_output_name+'_deformation.csv','w') as outh:
        csvh=csv.DictWriter(outh,attributeDefinitionFields)
        csvh.writeheader()
        for a in deformation_attributes:
            csvh.writerow(a._asdict())

    with open(args.root_output_name+'_unused.csv','w') as outh:
        csvh=csv.DictWriter(outh,attspecFields)
        csvh.writeheader()
        for a in unused_attributes:
            csvh.writerow(a)        