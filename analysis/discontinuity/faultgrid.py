from importlib.abc import Loader
import os.path
import csv
import numpy as np
import numpy.linalg as npla
import math
# from okada85 import okada85
import yaml
import argparse
import subprocess
from shapely import geometry, wkt
# from shapely.geometry import MultiLineString
import shapely
import os
import sys

calc_okada=os.path.join(os.path.dirname(os.path.realpath(__file__)),"calc_okada")
if sys.platform.startswith('win'):
    calc_okada=calc_okada+'.exe'
if not os.path.exists(calc_okada):
    raise RuntimeError(f"Cannot find calc_okada at {calc_okada}")


def main():
    parser=argparse.ArgumentParser(description="Construct a fault model for assessing discontinuity options")
    parser.add_argument("fault_def",help="calc_okada fault definition file")
    parser.add_argument("model_def",help="JSON formatted model file")
    args=parser.parse_args()

    if not os.path.exists(args.fault_def):
        raise RuntimeError(f"Cannot find {args.fault_def}")


    with open(args.model_def) as mdh:
        model=yaml.load(mdh.read(),Loader=yaml.Loader)
    rootname=os.path.splitext(args.model_def)[0]

    wktfile=f"{rootname}_fault.wkt"
    subprocess.run([calc_okada,"-w",wktfile,args.fault_def])
    faultgeom=readFaultWkt(wktfile)

    grid,subgrid,interp,spacing=construct_grid(model["grid"])
    gdisp=calcOkada(grid,args.fault_def)
    np.savetxt(f"{rootname}_grid.csv",
        np.hstack([gridPoints(grid),gridPoints(gdisp)]),
        "%.4f",",",header="x,y,dx,dy,dz",comments="")
    with open(f"{rootname}_grid.wkt","w") as wkth:
        wkth.write("wkt\n")
        for i in range(grid.shape[0]):
            wkth.write(geometry.LineString(grid[i,:,:]).wkt)
            wkth.write("\n")
        for i in range(grid.shape[1]):
            wkth.write(geometry.LineString(grid[:,i,:]).wkt)
            wkth.write("\n")

    finterp=lambda xy: interp(xy,gdisp)
    sdata=calcOkada(subgrid,args.fault_def)
    idata=gridMap(subgrid,finterp,3)
    sdiff=(sdata-idata)[:,:,:2]
    difflen=npla.norm(sdiff,axis=2)
    difflen=difflen.reshape(difflen.shape[0],difflen.shape[1],1)

    np.savetxt(f"{rootname}_interp.csv",
        np.hstack([gridPoints(subgrid),gridPoints(sdata),gridPoints(idata),gridPoints(difflen)]),
    "%.4f",",",header="x,y,dx,dy,dz,dxi,dyi,dzi,eh",comments="")

    nodes,cells=intersectingCells(grid,faultgeom)
    cellbuffer=geometry.MultiPolygon(cells).buffer(0)
    with open(f"{rootname}_flagged_cells.wkt","w") as wkth:
        wkth.write("id\tshape\n")
        wkth.write(f"flagged cells\t{cellbuffer.wkt}\n")

    buffer=npla.norm(spacing)/2.0
    inrange=gridPointsInRange(grid,faultgeom,buffer)
    err=1.0*inrange+0.1
    err=err.reshape((err.shape[0],err.shape[1],1))
    edata=gridMap(subgrid,lambda xy: interp(xy,err),1)
    np.savetxt(f"{rootname}_err_interp.csv",
        np.hstack([gridPoints(subgrid),gridPoints(edata)]),
    "%.4f",",",header="x,y,err",comments="")
    np.savetxt(f"{rootname}_err.csv",
        np.hstack([gridPoints(grid),gridPoints(err)]),
    "%.4f",",",header="x,y,err",comments="")
    flagged=np.zeros((grid.shape[0],grid.shape[1],1))
    flagged[nodes[:,0],nodes[:,1],0]=1
    flagged[nodes[:,0]+1,nodes[:,1],0]=1
    flagged[nodes[:,0]+1,nodes[:,1]+1,0]=1
    flagged[nodes[:,0],nodes[:,1]+1,0]=1

    err=1.0*flagged+0.1
    edata=gridMap(subgrid,lambda xy: interp(xy,err),1)
    np.savetxt(f"{rootname}_err2_interp.csv",
        np.hstack([gridPoints(subgrid),gridPoints(edata)]),
    "%.4f",",",header="x,y,err",comments="")
    np.savetxt(f"{rootname}_err2.csv",
        np.hstack([gridPoints(grid),gridPoints(err)]),
    "%.4f",",",header="x,y,err",comments="")

    # with open(f"{rootname}_data.csv","w") as fh:
    #     csvw=csv.writer(fh)
    #     csvw.writerow(("id","shape"))
    #     csvw.writerow(("fault",fgeom.wkt))


# def construct_fault(fault_def):
#     top=np.array(fault_def["top"],dtype=np.float) # {"top": [[e1,n1],[e2,n2],[...]]}
#     if len(top.shape) != 2 or top.shape[0] < 2 or top.shape[1] != 2:
#         raise RuntimeError("Invalid definition of top of fault")
#     depth=float(fault_def.get("depth",0.0))
#     width=float(fault_def["width"])
#     dip=float(fault_def.get("dip",90.0))
#     dislocation=float(fault_def["dislocation"])
#     rake=float(fault_def.get("rake",0.0))
#     rrake=math.radians(rake)
#     u1=dislocation*math.cos(rrake)
#     u2=dislocation*math.sin(rrake)
#     nsegments=top.shape[0]-1
#     funcs=[]
#     for xy0,xy1 in zip(top[:-1],top[1:]):
#         diff=xy1-xy0
#         strike=math.degrees(math.atan2(diff[0],diff[1]))
#         length=math.sqrt(diff[0]*diff[0]+diff[1]*diff[1])
#         func=lambda xy: np.array(okada85(xy[0]-xy0[0],xy[1]-xy0[1],depth,strike,dip,length,width,u1,u2,0.0)[:3])
#         funcs.append(func)
#     calc=lambda xy: sum((func(xy) for func in funcs))
#     geom=LineString(top)
#     return (calc,geom)

def construct_grid(grid_def):
    size=np.array(grid_def["size"])
    spacing=np.array(grid_def["spacing"])
    se=np.array(grid_def["se"])
    nw=se+(size-1)*spacing
    xvalues=np.linspace(se[0],nw[0],size[0])
    yvalues=np.linspace(se[1],nw[1],size[1])
    xv,yv=np.meshgrid(xvalues,yvalues)
    grid=np.array([xv,yv]).T

    subgridsize=(size-1)*int(grid_def.get("subdivision",5))+1
    xvalues=np.linspace(se[0],nw[0],subgridsize[0])
    yvalues=np.linspace(se[1],nw[1],subgridsize[1])
    xv,yv=np.meshgrid(xvalues,yvalues)
    subgrid=np.array([xv,yv]).T

    maxij=np.array([size[0]-2,size[1]-2])
    def interp(xy,data):
        wij=(xy-se)/spacing
        cij=np.maximum([0,0],np.minimum(maxij,np.floor(wij).astype(np.int)))
        wij -= cij
        return sum((
            data[cij[0],cij[1]]*(1-wij[0])*(1-wij[1]),
            data[cij[0],cij[1]+1]*(1-wij[0])*(wij[1]),
            data[cij[0]+1,cij[1]+1]*(wij[0])*(wij[1]),
            data[cij[0]+1,cij[1]]*(wij[0])*(1-wij[1]),
        ))
    return grid,subgrid,interp,spacing

def gridPoints(grid):
    shape=grid.shape
    return grid.reshape((shape[0]*shape[1],shape[2]))

def gridMap( grid, func, nfunc=3 ):
    data=np.zeros((grid.shape[0],grid.shape[1],nfunc))
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            data[i,j]=func(grid[i,j])
    return data

def calcOkada( grid, fault_def):
    tmpin="calc_okada_in.tmp"
    tmpout="calc_okada_out.tmp"
    np.savetxt(tmpin,gridPoints(grid),"%.4f")
    subprocess.run([calc_okada,fault_def,tmpin,tmpout])
    defm=np.loadtxt(tmpout,skiprows=1)
    os.remove(tmpin)
    os.remove(tmpout)
    defm=defm[:,2:]
    return defm.reshape(grid.shape[0],grid.shape[1],3)


def calcFault(fault,points):
    calc=fault["calc"]
    return np.array([calc(p) for p in points])

def readFaultWkt( wktfile ):
    geoms=[]
    with open(wktfile) as wkth:
        for l in wkth:
            fwkt=l.strip().split('\t',)[-1]
            if 'MULTILINESTRING' in fwkt:
                geom=wkt.loads(fwkt)
                tracegeom=geom[1]
                geoms.append(tracegeom)
    return geometry.MultiLineString(geoms)

def intersectingCells( grid, geom ):
    indices=[]
    cells=[]
    for i,ci in enumerate(zip(grid[:-1],grid[1:])):
        ci0,ci1=ci
        for j,pts in enumerate(zip(ci0[:-1],ci0[1:],ci1[1:],ci1[:-1],ci0[:-1])):
            p=geometry.Polygon(pts)
            if p.intersects(geom):
                indices.append([i,j])
                cells.append(p)
    return np.array(indices),cells

def gridPointsInRange( grid, geom, buffer ):
    bgeom=geom.buffer(buffer)
    inrange=np.zeros((grid.shape[0],grid.shape[1],1))
    for i,col in enumerate(grid):
        for j,pt in enumerate(col):
            if geometry.Point(pt).intersects(bgeom):
                inrange[i,j]=1
    return inrange


    



if __name__ == "__main__":
    main()
