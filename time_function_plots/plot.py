import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from TimeFunction import BaseTimeFunction
from Constants import *
import yaml
import argparse

parser=argparse.ArgumentParser(description="Plot time functions defined in YAML file to SVG, PNG")
parser.add_argument("yaml_file",help="YAML file defining plots and time functions")
parser.add_argument("-p","--png",action="store_true",help="Create PNG plot instead of SVG")
parser.add_argument("-b","--plot-base-name",default="",help="prefix for names of plot files")
parser.add_argument("-a", "--annotate", action="store_true", help="Add reference date annotations")
args=parser.parse_args()


allplots=yaml.load(open(args.yaml_file).read(),Loader=yaml.Loader)

minepoch=float(allplots.get("min_epoch","2012.0"))
maxepoch=float(allplots.get("max_epoch","2012.0"))
nepoch=int(allplots.get("intervals",400))
plotepochs=np.linspace(2012.0,2016.0,400,True)

format="png" if args.png else "svg"
fileprefix=args.plot_base_name
filesuffix="."+format

annotate=args.annotate or allplots.get('annotate', False)

epochlabels={
    TIME_PARAM_START_EPOCH: '$t_s$',
    TIME_PARAM_END_EPOCH: '$t_e$',
    TIME_PARAM_FUNCTION_REFERENCE_EPOCH: '$t_0$',
    TIME_PARAM_EVENT_EPOCH: '$t_v$',
}

plots=allplots.get("plots")
for plot in plots:
    plotdef = plots[plot]
    plotfile=fileprefix+plot+filesuffix

    timefunc=plotdef["time_function"]
    if not isinstance(timefunc,list):
        timefunc=[timefunc]

    try:
        xv=plotepochs.copy()
        functions=[]
        annotations=set()
        for fdef in timefunc:
            functions.append(BaseTimeFunction.Create(fdef))
            if TIME_PARAM_EVENT_EPOCH in fdef:
                epoch=float(fdef[TIME_PARAM_EVENT_EPOCH])
                np.append(xv,[epoch-0.0001,epoch+0.0001])
            if annotate:
                for param, label in epochlabels.items():
                    if param in fdef:
                        epoch=float(fdef[param])
                        annotations.add((epoch,label))
    except Exception as ex:
        print(f"Failed to load {plot}: {ex}")
        continue

    xv=np.sort(xv)
    yv=xv*0.0
    for func in functions:
        yv += [func.valueAt(e) for e in xv]

    title=plotdef.get("title","")

    fig1 = plt.figure()
    ax = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_xlim(minepoch, maxepoch)
    ax.set_ylim(-1.2, 1.2)
    plt.axhline(0.0,c='gray')

    if title:
        plt.title(title,fontsize=24)
    
    if annotations:
        epochs=set((a[0] for a in annotations))

        for epoch in sorted(epochs):
            text=', '.join(sorted(set((a[1] for a in annotations if a[0]==epoch))))
            ax.annotate(text,xy=(epoch,1.0),ha='center',va='bottom',fontsize=18)
            plt.axvline(epoch,0.0,2.2/2.4,c='gray')

    l1, = ax.plot(xv,yv, "b-",lw=5)

    # save the figure as a bytes string in the svg format.
    with open(plotfile,'wb') as f:
        plt.savefig(f, format=format)
