"""Six representative slices of the three-parameter steady COMSOL matrix.

Groups (a-c), one per row: T-strain, concentration-T, concentration-strain.
Each group contains two heatmaps. Columns: low and
high values of the third parameter. All panels share one linear viridis
color scale for R=c(0)/c(L), displayed horizontally below the six heatmaps.
Cells are actual cases, without interpolation.
The unequal concentration levels are categorical, equally spaced cells.

Suggested caption: Joint effects of temperature, strain and initial
concentration on the steady endpoint concentration ratio. Group (a) shows
T-strain slices at c0/c_ref=0.04 (left) and 1.00 (right); group (b) shows
concentration-T slices at absolute strain=5% (left) and 8% (right);
group (c) shows concentration-strain slices at T=450 K (left) and 600 K
(right). Each cell represents one
COMSOL simulation. The shared color scale and annotated ratios permit
comparisons across all six slices; intermediate fixed levels are omitted.

Redraw the included numerical data without rerunning COMSOL.
Exports a 600-dpi PNG and a fully vector PDF with embedded TrueType fonts.
Intended for full-width placement, approximately 170 mm wide. The cell
numbers then remain approximately 11 pt; single-column placement is too
small for six heatmaps. Reference notebooks use Times New Roman, larger
tick labels, inward ticks and heavier axes; no reference code is executed.
"""
from pathlib import Path
import argparse
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import FixedLocator
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SNAPSHOT=ROOT/'data'/'figures'/'figure9_data.csv'
LEVELS={'temperature_K':(450,500,550,600),'abs_strain':(.05,.06,.07,.08),'c0_ratio':(.04,.16,.4,1.)}
PANELS=(
    ('a','c0_ratio',.04,'temperature_K','abs_strain',r'$c_0/c_{\rm ref}=0.04$'),
    ('b','c0_ratio',1.,'temperature_K','abs_strain',r'$c_0/c_{\rm ref}=1.00$'),
    ('c','abs_strain',.05,'c0_ratio','temperature_K',r'$|\varepsilon_0|=5\%$'),
    ('d','abs_strain',.08,'c0_ratio','temperature_K',r'$|\varepsilon_0|=8\%$'),
    ('e','temperature_K',450,'c0_ratio','abs_strain',r'$T=450\,\mathrm{K}$'),
    ('f','temperature_K',600,'c0_ratio','abs_strain',r'$T=600\,\mathrm{K}$'),
)
LABELS={'temperature_K':r'$T$ (K)','abs_strain':r'$|\varepsilon_0|$ (%)','c0_ratio':r'$c_0/c_{\rm ref}$'}
# Large vector master, comparable to the supplied thesis reference PDFs.
# Scaling all physical dimensions together preserves the checked layout.
MASTER_SCALE=2.0

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))

def ticklabels(key):
    if key=='abs_strain':return [f'{100*v:.0f}' for v in LEVELS[key]]
    if key=='c0_ratio':return [f'{v:.2f}' for v in LEVELS[key]]
    return [f'{v:.0f}' for v in LEVELS[key]]

def draw(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    data=read(SNAPSHOT);assert len(data)==96
    assert {r['panel'] for r in data}==set('abcdef')
    norm=Normalize(float(data[0]['color_min']),float(data[0]['color_max']))
    cmap=matplotlib.colormaps['viridis']
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],'mathtext.fontset':'stix',
        'axes.unicode_minus':False,'axes.linewidth':1.65*MASTER_SCALE,'font.size':13.5*MASTER_SCALE,
        'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42,'ps.fonttype':42})
    # Preserve heatmap height while adding row separation and a bottom colorbar.
    fig,axes=plt.subplots(3,2,figsize=(8.8*MASTER_SCALE,11.5*MASTER_SCALE))
    for index,(ax,(code,fixed,value,xkey,ykey,title)) in enumerate(zip(axes.flat,PANELS)):
        cells=sorted([r for r in data if r['panel']==code],key=lambda r:(int(r['cell_row']),int(r['cell_column'])))
        assert len(cells)==16 and len({r['case'] for r in cells})==16
        assert all(float(r[fixed])==value and r['x_parameter']==xkey and r['y_parameter']==ykey for r in cells)
        assert all(r['panel_group']=='abc'[index//2] and r['slice_position']==('left' if index%2==0 else 'right') for r in cells)
        values=np.array([float(r['endpoint_ratio']) for r in cells]).reshape(4,4)
        # pcolormesh writes each discrete cell as a vector polygon in PDF.
        edges=np.arange(5)-.5
        plot=ax.pcolormesh(edges,edges,values,cmap=cmap,norm=norm,shading='flat',rasterized=False)
        ax.set(xlim=(-.5,3.5),ylim=(-.5,3.5),aspect='auto')
        ax.set(xticks=range(4),xticklabels=ticklabels(xkey),yticks=range(4),yticklabels=ticklabels(ykey))
        ax.set_xlabel(LABELS[xkey],fontsize=18*MASTER_SCALE,labelpad=3*MASTER_SCALE)
        ax.set_ylabel(LABELS[ykey],fontsize=18*MASTER_SCALE,labelpad=4*MASTER_SCALE)
        ax.set_title(title,fontsize=17*MASTER_SCALE,pad=8*MASTER_SCALE)
        ax.tick_params(which='major',length=4.5*MASTER_SCALE,width=1.35*MASTER_SCALE,labelsize=13.5*MASTER_SCALE,pad=3*MASTER_SCALE)
        ax.xaxis.set_minor_locator(FixedLocator(np.arange(-.5,4,1)))
        ax.yaxis.set_minor_locator(FixedLocator(np.arange(-.5,4,1)))
        ax.grid(which='minor',color='white',lw=.95*MASTER_SCALE);ax.tick_params(which='minor',length=0)
        if index%2==0:
            ax.text(-.20,1.035,f'({"abc"[index//2]})',transform=ax.transAxes,
                    fontsize=18*MASTER_SCALE,ha='left',va='bottom',clip_on=False)
        for j in range(4):
            for i in range(4):
                rgb=np.array(cmap(norm(values[j,i]))[:3])
                linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
                color='#171717' if float(linear@np.array([.2126,.7152,.0722]))>.179 else 'white'
                ax.text(i,j,f'{values[j,i]:.3f}',ha='center',va='center',fontsize=13*MASTER_SCALE,color=color)
    fig.subplots_adjust(left=.105,right=.945,bottom=.165,top=.955,hspace=.46,wspace=.31)
    bar_ax=fig.add_axes([.240,.088,.570,.016])
    bar=fig.colorbar(plot,cax=bar_ax,orientation='horizontal',ticks=[1.21,2,3,4,5,6],format='%.2f')
    bar.solids.set_rasterized(False)
    bar.solids.set_edgecolor('face')
    bar.set_label(r'$c(0)/c(L)$',fontsize=18*MASTER_SCALE,labelpad=8*MASTER_SCALE)
    bar.ax.tick_params(direction='in',length=4.5*MASTER_SCALE,width=1.35*MASTER_SCALE,labelsize=13.5*MASTER_SCALE)
    fig.savefig(output_dir/'figure9.png',dpi=600,bbox_inches='tight',facecolor='white')
    fig.savefig(output_dir/'figure9.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)
    print(f'Figure 9: 3 labelled groups, 6 heatmaps, 96 displayed cells, {len({r["case"] for r in data})} unique cases; shared scale {norm.vmin:.6f} to {norm.vmax:.6f}')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'reproduced_figures')
    args=parser.parse_args()
    draw(args.output_dir)
