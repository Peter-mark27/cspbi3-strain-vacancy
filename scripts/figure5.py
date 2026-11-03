"""Baseline concentration evolution: spatial snapshots and left-end history.

(a) Original four COMSOL concentration profiles; all 460 y=0 nodes retained.
(b) Direct stored c/c0 at x=0, y=0, shown over 0--12 ms. Colored markers
identify the four times in (a) and connect by arrows to native 2D concentration maps. Endpoint location is stated in the caption,
not annotated inside the panel. Both y axes are labeled c/c0.

One frozen CSV retains the four spatial profiles (series 0--3) and all 800
left-end history states (series 4), including the full-profile residual,
and native 2D concentration fields (series 5--8; 70840 nodes per field).
No temporal resampling or forced initial values is used. Only the inset
images use spatial linear interpolation of the native 2D domain nodes.
All insets share a c/c0 scale of 0.84--1.16 and retain 3:1 geometry.
Viridis matches the supplied reference style; no colorbar is displayed.
Exact time metadata are used where available; remaining times retain the
COMSOL concentration-export header precision. The selected settling state
is verified using the persistent relative Euclidean profile criterion 1e-4.

Redraws the packaged frozen CSV without modifying it.
Exports a vector PDF and 600 dpi PNG matching Figure 10's output dimensions.

Suggested caption: Evolution of the normalized vacancy concentration along
y=0 at T=600 K, epsilon0=-0.08, and c0=2.002 mol/m3. (a) Spatial profiles at
the indicated times. (b) Concentration at the compressed end (x=0) as a
function of time. Colored markers and arrows link the curve to concentration maps at the times in (a); all maps share the same color scale.
"""
from __future__ import annotations
import argparse,csv,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator,FormatStrFormatter
from matplotlib.tri import Triangulation, LinearTriInterpolator
from matplotlib.colors import Normalize

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'figures'/'figure5_data.csv'
DEFAULT_OUTPUT=ROOT/'reproduced_figures'
# Measured from Figure 10's existing PDF MediaBox and 600-dpi PNG.
PDF_SIZE_IN=(1491.63953125/72,638.491803125/72)
PNG_SIZE_PX=(12429,5320)
PALETTE=('#3C5488','#00A087','#E69F00','#E64B35')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def table(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def add_clouds(ax, rows, history, clouds):
    # Actual 2D native-node results, interpolated spatially for rendering only.
    groups=[sorted((r for r in clouds if int(r['series'])==5+i),key=lambda r:int(r['node'])) for i in range(4)]
    xy=np.array([[float(r['x_m'])*1e9,float(r['y_m'])*1e9] for r in groups[0]])
    tri=Triangulation(xy[:,0],xy[:,1])
    xx,yy=np.meshgrid(np.linspace(0,300,901),np.linspace(0,100,301))
    norm=Normalize(vmin=.84,vmax=1.16)
    # Preserve the 300 nm : 100 nm specimen aspect ratio in all inset maps.
    box=ax.get_position();fig=ax.figure
    height=.40*(box.width*fig.get_figwidth())/(box.height*fig.get_figheight())/3
    early=next(r for r in rows if int(r['series'])==1 and int(r['node'])==0)
    ymin,ymax=ax.get_ylim()
    early_y=(float(early['c_over_c0'])-ymin)/(ymax-ymin)
    initial=next(r for r in rows if int(r['series'])==0 and int(r['node'])==0)
    initial_bottom=(float(initial['c_over_c0'])-ymin)/(ymax-ymin)
    early_bottom=(initial_bottom+.59)/2
    late=next(r for r in rows if int(r['series'])==2 and int(r['node'])==0)
    steady=next(r for r in rows if int(r['series'])==3 and int(r['node'])==0)
    xmin,xmax=ax.get_xlim()
    late_x=(float(late['time_s'])*1e3-xmin)/(xmax-xmin)
    steady_x=(float(steady['time_s'])*1e3-xmin)/(xmax-xmin)
    # Lower only the initial map; keep the horizontal arrow at its data level.
    anchors=[(.12,initial_bottom-.0075),(.14,early_bottom),(.15,.59),(.58,.65)]
    targets=[(.12,initial_bottom),(.14,early_y),(late_x,.59+height),(steady_x,.65+height)]
    for i,(group,anchor,color) in enumerate(zip(groups,anchors,PALETTE)):
        thisxy=np.array([[float(r['x_m'])*1e9,float(r['y_m'])*1e9] for r in group])
        np.testing.assert_array_equal(thisxy,xy)
        values=np.array([float(r['c_over_c0']) for r in group])
        image=LinearTriInterpolator(tri,values)(xx,yy)
        assert np.isfinite(image).all()
        inset=ax.inset_axes([*anchor,.40,height],zorder=4)
        im=inset.imshow(image,origin='lower',extent=(0,300,0,100),cmap='viridis',norm=norm,aspect='auto',interpolation='bilinear')
        inset.set_xticks([]);inset.set_yticks([])
        for spine in inset.spines.values():spine.set_linewidth(.8);spine.set_edgecolor('#777777')
        ref=next(r for r in rows if int(r['series'])==i and int(r['node'])==0)
        point=(float(ref['time_s'])*1e3,float(ref['c_over_c0']))
        ax.annotate('',xy=targets[i],xycoords='axes fraction',xytext=point,textcoords='data',
                    arrowprops=dict(arrowstyle='->',color=color,lw=2.4,mutation_scale=20,shrinkA=7.5,shrinkB=0),zorder=5)
    # Labels match the mechanism schematic; annotate only the steady-state map.
    left,bottom=anchors[3]
    ax.text(left,bottom-.018,'Accumulation',transform=ax.transAxes,
            ha='left',va='top',fontsize=20,color='#333333',zorder=5)
    ax.text(left+.40,bottom-.018,'Depletion',transform=ax.transAxes,
            ha='right',va='top',fontsize=20,color='#333333',zorder=5)

def build(rows,history,clouds):
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif','DejaVu Serif'],'mathtext.fontset':'stix','axes.unicode_minus':False,
        'axes.linewidth':2.6,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':600})
    fig,axes=plt.subplots(1,2,figsize=PDF_SIZE_IN)
    ax=axes[0]
    for series,color in enumerate(PALETTE):
        rr=sorted((r for r in rows if int(r['series'])==series),key=lambda r:int(r['node']))
        t=float(rr[0]['time_s'])*1e3;label=r'$t=0$' if series==0 else rf'$t={t:.3f}\,\mathrm{{ms}}$'
        ax.plot([float(r['x_m'])*1e9 for r in rr],[float(r['c_over_c0']) for r in rr],color=color,lw=4.0,label=label)
    ax.set_xlim(-10,310);ax.set_ylim(.84,1.18)
    assert all(.84<float(r['c_over_c0'])<1.18 for r in rows)
    ax.set_xlabel(r'$x$ (nm)',fontsize=36,labelpad=6)
    ax.xaxis.set_major_locator(MultipleLocator(50));ax.yaxis.set_major_locator(MultipleLocator(.05))
    ax.legend(loc='upper right',bbox_to_anchor=(.965,.965),fontsize=26,ncol=1,frameon=False,
        handlelength=36.36/17,handletextpad=.35,labelspacing=.25,borderaxespad=0,borderpad=0)
    ax=axes[1]
    time=np.array([float(r['time_s'])*1e3 for r in history])
    concentration=np.array([float(r['c_over_c0']) for r in history])
    ax.plot(time,concentration,color='#555555',lw=4.0,zorder=2)
    for series,color in enumerate(PALETTE):
        r=next(r for r in rows if int(r['series'])==series and int(r['node'])==0)
        k=int(r['index_zero_based'])
        assert np.isclose(concentration[k],float(r['c_over_c0']),rtol=1e-9,atol=1e-10)
        ax.plot(time[k],concentration[k],'o',color=color,ms=14,mec='white',mew=1.2,zorder=3)
    ax.set_xlim(-.4,12);ax.set_ylim(.99,1.18)
    ax.set_xlabel(r'$t$ (ms)',fontsize=36,labelpad=6)
    ax.xaxis.set_major_locator(MultipleLocator(2));ax.yaxis.set_major_locator(MultipleLocator(.04))
    for code,ax in zip('ab',axes):
        ax.set_ylabel(r'$c/c_0$',fontsize=36,labelpad=6)
        ax.tick_params(which='major',direction='in',labelsize=25,width=2.4,length=9,pad=6)
        ax.tick_params(which='minor',direction='in',width=1.7,length=4.4)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2));ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.text(-.16,1.035,f'({code})',transform=ax.transAxes,fontsize=31,ha='left',va='bottom')
    fig.subplots_adjust(left=.095,right=.98,bottom=.165,top=.925,wspace=.30)
    add_clouds(axes[1],rows,history,clouds)
    return fig
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,default=DEFAULT_OUTPUT);args=p.parse_args()
    data=table(DATA)
    rows=[r for r in data if int(r['series'])<4]
    history=[r for r in data if int(r['series'])==4]
    clouds=[r for r in data if int(r['series'])>=5]
    assert len(clouds)==4*70840
    assert len(rows)==1840 and len(history)==800
    history=sorted(history,key=lambda r:int(r['index_zero_based']))
    frozen=sha(DATA)
    fig=build(rows,history,clouds);args.output_dir.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output_dir/'figure5.pdf',facecolor='white',metadata={'Title':'Baseline vacancy concentration spatial and temporal evolution'})
    # Small subpixel allowance avoids truncation to one pixel below the target.
    fig.set_size_inches(*((n+.25)/600 for n in PNG_SIZE_PX))
    fig.savefig(args.output_dir/'figure5.png',dpi=600,facecolor='white')
    plt.close(fig)
    assert sha(DATA)==frozen
    print('Figure5:',len(rows),'spatial rows and',len(history),'history rows; combined CSV SHA256',frozen)
if __name__=='__main__':main()
