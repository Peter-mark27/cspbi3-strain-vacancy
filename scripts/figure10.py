"""COMSOL dimensionless response and short-screening-length prediction.

(a) Actual ratios at all 64 simulated Pi,delta pairs, with constant-alpha
geometric guides. There is no filled theoretical surface or interpolation.
(b) Data and one theoretical curve R=exp(2*alpha).
The finite-length approximation is retained in the case data for comparison;
short_screening_error_percent is the absolute relative error against COMSOL.

Suggested caption: Dimensionless organization of the steady endpoint ratio.
(a) COMSOL results in the Pi-delta parameter space; dashed lines indicate
constant alpha=Pi*delta. The delta axis is logarithmic and colors denote
actual endpoint concentration ratios. (b) COMSOL ratios versus alpha,
compared with the short-screening-length prediction (one dashed curve).
The maximum relative deviation from COMSOL is approximately 3.69%.

Redraw the included numerical data without rerunning COMSOL.
Exports a 600-dpi PNG and a fully vector PDF with embedded TrueType fonts.
"""
from pathlib import Path
import argparse
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator,NullFormatter
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SNAPSHOT=ROOT/'data'/'figures'/'figure10_data.csv'
MASTER_SCALE=2.0
RATIO_LABEL_SIZE=16*MASTER_SCALE
# Preserve the original panel proportions with modestly larger typography.
# Scale geometry and physical style together for a large vector master.

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as stream:return list(csv.DictReader(stream))

def style(ax,code):
    ax.tick_params(which='major',length=4.5*MASTER_SCALE,width=1.2*MASTER_SCALE,
                   labelsize=12.5*MASTER_SCALE,pad=3*MASTER_SCALE)
    ax.tick_params(which='minor',length=2.2*MASTER_SCALE,width=.85*MASTER_SCALE)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.text(-.16,1.035,f'({code})',transform=ax.transAxes,fontsize=15.5*MASTER_SCALE,
            ha='left',va='bottom',clip_on=False)

def draw(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    data=read(SNAPSHOT)
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],'mathtext.fontset':'stix',
        'axes.unicode_minus':False,'axes.linewidth':1.3*MASTER_SCALE,'font.size':12.5*MASTER_SCALE,
        'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42,'ps.fonttype':42})
    fig=plt.figure(figsize=(10.5*MASTER_SCALE,4.35*MASTER_SCALE))
    left=fig.add_axes([.078,.17,.342,.775]);right=fig.add_axes([.603,.17,.375,.775])
    points=[r for r in data if r['panel']=='a' and r['type']=='COMSOL']
    assert len(points)==64 and len({r['case'] for r in points})==64
    pi=np.array([float(r['Pi']) for r in points]);delta=np.array([float(r['delta']) for r in points])
    norm=Normalize(float(points[0]['color_min']),float(points[0]['color_max']))
    left.set(xlim=(pi.min()*.96,pi.max()*1.04),ylim=(delta.min()*.94,delta.max()*1.06))
    left.set_yscale('log');left.set_yticks([.02,.04,.06,.08,.1],labels=['0.02','0.04','0.06','0.08','0.10'])
    left.set_xticks([5,6,7,8,9])
    left.yaxis.set_minor_formatter(NullFormatter())
    pg=np.linspace(*left.get_xlim(),500)
    for a in (.2,.4,.6,.8):
        dg=a/pg;take=(dg>=left.get_ylim()[0])&(dg<=left.get_ylim()[1])
        left.plot(pg[take],dg[take],color='#909090',ls=(0,(3,2.5)),lw=1.0*MASTER_SCALE,zorder=1)
        label_pi=8.7
        left.text(label_pi,a/label_pi,rf'$\alpha={a:.1f}$',color='#666666',fontsize=10.5*MASTER_SCALE,
                  ha='center',va='center',bbox=dict(facecolor='white',edgecolor='none',pad=.8*MASTER_SCALE),zorder=4)
    plot=left.scatter(pi,delta,c=[float(r['endpoint_ratio']) for r in points],cmap='viridis',norm=norm,
                      s=81*MASTER_SCALE**2,edgecolors='white',linewidths=.7*MASTER_SCALE,zorder=3)
    left.set_xlabel(r'$\Pi$',fontsize=18*MASTER_SCALE,labelpad=3*MASTER_SCALE)
    left.set_ylabel(r'$\delta$',fontsize=20*MASTER_SCALE,labelpad=4*MASTER_SCALE)
    style(left,'a')
    # Original major-tick density; align end labels inward without changing normalization.
    bar_ticks=[1.21,2,3,4,5,6]
    bar_ax=fig.add_axes([.443,.17,.017,.775]);bar=fig.colorbar(plot,cax=bar_ax,ticks=bar_ticks,format='%.2f')
    bar.solids.set_rasterized(False)
    bar.solids.set_edgecolor('face')
    bar.set_label(r'$c(0)/c(L)$',fontsize=RATIO_LABEL_SIZE,labelpad=5*MASTER_SCALE)
    bar.ax.tick_params(direction='in',length=3.5*MASTER_SCALE,width=1.1*MASTER_SCALE,labelsize=11.5*MASTER_SCALE)
    bar.ax.get_yticklabels()[0].set_verticalalignment('bottom')
    bar.ax.get_yticklabels()[-1].set_verticalalignment('top')
    # At the upper limit the frame already supplies the tick position.
    for line in (bar.ax.yaxis.get_major_ticks()[-1].tick1line,bar.ax.yaxis.get_major_ticks()[-1].tick2line):
        line.set_visible(False)
    actual=[r for r in data if r['panel']=='b' and r['type']=='COMSOL']
    assert len(actual)==64 and {r['case'] for r in actual}=={r['case'] for r in points}
    curve=[r for r in data if r['panel']=='b' and r['type']=='short_screening']
    assert len(curve)==401
    curve_alpha=np.array([float(r['alpha']) for r in curve])
    curve_ratio=np.array([float(r['short_screening_ratio']) for r in curve])
    assert np.allclose(curve_ratio,np.exp(2*curve_alpha),rtol=1e-12,atol=0)
    right.plot(curve_alpha,curve_ratio,
               color='#4D4D4D',ls=(0,(3.2,2.2)),lw=1.9*MASTER_SCALE,zorder=1)
    alpha=np.array([float(r['alpha']) for r in actual]);ratios=np.array([float(r['endpoint_ratio']) for r in actual])
    right.scatter(alpha,ratios,s=81*MASTER_SCALE**2,color='#3C5488',edgecolors='white',linewidths=.7*MASTER_SCALE,zorder=3)
    xspan=float(alpha.max()-alpha.min());ymin=min(float(r['linear_pb_ratio']) for r in actual);ymax=float(ratios.max());yspan=ymax-ymin
    right.set(xlim=(alpha.min()-.05*xspan,alpha.max()+.05*xspan),ylim=(ymin-.055*yspan,ymax+.055*yspan))
    right.set_xticks([.2,.4,.6,.8])
    right.set_xlabel(r'$\alpha=\Pi\delta$',fontsize=18*MASTER_SCALE,labelpad=3*MASTER_SCALE)
    right.set_ylabel(r'$c(0)/c(L)$',fontsize=RATIO_LABEL_SIZE,labelpad=4*MASTER_SCALE)
    right.yaxis.set_minor_locator(AutoMinorLocator(2));style(right,'b')
    handles=[Line2D([],[],marker='o',ls='None',markersize=9*MASTER_SCALE,markerfacecolor='#3C5488',markeredgecolor='white',label='Finite-element solution'),
             Line2D([],[],color='#4D4D4D',ls=(0,(2.4,1.45)),lw=1.9*MASTER_SCALE,
                    label='Analytical solution')]
    right.legend(handles=handles,loc='upper left',fontsize=11*MASTER_SCALE,frameon=False,
                 labelspacing=.65,handlelength=1.75,handletextpad=.38,borderaxespad=.35)
    fig.savefig(output_dir/'figure10.png',dpi=600,bbox_inches='tight',pad_inches=.14*MASTER_SCALE,facecolor='white')
    fig.savefig(output_dir/'figure10.pdf',bbox_inches='tight',pad_inches=.14*MASTER_SCALE,facecolor='white');plt.close(fig)
    error=max(float(r['short_screening_error_percent']) for r in actual)
    print(f'Figure 10: 2 panels, 64 cases in each; one R=exp(2*alpha) curve; max deviation relative to COMSOL {error:.6f}%')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'reproduced_figures')
    args=parser.parse_args()
    draw(args.output_dir)
