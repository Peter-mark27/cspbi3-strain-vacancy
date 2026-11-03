"""Five-panel baseline flux competition figure from direct COMSOL exports.

Panels (a-d): complete y=0 baseline fluxes at initial time, saved times nearest
0.25/2.00 ms, and the first persistently converged state. Panel (e): spatial mean-absolute
flux shares up to the first persistently converged concentration profile.
No endpoints are discarded, no flux signs are reversed, and no time states
are interpolated. Small exported initialization/terminal residuals remain.

Suggested caption: Evolution of strain-driven transport and its opposing
electrostatic and diffusive fluxes at T=600 K, |epsilon0|=8%, c0=2.002 mol/m3,
L=300 nm. (a-d) Signed spatial flux components and their directly exported
sum J_x at t=0, 0.253 ms, 1.992 ms, and t_ss=9.203 ms. (e) Spatially averaged
absolute component shares, normalized by their sum. t_ss=9.203 ms is the
first saved time after which the discrete normalized full-profile residual
stays <=1e-4. Component shares are magnitudes, not signed net transport or
fractions of transferred vacancies. Steady component balance is demonstrated
by the signed spatial curves and small net-flux residual in panel (d).

Uses the packaged frozen CSV without modifying it.
Connecting curves retain all saved spatial nodes/time points.
All panels use the original benchmark MPH and the same y=0 cut line (460
nodes). Component lines are solid; total J_x is thin dark gray dashed.
Legends use a single column: upper left in (a-d), upper right in (e).
Spatial samples have three
full gray dashes. Titles give only the actual snapshot times, including
t=9.203 ms for panel (d). Baseline conditions and the steady criterion are
given in the caption rather than as a global plot title.
All five plotting frames have equal width and height. Snapshot titles sit
closer to the upper frame; moderate vertical separation keeps adjacent
labels comfortably apart without changing the equal plotting-frame sizes.
No horizontal zero reference line is drawn. J_x corresponds to J in the
one-dimensional manuscript notation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FuncFormatter
from matplotlib.transforms import Bbox
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data'/'figures'/'figure7_data.csv'
DEFAULT_OUTPUT = ROOT/'reproduced_figures'
WIDTH = 20.717215711805554
EXPORT_HEIGHT = 21.142589301215285  # Preserve the original canvas during label edits.
PNG_WIDTH = (12429+.25)/600
Q = 2.0
COLORS = {'Jdiff':'#3C5488', 'Jstrain':'#E64B35', 'Jelec':'#00A087', 'Jx':'#555555'}
NAMES = ('Jdiff','Jstrain','Jelec','Jx')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read_csv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as stream:
        return list(csv.DictReader(stream))

def values(rows,panel,name):
    group=sorted((r for r in rows if r['panel']==panel and r['component']==name),key=lambda r:int(r['sample']))
    assert group
    return np.array([float(r['x_SI']) for r in group]),np.array([float(r['value_SI']) for r in group]),group[0]

def style():
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif','DejaVu Serif'],
        'mathtext.fontset':'stix','pdf.fonttype':42,'ps.fonttype':42,
        'axes.labelsize':18*Q,'xtick.labelsize':14*Q,'ytick.labelsize':14*Q,
        'legend.fontsize':11.5*Q,'axes.linewidth':1.6*Q,'axes.unicode_minus':True,
        'axes.grid':False,'savefig.dpi':600,'lines.dash_capstyle':'butt'})

def build(rows):
    style()
    fig=plt.figure(figsize=(21,21.855))
    gs=fig.add_gridspec(3,4,left=.095,right=.985,bottom=.070,top=.945,
        wspace=.90,hspace=.36,height_ratios=[1,1,1])
    axes=[fig.add_subplot(gs[0,:2]),fig.add_subplot(gs[0,2:]),
          fig.add_subplot(gs[1,:2]),fig.add_subplot(gs[1,2:]),fig.add_subplot(gs[2,1:3])]
    labels={'Jdiff':r'$J_{\mathrm{diff}}$','Jstrain':r'$J_{\mathrm{strain}}$',
            'Jelec':r'$J_{\mathrm{elec}}$','Jx':r'$J$'}
    for panel,ax in zip('abcde',axes):
        ax.tick_params(which='major',direction='in',length=4*Q,width=1.2*Q,pad=4*Q)
        ax.tick_params(which='minor',direction='in',length=2.5*Q,width=.9*Q)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2));ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        label_y = 0.998 if panel == 'e' else 0.983
        ax.text(-.16,label_y,f'({panel})',transform=ax.transAxes,fontsize=17*Q,fontweight='normal',ha='left',va='bottom')
    for panel,ax in zip('abcd',axes[:4]):
        handles=[]
        for name in NAMES:
            x,y,meta=values(rows,panel,name)
            # Dash net flux so coincident strain/net curves remain distinguishable at t=0.
            line=(0,(5,2.6)) if name=='Jx' else '-'
            curve,=ax.plot(x*1e9,y,color=COLORS[name],lw=2.88 if name=='Jx' else 3.52,
                    ls=line,label=labels[name],zorder=5 if name=='Jx' else 3)
            handles.append(curve)
        ax.set_xlim(-6,306);ax.set_ylim(-8.6e-6,8.6e-6)
        ax.xaxis.set_major_locator(MultipleLocator(50));ax.yaxis.set_major_locator(MultipleLocator(2e-6))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y,_:f'{y/1e-6:g}'))
        ax.set_xlabel(r'$x$ (nm)',labelpad=3*Q)
        ax.set_ylabel(r'$J_i$ (mol m$^{-2}$ s$^{-1}$)',fontsize=16*Q,labelpad=4*Q)
        ax.text(0,1.012,r'$\times10^{-6}$',transform=ax.transAxes,ha='left',va='bottom',fontsize=12*Q,clip_on=False)
        t=float(meta['time_s'])*1e3
        title=r'$t=0$' if panel=='a' else rf'$t={t:.3f}\,\mathrm{{ms}}$'
        ax.set_title(title,fontsize=13*Q,pad=3*Q)
        # Exactly three full gray dashes in the sample, scaled with line width.
        ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(.035,.97),
                  ncol=1,frameon=False,handlelength=(3*5+2*2.6)*2.88/(11.5*Q),
                  handletextpad=.35,labelspacing=.25,borderpad=0,borderaxespad=0)
    ax=axes[4]
    for name in ('Jstrain','Jelec','Jdiff'):
        time,share,meta=values(rows,'e',name)
        ls='-'
        label={'Jstrain':r'$\chi_{\mathrm{strain}}$','Jelec':r'$\chi_{\mathrm{elec}}$',
               'Jdiff':r'$\chi_{\mathrm{diff}}$'}[name]
        ax.plot(time/float(meta['tss_s']),share,color=COLORS[name],lw=2.0*Q,ls=ls,label=label)
    ax.set_xlim(-.025,1.025);ax.set_ylim(-.025,1.04)
    ax.xaxis.set_major_locator(MultipleLocator(.2));ax.yaxis.set_major_locator(MultipleLocator(.2))
    ax.set_xlabel(r'$t/t_{\mathrm{ss}}$',labelpad=3*Q);ax.set_ylabel(r'$\chi_i(t)$',labelpad=4*Q)
    ax.legend(loc='upper right',bbox_to_anchor=(.965,.965),ncol=1,frameon=False,
              handlelength=(3*5+2*2.6)*2.88/(11.5*Q),handletextpad=.35,
              labelspacing=.25,borderpad=0,borderaxespad=0)
    return fig,axes

def crop_boxes(fig):
    pad=.12*Q
    for _ in range(8):
        fig.canvas.draw();b=fig.get_tightbbox(fig.canvas.get_renderer())
        correction=WIDTH-(b.width+2*pad)
        if abs(correction)<1e-7: break
        fig.set_size_inches(fig.get_figwidth()+correction/.92,fig.get_figheight(),forward=True)
    fig.canvas.draw();b=fig.get_tightbbox(fig.canvas.get_renderer())
    return Bbox.from_bounds(b.x0-pad,b.y0-pad,PNG_WIDTH,EXPORT_HEIGHT),Bbox.from_bounds(b.x0-pad,b.y0-pad,WIDTH,EXPORT_HEIGHT)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    rows=read_csv(DATA);frozen=sha(DATA)
    assert set(r['panel'] for r in rows)==set('abcde')
    fig,axes=build(rows);pngbox,pdfbox=crop_boxes(fig)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output_dir/'figure7.pdf',bbox_inches=pdfbox,facecolor='white',
                metadata={'Title':'Baseline spatial flux evolution and competition shares'})
    fig.savefig(args.output_dir/'figure7.png',dpi=600,bbox_inches=pngbox,facecolor='white')
    plt.close(fig);assert sha(DATA)==frozen
    print(f'Figure 7: five panels; {len(rows)} data rows; CSV SHA256 {frozen}')
    for panel in 'abcd':
        _,_,meta=values(rows,panel,'Jx');print(panel,'actual_time_s',meta['time_s'])
    print('tss_s',meta['tss_s'],'All 460 benchmark nodes retained per spatial curve. J_x uses thin gray dashes.')
    for name in ('Jstrain','Jelec','Jdiff'):
        _,y,_=values(rows,'e',name);print(name,'share at tss',y[-1])

if __name__=='__main__': main()
