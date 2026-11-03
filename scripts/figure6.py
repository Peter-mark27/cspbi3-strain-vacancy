"""Baseline electrostatic-feedback formation at four saved times.

Suggested caption: Spatial charge density rho (a), potential phi (b), and
longitudinal electric field E_x (c) at t=0, 0.253, 1.992 and 9.203 ms for
T=600 K, |epsilon0|=8%, c0=2.002 mol/m3 and L=300 nm. All 460 y=0 cut-line
nodes are retained from the same baseline simulation as Figure 5.
rho=F(zV*c-cb), with cb=zV*c0; phi obeys the Poisson problem with its
dielectric coefficient and boundary conditions; E_x=-partial_x phi.
All three fields are direct COMSOL stored exports, not reconstructed from
smoothed concentrations or numerically differentiated plotting curves.
Small initialization and endpoint export residuals are preserved.

rho is an affine transformation of the Figure 5 concentration at fixed c0;
the purpose here is to show electrostatic-field formation, not to claim an
independent validation from rho alone. Input is the packaged frozen CSV,
which is read without modification.
The layout uses three vertically stacked axes of approximately
1.6:1 width-to-height ratio and a 12.8-inch PDF master width.
The shared four-time legend appears only in panel (a); all panels retain
the same time colors, spatial limits and source data.
"""
from __future__ import annotations
import argparse,csv,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
from matplotlib.transforms import Bbox

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'figures'/'figure6_data.csv'
DEFAULT_OUTPUT=ROOT/'reproduced_figures'
PALETTE=('#3C5488','#00A087','#E69F00','#E64B35')
WIDTH=12.80;PNG_WIDTH=(7680+.25)/600;HANDLE_PT=(3*5+2*2.6)*2.88
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def table(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def build(rows):
    plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif','DejaVu Serif'],'mathtext.fontset':'stix','axes.unicode_minus':False,
        'axes.linewidth':3.2,'axes.labelsize':36,'xtick.labelsize':28,'ytick.labelsize':28,
        'legend.fontsize':23,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':600})
    fig=plt.figure(figsize=(12.9,26.3))
    axes=[fig.add_axes([.145,bottom,.835,.254]) for bottom in (.715,.392,.069)]
    for panel,ax in zip('abc',axes):
        ax.tick_params(which='major',direction='in',length=8,width=2.4,pad=8)
        ax.tick_params(which='minor',direction='in',length=5,width=1.8)
        ax.xaxis.set_minor_locator(AutoMinorLocator(2));ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.text(-.065,1.018,f'({panel})',transform=ax.transAxes,fontsize=31,ha='left',va='bottom')
        for series,color in enumerate(PALETTE):
            rr=sorted((r for r in rows if r['panel']==panel and int(r['series'])==series),key=lambda r:int(r['node']))
            t=float(rr[0]['time_s'])*1e3;label=r'$t=0$' if series==0 else rf'$t={t:.3f}\,\mathrm{{ms}}$'
            scale={'a':1e-3,'b':1.,'c':1e-6}[panel]
            ax.plot([float(r['x_SI'])*1e9 for r in rr],[float(r['value_SI'])*scale for r in rr],color=color,lw=3.52,label=label)
        ax.set_xlim(-6,306);ax.xaxis.set_major_locator(MultipleLocator(50));ax.set_xlabel(r'$x$ (nm)',labelpad=6)
        ax.set_ylabel({'a':r'$\rho$ (kC m$^{-3}$)','b':r'$\phi$ (V)','c':r'$E_x$ (MV m$^{-1}$)'}[panel],fontsize=32,labelpad=8)
        ax.set_ylim({'a':(-34,35),'b':(-.365,.012),'c':(-.055,1.34)}[panel])
        ax.yaxis.set_major_locator(MultipleLocator({'a':10,'b':.05,'c':.25}[panel]))
        if panel=='a':
            ax.legend(loc='upper right',bbox_to_anchor=(.97,.965),ncol=1,frameon=False,
                handlelength=HANDLE_PT/23,handletextpad=.35,labelspacing=.25,borderaxespad=0,borderpad=0)
    return fig
def boxes(fig):
    pad=.24
    for _ in range(8):
        fig.canvas.draw();b=fig.get_tightbbox(fig.canvas.get_renderer());correction=WIDTH-b.width-2*pad
        if abs(correction)<1e-7:break
        fig.set_size_inches(fig.get_figwidth()+correction/.92,fig.get_figheight(),forward=True)
    fig.canvas.draw();b=fig.get_tightbbox(fig.canvas.get_renderer())
    return Bbox.from_bounds(b.x0-pad,b.y0-pad,PNG_WIDTH,b.height+2*pad),Bbox.from_bounds(b.x0-pad,b.y0-pad,WIDTH,b.height+2*pad)
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,default=DEFAULT_OUTPUT);args=p.parse_args()
    rows=table(DATA);frozen=sha(DATA);fig=build(rows);pngbox,pdfbox=boxes(fig);args.output_dir.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output_dir/'figure6.pdf',bbox_inches=pdfbox,facecolor='white',metadata={'Title':'Baseline electrostatic-feedback formation'})
    fig.savefig(args.output_dir/'figure6.png',dpi=600,bbox_inches=pngbox,facecolor='white');plt.close(fig)
    assert sha(DATA)==frozen;print('Figure6: three electrostatic fields at four states,',len(rows),'rows; CSV SHA256',frozen)
if __name__=='__main__':main()
