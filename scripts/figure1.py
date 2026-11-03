"""Figure 1: strain-gradient loading and competing vacancy-transport mechanisms.

Run: python figure1.py
Dependencies: numpy, matplotlib.
Writes figure1.pdf and figure1.png to reproduced_figures/ by default.
Use --output-dir to select another output directory. PNG resolution is 600 dpi.

This is a qualitative schematic, not numerical data. Both panels use the
same specimen geometry and drawing scale. Vertical compression is illustrated
schematically, while the right end remains strain-free.
Vacancies are positively charged everywhere.
Counts are 8 / 16 / 5 across left / middle / right regions,
whose reference-width fractions are .16 / .68 / .16. Symbol counts
and positions are qualitative, not quantitative concentration data.
Arrows communicate directions, not quantitative spatial flux profiles.
"""
from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "reproduced_figures"
INK = "#262B34"
GRAY = "#92969B"
BLUE = "#4D79C2"
VACANCY = "#E54A6B"
PLUS = "#284A72"
STRAIN = "#E84E68"
DIFFUSION = "#00A9D7"
ELECTRIC = "#00AF69"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
    "font.size": 11,
    "mathtext.fontset": "stix", "text.color": INK,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": "white",
})

X0, X1 = 1.40, 10.30
L = X1-X0
STRAIN_VISUAL = 0.22
H = 3.60
REFERENCE_PADDING = .28
DIMENSION_GAP = .30
PANEL_LABEL_GAP_A = .30
PANEL_LABEL_GAP_B = .42


def low(x, bottom):
    q = (np.asarray(x)-X0)/L
    # Shared outline: preserve the accepted compressed left-end height and
    # expanded strain-free right-end height from panel (a).
    return bottom+H*STRAIN_VISUAL*(1-q)/2-REFERENCE_PADDING*q


def high(x, bottom):
    return bottom+H-(low(x, bottom)-bottom)


def region_vertices(left, right, bottom):
    return [(left, low(left, bottom)), (right, low(right, bottom)),
            (right, high(right, bottom)), (left, high(left, bottom))]


def text(ax,x,y,s,fs=11,ha="center",**kw):
    return ax.text(x,y,s,ha=ha,va="center",fontsize=fs,**kw)


def arrow(ax,start,end,color=INK,lw=1.0,head=7):
    obj=FancyArrowPatch(start,end,arrowstyle="-|>",mutation_scale=head,
                        color=color,linewidth=lw,shrinkA=0,shrinkB=0,zorder=5)
    ax.add_patch(obj)


def dimension(ax,start,end,word,word_xy):
    obj=FancyArrowPatch(start,end,arrowstyle="<|-|>",mutation_scale=8.0,
                        color=GRAY,linewidth=.9,shrinkA=0,shrinkB=0)
    ax.add_patch(obj)
    text(ax,*word_xy,word,fs=12.5,
         bbox={"facecolor":"white","edgecolor":"none","pad":1.0},zorder=6)


def vacancy(ax,x,y):
    r=.132
    ax.add_patch(Circle((x,y),r,facecolor=VACANCY,edgecolor="#F28DA2",lw=.45,zorder=4))
    ax.plot([x-.073,x+.073],[y,y],color=PLUS,lw=.95,zorder=5)
    ax.plot([x,x],[y-.073,y+.073],color=PLUS,lw=.95,zorder=5)


def brace(ax,a,b,y,depth=.24,end_y=None):
    m=(a+b)/2
    width=b-a
    end_y=y if end_y is None else end_y
    def point(f,d):
        return (a+f*width, y+f*(end_y-y)-d*depth)
    v=[point(0,0),
       point(0,.52),point(.10,.55),point(.25,.55),
       point(.40,.55),point(.50,.60),point(.50,1),
       point(.50,.60),point(.60,.55),point(.75,.55),
       point(.90,.55),point(1,.52),point(1,0)]
    ax.add_patch(PathPatch(MplPath(v,[MplPath.MOVETO]+[MplPath.CURVE4]*12),
                           facecolor="none",edgecolor="#456BA6",lw=.9,
                           capstyle="round",joinstyle="round"))


def panel_a(ax):
    original_bottom=1.90
    b=original_bottom-REFERENCE_PADDING
    t=original_bottom+H+REFERENCE_PADDING
    # Expand the undeformed reference vertically while retaining the
    # compressed left edge. The right edge still coincides with the
    # reference boundaries, so it remains the strain-free end.
    def lower_a(x):
        return low(x,original_bottom)
    def upper_a(x):
        return high(x,original_bottom)
    text(ax,.73,t+PANEL_LABEL_GAP_A,"(a)",fs=13.5,ha="left",fontweight="normal")
    ax.add_patch(Polygon([(X0,lower_a(X0)),(X1,lower_a(X1)),
                         (X1,upper_a(X1)),(X0,upper_a(X0))],closed=True,
                         facecolor="#B7C8E7",edgecolor=INK,lw=1.25,zorder=4))
    # Reference lines stay at the undeformed height. The solid specimen
    # outline is drawn above them, avoiding gray dashes across the corners.
    ax.plot([X0,X1],[t,t],color=GRAY,lw=.85,ls=(0,(2.676470588,2.294117647)),zorder=3)
    ax.plot([X0,X1],[b,b],color=GRAY,lw=.85,ls=(0,(2.676470588,2.294117647)),zorder=3)
    ax.plot([X0,X0],[b,lower_a(X0)],color=GRAY,lw=.65,ls=(0,(2,2)))
    ax.plot([X0,X0],[upper_a(X0),t],color=GRAY,lw=.65,ls=(0,(2,2)))
    for q in [0,.22,.44,.66]:
        x=X0+q*L
        # Begin at the undeformed top/bottom boundaries and end at the
        # deformed surface. The first pair is centered exactly at x=0.
        for start,end in [(t,upper_a(x)),(b,lower_a(x))]:
            head_length=min(.13,.65*abs(end-start))
            ax.arrow(x,start,0,end-start,width=.038,
                     head_width=.18,head_length=head_length,
                     length_includes_head=True,facecolor="#D9364E",
                     edgecolor="none",zorder=5)

    height_dimension_x=X0-DIMENSION_GAP
    ax.plot([height_dimension_x-.09,X0-.045],[t,t],color=GRAY,lw=.65)
    ax.plot([height_dimension_x-.09,X0-.045],[b,b],color=GRAY,lw=.65)
    dimension(ax,(height_dimension_x,b),(height_dimension_x,t),
              r"$H$",(height_dimension_x,(b+t)/2))
    length_dimension_y=b-DIMENSION_GAP
    ax.plot([X0,X0],[length_dimension_y-.09,b-.045],color=GRAY,lw=.65)
    ax.plot([X1,X1],[length_dimension_y-.09,b-.045],color=GRAY,lw=.65)
    dimension(ax,(X0,length_dimension_y),(X1,length_dimension_y),
              r"$L$",((X0+X1)/2,length_dimension_y))

    arrow(ax,(10.61,3.08),(11.35,3.08),color=BLUE,lw=1.1,head=7)
    arrow(ax,(10.61,3.08),(10.61,3.86),color=BLUE,lw=1.1,head=7)
    text(ax,11.47,3.08,"x",fs=11.5,color=BLUE,fontstyle="italic",fontweight="bold")
    text(ax,10.61,4.01,"y",fs=11.5,color=BLUE,fontstyle="italic",fontweight="bold")

    # Original-style left-pointing arrow, red to blue with a pale blue tail.
    cmap=LinearSegmentedColormap.from_list("strain",[
        (0.00,"#FF0000"),(.11,"#FF454A"),(.27,"#E82033"),
        (.45,"#A63151"),(.62,"#48628C"),(.77,"#0079BC"),
        (1.00,"#94AFE0")],N=4096)
    yc=.69
    head_base=X0+.98
    clip=Polygon([(X0,yc),(head_base,yc+.37),(head_base,yc+.18),
                  (X1,yc+.18),(X1,yc-.18),(head_base,yc-.18),
                  (head_base,yc-.37)],closed=True,transform=ax.transData)
    # A single continuous image avoids PDF-viewer seams between narrow
    # vector strips. The outline, labels, and other artwork remain vector.
    gradient=np.tile(np.linspace(0,1,4096,dtype=np.float32),(64,1))
    band=ax.imshow(gradient,extent=(X0,X1,yc-.37,yc+.37),
                   cmap=cmap,vmin=0,vmax=1,origin="lower",
                   interpolation="bilinear",resample=True,zorder=1)
    band.set_clip_path(clip)
    text(ax,X0+.68,yc,"Compressive",fs=11.5,ha="left",color="black")
    # Optical centering: this label has no descenders, unlike Compressive.
    text(ax,X1-.10,yc-.035,"Strain-free",fs=11.5,ha="right",color="black")


def panel_b(ax):
    b=1.1
    text(ax,.73,high(X1,b)+PANEL_LABEL_GAP_B,"(b)",fs=13.5,ha="left",fontweight="normal")
    edges=[X0,X0+.16*L,X0+.84*L,X1]
    left_center=(edges[0]+edges[1])/2
    right_center=(edges[2]+edges[3])/2
    text(ax,left_center,high(left_center,b)+.28,"Accumulation",fs=11)
    text(ax,right_center,high(right_center,b)+.28,"Depletion",fs=11)
    for left,right,col in zip(edges[:-1],edges[1:],["#F1BB8D","#D4EEEE","#94AED9"]):
        ax.add_patch(Polygon(region_vertices(left,right,b),closed=True,
                             facecolor=col,edgecolor="none",zorder=1))
    ax.add_patch(Polygon(region_vertices(X0,X1,b),closed=True,
                         facecolor="none",edgecolor=INK,lw=1.25,zorder=3))

    # Accumulation: a balanced 2-1-2-1-2 stagger, rather than two rails.
    left_positions=[(.037,.88),(.124,.86),(.080,.69),
                    (.040,.49),(.122,.47),(.077,.30),
                    (.037,.12),(.122,.13)]
    for q,f in left_positions:
        x=X0+q*L
        vacancy(ax,x,low(x,b)+(high(x,b)-low(x,b))*f)
    # Interior: a nearly uniform staggered array with the same symbol size
    # as both ends. Counts are illustrative, not calibrated concentrations.
    # Mechanism annotations occupy row gaps; markers are not displaced to
    # surround the annotations. End-region widths and arrangements are fixed.
    row_a=np.linspace(3.30,8.20,4)
    row_b=row_a+.20
    # Respace the interior rows and mechanisms about the unchanged midline,
    # while keeping all marker sizes and the four-column arrangement.
    midline=b+H/2
    old_mid_height=H*(1-STRAIN_VISUAL/2)
    content_scale=(old_mid_height+REFERENCE_PADDING)/old_mid_height
    def interior_y(y):
        return midline+(y-midline)*content_scale
    middle_positions=[(x,high(x,b)-.20) for x in row_a]
    middle_positions += [(x,interior_y(3.36)) for x in row_b]
    middle_positions += [(x,interior_y(2.40)) for x in row_a]
    middle_positions += [(x,low(x,b)+.20) for x in row_b]
    for x,y in middle_positions:
        vacancy(ax,x,y)
    # Depletion: five staggered symbols with an even vertical rhythm.
    for q,f in [(.882,.86),(.952,.68),(.887,.50),(.955,.32),(.884,.14)]:
        x=X0+q*L
        vacancy(ax,x,low(x,b)+(high(x,b)-low(x,b))*f)

    # Preserve the original concise mechanism labels and central arrow layout.
    center=(edges[1]+edges[2])/2
    left,right=center-1.82,center+1.82
    for y,name,color,direction in [
        (3.65,"Strain-driven drift",STRAIN,-1),
        (2.72,"Concentration-driven diffusion",DIFFUSION,1),
        (1.79,"Electrostatic feedback",ELECTRIC,1),
    ]:
        y=interior_y(y)
        text(ax,center,y+.25,name,fs=10.5)
        a,z=(right,left) if direction<0 else (left,right)
        ax.arrow(a,y,z-a,0,width=.15,head_width=.33,head_length=.34,
                 length_includes_head=True,facecolor=color,edgecolor="none",zorder=5)

    brace(ax,edges[0],edges[1],low(edges[0],b)-.035,
          end_y=low(edges[1],b)-.035)
    brace(ax,edges[2],edges[3],low(edges[2],b)-.035,
          end_y=low(edges[3],b)-.035)
    for center,label in [(left_center,"Positive space charge"),
                         (right_center,"Negative space charge")]:
        text(ax,center,low(center,b)-.035-.24-.25,label,fs=10.3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    fig_width,fig_height=6.10,5.66
    fig=plt.figure(figsize=(fig_width,fig_height),facecolor="white")
    # Keep the previous drawing scale, reserve space above (a), and reduce
    # the whitespace between the gradient arrow and the second panel.
    unit=.48*6.00/6.15
    left_margin=(fig_width-12*unit)/2
    ax_a=fig.add_axes([left_margin/fig_width,2.616/fig_height,
                      12*unit/fig_width,6.45*unit/fig_height])
    ax_b=fig.add_axes([left_margin/fig_width,.06/fig_height,
                      12*unit/fig_width,5.70*unit/fig_height])
    for ax,ymax in [(ax_a,6.45),(ax_b,5.70)]:
        ax.set(xlim=(0,12),ylim=(0,ymax))
        ax.set_aspect("equal")
        ax.axis("off")
    panel_a(ax_a)
    panel_b(ax_b)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    outside=[]
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text():
            box=t.get_window_extent(renderer)
            if box.x0<fig.bbox.x0 or box.x1>fig.bbox.x1 or box.y0<fig.bbox.y0 or box.y1>fig.bbox.y1:
                outside.append(t.get_text())
    if outside:
        raise RuntimeError(outside)
    fig.savefig(out/"figure1.pdf",dpi=600,metadata={"Title":"Figure 1 - consistent specimen geometry","Subject":"Qualitative schematic; no numerical concentration data"})
    fig.savefig(out/"figure1.png",dpi=600)
    plt.close(fig)
    print("Created:",out/"figure1.pdf")


if __name__=="__main__":
    main()
