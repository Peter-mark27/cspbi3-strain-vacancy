# Figure 2: final computational workflow.
# Requires Matplotlib and a LaTeX installation with amsmath and bm.
# Outputs figure2.pdf and figure2.png to reproduced_figures/ by default.
# Use --output-dir to choose another directory. Serif fonts have a fallback.
from pathlib import Path
import argparse
import os, tempfile
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'figure2_matplotlib_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description="Reproduce the computational-workflow schematic.")
parser.add_argument("--output-dir", type=Path, default=ROOT / "reproduced_figures")
args = parser.parse_args()
OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)
FONT_NAME = font_manager.FontProperties(
    family=["Times New Roman", "Liberation Serif", "DejaVu Serif"]
).get_name()
plt.rcParams.update({'font.family': FONT_NAME, 'mathtext.fontset': 'stix',
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(14.6, 13.15))
fig.subplots_adjust(0, 0, 1, 1)
ax.set(xlim=(0,1460), ylim=(1315,0)); ax.axis('off')
INK='#191D20'; PEACH='#F7E1D5'; GREEN='#B6D7CE'; BLUE='#B6D0E7'
RED='#B9463B'; ARROW='#396D94'

def text(x,y,s,size=17,weight='normal',color=INK):
    return ax.text(x,y,s,ha='center',va='center',fontsize=size,fontweight=weight,
                   color=color,linespacing=1.26)

def box(x,y,w,h,s='',fill=PEACH,size=17):
    ax.add_patch(Rectangle((x,y),w,h,facecolor=fill,edgecolor=INK,lw=1.6,zorder=2))
    if s: text(x+w/2,y+h/2,s,size)

def panel(x,y,w,h,title,fill):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=23',
        facecolor=fill,edgecolor='#99A4A2',lw=1.5,linestyle=(0,(7,4)),zorder=0))
    text(x+w/2,y+34,title,20,'bold')

def arrow(x1,y1,x2,y2,color=INK,lw=1.55):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',
                 mutation_scale=14,lw=lw,color=color,shrinkA=0,shrinkB=0,zorder=3))

def poly(points,color=ARROW,arrowend=True):
    for a,b in zip(points[:-2],points[1:-1]):
        ax.plot([a[0],b[0]],[a[1],b[1]],color=color,lw=1.5,zorder=1)
    a,b=points[-2:]
    if arrowend: arrow(*a,*b,color)
    else: ax.plot([a[0],b[0]],[a[1],b[1]],color=color,lw=1.5,zorder=1)



plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'
plain_text = text
def text(x,y,s,size=17,weight='normal',color=INK):
    if s.startswith('$') and s.endswith('$'):
        return ax.text(x,y,r'$\displaystyle '+s[1:],ha='center',va='center',
            fontsize=20,fontfamily='serif',fontweight='normal',color=color,
            usetex=True)
    return plain_text(x,y,s,size,weight,color)



def step(y,h,title,detail=None,text_offset=0,extra_gap=0):
    box(42,y,330,h,fill=PEACH)
    if detail:
        text(207,y+29+text_offset-extra_gap/2,title,18,'bold')
        text(207,y+h-31+text_offset+extra_gap/2,detail,17)
    else: text(207,y+h/2+text_offset,title,18,'bold')

step(65,110,'Model setup','Geometry, material properties\nand loading',text_offset=-5.5)
step(261.666667,90,'Atomistic\nparameterization')
step(438.333333,95,'Mechanical equilibrium','Plane stress; strain field')
step(620,95,'Local transport parameters',r'$E_f,\ E_m,\ D$',extra_gap=3)
step(801.666667,95,'Coupled vacancy transport\nand electrostatics')
step(983.333333,95,'Convergence verification','Mesh and time step')
step(1155,110,'Results and analysis','Evolution, steady states\nand parameter effects',text_offset=-5.5)
for a,b in [(175,261.666667),(351.666667,438.333333),(533.333333,620),(715,801.666667),(896.666667,983.333333),(1078.333333,1155)]:
    arrow(207,a,207,b)

panel(465,65,950,340,'Atomistic parameterization','#EFF6F3')
xs=[492,802,1112]
for x in xs: box(x,140,275,130,fill=GREEN)
text(629.5,174,'Formation energy',18,'bold')
text(629.5,226,'Strain-gradient\ncalculations',17)
text(939.5,174,'Migration barrier',18,'bold')
text(939.5,226,'CI-NEB\ncalculations',17)
text(1249.5,174,'Vacancy diffusion',18,'bold')
text(1249.5,226,'MD\nArrhenius analysis',17)
for x,s in zip(xs,[r'$E_f^0,\ A_f$',r'$E_m^0,\ A_m$',r'$D_0$']):
    arrow(x+137.5,270,x+137.5,305,'#3D786B')
    box(x,305,275,62,s,GREEN,23)
arrow(372,306.666667,465,306.666667,RED,2)

# Explicit two-source transfer: atomistic coefficients and the mechanical field.
panel(465,485,950,200,'Parameter transfer','#FAF4EC')
text(940,554,'Atomistic parameters + computed strain field',17)
text(629.5,616,r'$E_f=E_f^0+A_f\varepsilon_{yy}$',21)
text(939.5,616,r'$E_m=E_m^0+A_m\varepsilon_{yy}$',21)
text(1249.5,616,r'$D=D_0\exp[-E_m/(RT)]$',19.5)
arrow(940,405,940,485,ARROW)

poly([(372,485.833333),(416,485.833333),(416,555),(465,555)],color=ARROW)

# The main-list entry and expanded transfer region are linked without adding
# a second physical dependency arrow.
ax.plot([372,465],[667.5,667.5],color='#929B9A',lw=1.3,linestyle=(0,(3,3)))

panel(465,755,950,510,'Coupled transport and electrostatics','#EDF3F9')
arrow(940,685,940,755,ARROW)

box(495,815,890,115,fill=BLUE)
text(940,841,'Initial and boundary conditions',18,'bold')
text(940,873,'Uniform, charge-neutral initial state',16)
text(940,904,'Blocking and electrically insulating boundaries',16)
box(495,968,470,190,fill=BLUE)
text(730,995,'Vacancy transport',18,'bold')
text(730,1045,r'$J=-D\left[\frac{\partial c}{\partial x}+\frac{c}{RT}\frac{\partial E_f}{\partial x}'
     r'+\frac{z_VFc}{RT}\frac{\partial\phi}{\partial x}\right]$',20)
text(730,1115,r'$\frac{\partial c}{\partial t}=-\frac{\partial J}{\partial x}$',23)
box(1085,968,300,190,fill=BLUE)
text(1235,995,'Electrostatics',18,'bold')
text(1235,1045,r'$-\nabla\!\cdot(\epsilon\nabla\phi)=\rho$',23)
text(1235,1115,r'$\rho=F(z_Vc-c_b)$',23)
arrow(730,930,730,968,ARROW)
arrow(1235,930,1235,968,ARROW)
arrow(965,1045,1085,1045,ARROW)
text(1025,1024,'Space charge',14,color=ARROW)
arrow(1085,1097,965,1097,ARROW)
text(1025,1127,'Electrostatic\ndrift',14,color=ARROW)
arrow(372,849.166667,465,849.166667,RED,2)
box(650,1200,585,43,fill=BLUE)
# Small optical adjustment for the predominantly upper-weighted glyphs.
text(942.5,1223,'Fully coupled time-dependent solution',18,'bold')
poly([(730,1158),(730,1179),(852,1179),(852,1200)])
poly([(1235,1158),(1235,1179),(1030,1179),(1030,1200)])


for artist in ax.texts:
    label=artist.get_text()
    x,y=artist.get_position()
    if y == 616 and artist.get_usetex():
        artist.set_fontsize(22)
        artist.set_y(623)
        if x > 1200: artist.set_x(1235)
fig.canvas.draw()
for artist in ax.texts:
    if artist.get_position()[0] > 465 and 485 < artist.get_position()[1] < 685:
        bounds=artist.get_window_extent(fig.canvas.get_renderer()).transformed(ax.transData.inverted())
        assert bounds.x0 > 480 and bounds.x1 < 1400, (artist.get_text(),bounds)


for patch in ax.patches:
    if isinstance(patch, Rectangle) and abs(patch.get_x()-42)<0.01:
        patch.set_linewidth(1.3)
for artist in ax.texts:
    label=artist.get_text()
    if label == 'Atomistic parameters + computed strain field':
        artist.set_y(565)
    elif label == 'Uniform, charge-neutral initial state':
        artist.set_text('Uniform, charge-neutral state')
    elif label == 'Blocking and electrically insulating boundaries':
        artist.set_text('Blocking and insulating boundaries')

# Balance the transfer-panel content without changing its frame or type sizes.
transfer_formulas = []
for artist in ax.texts:
    x, y = artist.get_position()
    if artist.get_text() == 'Atomistic parameters + computed strain field':
        artist.set_y(573)
    elif artist.get_usetex() and x > 465 and 600 < y < 650:
        artist.set_y(634)
        transfer_formulas.append(artist)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
transfer_formulas.sort(key=lambda artist: artist.get_position()[0])
widths = [artist.get_window_extent(renderer).transformed(
    ax.transData.inverted()).width for artist in transfer_formulas]
assert len(widths) == 3
left, right = 505, 1375
gap = (right - left - sum(widths)) / 2
assert gap > 0
cursor = left
for artist, width in zip(transfer_formulas, widths):
    artist.set_x(cursor + width / 2)
    cursor += width + gap

for ext in ['pdf','png']:
    fig.savefig(OUT/f'figure2.{ext}',dpi=600,facecolor='white')
plt.close(fig)
print(OUT/'figure2.pdf')
