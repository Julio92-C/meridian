#!/usr/bin/env python3
"""
Figure 1 — MERIDIAN: workflow schematic
High-resolution (300 DPI) publication figure.

Layout (top -> bottom, diamond / hourglass):
  Input row (4 data sources + config.yaml)
       v
  Core pipeline (R/00-R/02, centred)
       /                       \\
Community (R/04-R/08)   Functional (R/03 TPM + R/09-R/11 + shared a/b/DAA)
       \\                       /
     R/12 Network (centred)
       /                       \\
Quarto Dashboard          R/14 Panels
       v                       v
  R/13 Manifest    <-     Publication outputs
  (also fed by R/14 + Quarto)

Every module box is two-line (bold id + description) so long labels
never overflow the frame. A shared-methods box under the Functional
group records that alpha / beta-diversity / differential-abundance are
also computed on the functional gene tables (R/08 is two-level taxa+gene;
R/09-R/11 each carry per-domain alpha + PCoA + PERMANOVA). TPM
normalisation (R/03) is applied on the functional branch only; community
analyses run on the cleaned, non-TPM tables.

IMPORTANT: run with the real Python, not Inkscape's bundled one:
  & "C:\\Users\\julio\\AppData\\Local\\Programs\\Python\\Python313\\python.exe" generate_figure1.py
Out:  figure1_pipeline_schematic.png  (300 DPI)
      figure1_pipeline_schematic.svg  (vector / editable)
"""

import os, pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# --- Figure geometry ---------------------------------------------------------
FIG_W = 7.20   # inches  (183 mm - Bioinformatics double-column)
FIG_H = 7.70   # inches
DPI   = 300

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(-0.43, 7.45)
ax.axis('off')
fig.patch.set_facecolor('white')

# --- Colour palette  (face, edge, text) --------------------------------------
PAL = {
    'input': ('#cee5f5', '#2070a0', '#0e3255'),
    'cfg':   ('#fef3c7', '#c08000', '#5c3900'),
    'core':  ('#e8e8e8', '#505050', '#1a1a1a'),
    'comm':  ('#c8edcc', '#1e7d34', '#0b3515'),
    'func':  ('#f9d4d5', '#c0302b', '#620f0e'),
    'net':   ('#ddd5f5', '#5528a0', '#280a50'),
    'out':   ('#fef5aa', '#9a8200', '#3d3200'),
}
BG_COMM = '#e8f8ea'
BG_FUNC = '#fdf0f0'

# --- Box dimensions ----------------------------------------------------------
IBW = 1.20    # input-row box width
IBH = 0.44    # input-row box height (room for two-line sublabels)
CBW = 1.58    # analysis box width  (community / functional domains)
CBH = 0.34    # analysis box height
COREW = 1.74  # core box width
NBW = 1.96    # network box width
OBW2 = 1.85   # bottom-row output box width (two columns side by side)
OBH2 = 0.42   # bottom-row output box height (room for two-line subs)
METHW = 2.04  # functional shared-methods box width
BRAD = 0.065  # corner-rounding radius

# --- Helpers -----------------------------------------------------------------
def box(cx, cy, top, sub=None, ckey='core', top_fs=6.0, sub_fs=5.3,
        w=None, h=None, sub_italic=False):
    """Two-line (bold id + description) rounded-rect box centred at (cx, cy)."""
    w = w or CBW
    h = h or CBH
    fc, ec, tc = PAL[ckey]
    # shadow
    ax.add_patch(FancyBboxPatch(
        (cx-w/2+.012, cy-h/2-.012), w, h,
        boxstyle=f'round,pad=0,rounding_size={BRAD}',
        fc='#b0b0b0', ec='none', alpha=.25, zorder=2))
    # box
    ax.add_patch(FancyBboxPatch(
        (cx-w/2, cy-h/2), w, h,
        boxstyle=f'round,pad=0,rounding_size={BRAD}',
        lw=.8, ec=ec, fc=fc, zorder=3))
    if sub:
        ax.text(cx, cy+h*0.20, top, ha='center', va='center',
                fontsize=top_fs, fontweight='bold', color=tc, zorder=4)
        ax.text(cx, cy-h*0.23, sub, ha='center', va='center',
                fontsize=sub_fs, color=tc, zorder=4,
                style='italic' if sub_italic else 'normal',
                multialignment='center')
    else:
        ax.text(cx, cy, top, ha='center', va='center',
                fontsize=top_fs, fontweight='bold', color=tc, zorder=4,
                multialignment='center')

def arr(x1, y1, x2, y2, color='#444444', rad=0.0, lw=0.80, alpha=1.0,
        ls='-'):
    """Directional arrow from (x1,y1) to (x2,y2)."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle='->', color=color, lw=lw,
            mutation_scale=8, linestyle=ls,
            connectionstyle=f'arc3,rad={rad}',
            alpha=alpha),
        zorder=1)

def grp_rect(x, y, w, h, title, ckey):
    """Dashed group-background rectangle with title above."""
    fc, ec, tc = PAL[ckey]
    ax.add_patch(FancyBboxPatch((x, y), w, h,
        boxstyle='round,pad=0,rounding_size=.12',
        lw=0.85, ec=ec, fc=BG_COMM if ckey == 'comm' else BG_FUNC,
        linestyle='--', alpha=.65, zorder=1))
    ax.text(x+w/2, y+h+.030, title, ha='center', va='bottom',
            fontsize=6.8, color=ec, fontweight='bold', zorder=4)

# --- Layout constants --------------------------------------------------------
CX      = FIG_W / 2     # horizontal centre = 3.60
CX_COMM = 1.52          # community column centre
CX_FUNC = 5.64          # functional column centre

# Y positions (bottom -> top = increasing y in matplotlib)
Y_INP   = 7.10
Y_CFG   = 6.58
Y_R00   = 6.18
Y_R01   = 5.70
Y_R02   = 5.22

COMM_YS = [4.25, 3.79, 3.33, 2.87, 2.41]   # R/04-R/08
FUNC_YS = [4.25, 3.79, 3.33, 2.87]         # R/03 (TPM) + R/09-R/11
Y_FMETH = 2.41                             # shared-methods box (aligns with R/08)

Y_NET   = 1.42
Y_QRTO  = 0.80   # Quarto dashboard (left) + R/14 panels (right)
Y_MFST  = 0.18   # R/13 manifest (left) + publication outputs (right)

X_OUTL  = 2.50   # left output column: Quarto -> Manifest
X_OUTR  = 5.40   # right output column: R/14 panels -> Publication outputs

PAD = 0.11

# --- 1. GROUP BACKGROUNDS  (z = 1) -------------------------------------------
grp_rect(
    x=CX_COMM - CBW/2 - PAD,
    y=COMM_YS[-1] - CBH/2 - PAD,
    w=CBW + 2*PAD,
    h=COMM_YS[0] - COMM_YS[-1] + CBH + 2*PAD,
    title='Community Analyses', ckey='comm')

grp_rect(
    x=CX_FUNC - METHW/2 - PAD,
    y=Y_FMETH - CBH/2 - PAD,
    w=METHW + 2*PAD,
    h=FUNC_YS[0] - Y_FMETH + CBH + 2*PAD,
    title='Functional Profiling', ckey='func')

# --- 2. INPUT ROW ------------------------------------------------------------
inp_xs = [0.82, 2.08, 3.34, 4.60]
inp_data = [
    ('Kraken2 / Bracken',  'taxonomy tables'),
    ('ABRicate',           'CARD · VFDB ·\nPlasmidFinder'),
    ('Re-centrifuge',      'contamination tags'),
    ('Sample metadata',    'groups · covariates'),
]
for (lbl, sub), xi in zip(inp_data, inp_xs):
    box(xi, Y_INP, lbl, sub, 'input', top_fs=5.9, sub_fs=5.0,
        w=IBW, h=IBH, sub_italic=True)

# config.yaml - upper-right, staggered down so arrow is clear
X_CFG = 5.45
box(X_CFG, Y_CFG, 'config.yaml', 'one YAML per study', 'cfg',
    top_fs=5.9, sub_fs=5.0, w=1.15, h=IBH-0.06, sub_italic=True)

# --- 3. CORE PIPELINE (centred) ----------------------------------------------
core_items = [
    ('R/00', 'setup'),
    ('R/01', 'load inputs'),
    ('R/02', 'clean & decontaminate'),
]
for (mid, desc), y in zip(core_items, [Y_R00, Y_R01, Y_R02]):
    box(CX, y, mid, desc, 'core', top_fs=6.0, sub_fs=5.4, w=COREW)

# --- 4. COMMUNITY ANALYSES (left column) -------------------------------------
comm_items = [
    ('R/04', 'taxonomy'),
    ('R/05', 'relative abundance'),
    ('R/06', 'alpha diversity'),
    ('R/07', 'beta diversity (PCoA)'),
    ('R/08', 'differential abundance (DAA)'),
]
for (mid, desc), y in zip(comm_items, COMM_YS):
    box(CX_COMM, y, mid, desc, 'comm', top_fs=5.9, sub_fs=5.3)

# --- 5. FUNCTIONAL PROFILING (right column) ----------------------------------
# R/03 TPM normalisation lives at the TOP of this block: it runs on the
# functional gene tables only (community analyses use the cleaned R/02 tables).
func_items = [
    ('R/03', 'normalise (TPM)',          'core'),
    ('R/09', 'resistome · CARD',         'func'),
    ('R/10', 'virulome · VFDB',          'func'),
    ('R/11', 'mobilome · PlasmidFinder', 'func'),
]
for (mid, desc, ck), y in zip(func_items, FUNC_YS):
    box(CX_FUNC, y, mid, desc, ck, top_fs=5.9, sub_fs=5.3)

# Shared-methods box: alpha / beta / DAA are also run on the functional genes.
fc, ec, tc = PAL['func']
ax.add_patch(FancyBboxPatch(
    (CX_FUNC-METHW/2, Y_FMETH-CBH/2), METHW, CBH,
    boxstyle=f'round,pad=0,rounding_size={BRAD}',
    lw=0.85, ec=ec, fc='#fce7e7', linestyle='--', zorder=3))
ax.text(CX_FUNC, Y_FMETH+CBH*0.20,
        'α · β-diversity · DAA',
        ha='center', va='center', fontsize=6.0, fontweight='bold',
        color=tc, zorder=4)
ax.text(CX_FUNC, Y_FMETH-CBH*0.25,
        'also computed on functional gene tables',
        ha='center', va='center', fontsize=4.8, style='italic',
        color=tc, zorder=4)

# --- 6. NETWORK (centred) ----------------------------------------------------
box(CX, Y_NET, 'R/12  network', 'tripartite · chord · Sankey',
    'net', top_fs=6.0, sub_fs=5.4, w=NBW)

# --- 7. OUTPUTS (two columns: dashboard+manifest / panels+publication) ------
box(X_OUTL, Y_QRTO, 'Quarto dashboard', 'HTML · per-domain panels',
    'out', top_fs=5.8, sub_fs=5.0, w=OBW2, h=OBH2)
box(X_OUTR, Y_QRTO, 'R/14  panels', 'composite figures +\nsupp_tables.xlsx',
    'core', top_fs=5.8, sub_fs=5.0, w=OBW2, h=OBH2)
box(X_OUTL, Y_MFST, 'R/13  manifest', 'manifest.json (v1.2) ·\nkind-tagged catalogue',
    'out', top_fs=5.8, sub_fs=5.0, w=OBW2, h=OBH2)
box(X_OUTR, Y_MFST, 'Publication outputs', 'PNG / TIFF panels +\ntables workbook',
    'out', top_fs=5.8, sub_fs=5.0, w=OBW2, h=OBH2)

# --- 8. ARROWS ---------------------------------------------------------------
# 8a. each data input -> R/00 (fan inward)
for xi in inp_xs:
    ax.annotate('', xy=(CX, Y_R00 + CBH/2), xytext=(xi, Y_INP - IBH/2),
        arrowprops=dict(arrowstyle='->', color=PAL['input'][1], lw=.70,
                        mutation_scale=7, connectionstyle='arc3,rad=0'),
        zorder=1)

# 8b. config.yaml -> R/00 (curved, from the right)
arr(X_CFG - 1.15/2, Y_CFG, CX + COREW/2, Y_R00,
    color=PAL['cfg'][1], rad=0.18, lw=0.75)

# 8c. core chain R/00 -> R/01 -> R/02
for y_top, y_bot in [(Y_R00, Y_R01), (Y_R01, Y_R02)]:
    arr(CX, y_top - CBH/2, CX, y_bot + CBH/2, color='#404040', lw=0.80)

# 8d. R/02 -> Community block (short; lands on the block's near top corner)
arr(CX - COREW/2, Y_R02 - CBH/2, CX_COMM + CBW/2 - 0.10, COMM_YS[0] + CBH/2,
    color=PAL['comm'][1], rad=0.12, lw=0.85)

# 8e. R/02 -> Functional block (short; R/03 TPM is the first step)
arr(CX + COREW/2, Y_R02 - CBH/2, CX_FUNC - CBW/2 + 0.10, FUNC_YS[0] + CBH/2,
    color=PAL['func'][1], rad=-0.12, lw=0.85)

# 8f. functional domains -> shared-methods box (faint, "feeds")
arr(CX_FUNC, FUNC_YS[-1] - CBH/2, CX_FUNC, Y_FMETH + CBH/2,
    color=PAL['func'][1], lw=0.65, alpha=0.6, ls=':')

# 8g. Community bottom -> R/12 Network (converge from left)
arr(CX_COMM, COMM_YS[-1] - CBH/2,
    CX - NBW/2 + 0.22, Y_NET + CBH/2,
    color=PAL['net'][1], rad=-0.18, lw=0.85)

# 8h. Functional methods -> R/12 Network (converge from right)
arr(CX_FUNC, Y_FMETH - CBH/2,
    CX + NBW/2 - 0.22, Y_NET + CBH/2,
    color=PAL['net'][1], rad=0.18, lw=0.85)

# 8i. R/12 -> Quarto (left) and R/12 -> R/14 panels (right): from R/12's
# own left/right bottom corners to the centre-top of each child, mirroring
# the R/02 -> Community/Functional corner fan-out (8d/8e).
arr(CX - NBW/2, Y_NET - CBH/2, X_OUTL, Y_QRTO + OBH2/2,
    color=PAL['out'][1], rad=-0.15, lw=0.85)
arr(CX + NBW/2, Y_NET - CBH/2, X_OUTR, Y_QRTO + OBH2/2,
    color=PAL['core'][1], rad=-0.15, lw=0.85)

# 8j. Quarto -> Manifest (same column) and R/14 -> Publication outputs
arr(X_OUTL, Y_QRTO - OBH2/2, X_OUTL, Y_MFST + OBH2/2, color=PAL['out'][1], lw=0.85)
arr(X_OUTR, Y_QRTO - OBH2/2, X_OUTR, Y_MFST + OBH2/2, color=PAL['core'][1], lw=0.85)

# 8k. R/14 panels -> Manifest (cross-feed; manifest also catalogues panels)
arr(X_OUTR - OBW2/2 + 0.08, Y_QRTO - OBH2/2, X_OUTL + OBW2/2 - 0.08, Y_MFST + OBH2/2,
    color=PAL['core'][1], rad=0.16, lw=0.65, alpha=0.7, ls=':')

# --- 9. VERTICAL FLOW LABELS (left margin) -----------------------------------
for y, lbl, col in [
    ((Y_INP + Y_R00)/2,  'Upstream\nInputs',  PAL['input'][1]),
    ((Y_R00 + Y_R02)/2,  'Core\nPipeline',     '#404040'),
    (3.33,               'Domain\nAnalyses',   '#6a3b6a'),
    ((Y_NET + Y_MFST)/2, 'Report\nOutputs',    PAL['out'][1]),
]:
    ax.text(0.12, y, lbl, ha='center', va='center', rotation=90,
            fontsize=5.8, fontweight='bold', color=col, alpha=0.78,
            clip_on=False)

# --- 10. COLOUR LEGEND -------------------------------------------------------
leg_items = [
    ('Input data',        'input'),
    ('Configuration',     'cfg'),
    ('Processing step',   'core'),
    ('Community module',  'comm'),
    ('Functional module', 'func'),
    ('Network module',    'net'),
    ('Report / output',   'out'),
]
SWATCH = 0.105               # smaller colour squares
LEG_FS = 4.8
GAP_ST = 0.05                # swatch -> text gap
GAP_IT = 0.22                # gap between legend items
ly = Y_MFST - CBH/2 - 0.34
ax.axhline(y=Y_MFST - CBH/2 - 0.16, xmin=0.02, xmax=0.98,
           color='#cccccc', lw=0.5, zorder=0)

# Measure the REAL rendered text width of each label (data units) so a label
# can never overlap the next swatch, then centre the packed row.
fig.canvas.draw()
_rend = fig.canvas.get_renderer()
_inv = ax.transData.inverted()
texts, twidths = [], []
for lbl, ckey in leg_items:
    t = ax.text(0, ly, lbl, ha='left', va='center', fontsize=LEG_FS,
                color='#333333', zorder=4, clip_on=False)
    bb = t.get_window_extent(renderer=_rend)
    p0 = _inv.transform((bb.x0, bb.y0))
    p1 = _inv.transform((bb.x1, bb.y0))
    texts.append(t)
    twidths.append(p1[0] - p0[0])

widths = [SWATCH + GAP_ST + tw for tw in twidths]
lx = (FIG_W - (sum(widths) + GAP_IT * (len(leg_items) - 1))) / 2
for (lbl, ckey), t, witem in zip(leg_items, texts, widths):
    fc, ec, tc = PAL[ckey]
    ax.add_patch(FancyBboxPatch(
        (lx, ly - SWATCH/2), SWATCH, SWATCH,
        boxstyle='round,pad=0,rounding_size=.025',
        lw=.55, fc=fc, ec=ec, zorder=3, clip_on=False))
    t.set_position((lx + SWATCH + GAP_ST, ly))
    lx += witem + GAP_IT

# --- 11. SAVE ----------------------------------------------------------------
out_dir = pathlib.Path(__file__).parent
png_out = out_dir / 'figure1_pipeline_schematic.png'
svg_out = out_dir / 'figure1_pipeline_schematic.svg'

fig.savefig(png_out, dpi=DPI, bbox_inches='tight', facecolor='white')
fig.savefig(svg_out,          bbox_inches='tight', facecolor='white')
plt.close(fig)

sz = os.path.getsize(png_out)
print(f'Saved {png_out}  ({DPI} DPI,  {sz:,} bytes)')
print(f'Saved {svg_out}  (vector)')
