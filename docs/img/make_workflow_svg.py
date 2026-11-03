"""Generate docs/img/workflow.svg, the MetaCheese workflow figure.

Redraw of page 1 of the project diagram, with the step/artifact IDs used in
the README and in the scripts. Edit this file and run
    python3 docs/img/make_workflow_svg.py docs/img/workflow.svg
"""
import sys
from xml.sax.saxutils import escape

W, H = 1100, 1510
out = []
add = out.append

STEP_FILL, STEP_TEXT = "#9fd0ec", "#0b3954"
ART_FILL, ART_STROKE = "#ffd43b", "#b08900"
STEP_ID, ART_ID = "#1c8fd6", "#e8962e"
PANEL = "#e3eef8"
PLANNED = "#ffffff"


def text(x, y, s, size=13, weight="normal", anchor="middle", fill="#1d2b36", family="Helvetica, Arial, sans-serif", style=""):
    add(f'<text x="{x}" y="{y}" font-family="{family}" font-size="{size}" font-weight="{weight}" '
        f'text-anchor="{anchor}" fill="{fill}" {style}>{escape(s)}</text>')


def badge(x, y, label, color):
    add(f'<rect x="{x-17}" y="{y-9}" width="34" height="18" rx="9" fill="{color}"/>')
    text(x, y + 4, label, size=10, weight="bold", fill="#fff")


def step(x, y, w, h, sid, lines, script=None, planned=False):
    dash = ' stroke-dasharray="6 4"' if planned else ""
    fill = PLANNED if planned else STEP_FILL
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{fill}" stroke="#4a90b8" stroke-width="1.2"{dash}/>')
    badge(x + 4, y + 2, sid, STEP_ID)
    n = len(lines) + (1 if script else 0)
    y0 = y + h / 2 - (n - 1) * 8 + 4
    for i, ln in enumerate(lines):
        text(x + w / 2, y0 + i * 16, ln, size=13, weight="bold" if i == 0 else "normal", fill=STEP_TEXT)
    if script:
        text(x + w / 2, y0 + len(lines) * 16, script, size=11, fill="#35617d",
             family="Menlo, Consolas, monospace")


def artifact(x, y, w, h, aid, lines, label=None, planned=False):
    # stacked document shape: two offset copies behind, wavy bottom edge
    for k in (2, 1, 0):
        dx = k * 4
        p = (f"M{x+dx},{y-dx} h{w} v{h-8} q{-w/4},{-8} {-w/2},0 t{-w/2},0 z")
        fill = PLANNED if planned else ART_FILL
        dash = ' stroke-dasharray="5 3"' if planned else ""
        add(f'<path d="{p}" fill="{fill}" stroke="{ART_STROKE}" stroke-width="1"{dash}/>')
    for i, ln in enumerate(lines):
        text(x + w / 2, y + h / 2 - (len(lines) - 1) * 7 + i * 14 - 2, ln, size=11, weight="bold", fill="#3d3000")
    badge(x + 2, y + h - 6, aid, ART_ID)
    if label:
        for i, ln in enumerate(label.split("\n")):
            text(x + w + 14, y + 16 + i * 15, ln, size=12, anchor="start", fill="#1d2b36")


def arrow(points, dashed=False):
    d = "M" + " L".join(f"{px},{py}" for px, py in points)
    dash = ' stroke-dasharray="5 4"' if dashed else ""
    add(f'<path d="{d}" fill="none" stroke="#1d2b36" stroke-width="1.4" marker-end="url(#ah)"{dash}/>')


REPORT = "#7b4fbf"  # lines into the integrated report (R.1)


def feed(points, end=False):
    """Purple line from an artifact the report reads to a collector line."""
    d = "M" + " L".join(f"{px},{py}" for px, py in points)
    mk = ' marker-end="url(#ahp)"' if end else ""
    add(f'<path d="{d}" fill="none" stroke="{REPORT}" stroke-width="1.3" stroke-dasharray="4 3"{mk}/>')
    if not end:
        px, py = points[-1]
        add(f'<circle cx="{px}" cy="{py}" r="2.6" fill="{REPORT}"/>')


def panel(x, y, w, h, title, anchor="start"):
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="14" fill="{PANEL}"/>')
    tx = x + 14 if anchor == "start" else x + w - 14
    text(tx, y + 20, title, size=14, weight="bold", anchor=anchor, fill="#1f6aa5")


add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
add('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
    '<path d="M0,0 L10,5 L0,10 z" fill="#1d2b36"/></marker>'
    '<marker id="ahp" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
    f'<path d="M0,0 L10,5 L0,10 z" fill="{REPORT}"/></marker></defs>')
add(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')

# ---------------------------------------------------------------- legend
add('<rect x="20" y="20" width="170" height="140" rx="10" fill="#f4f7fa" stroke="#d5dee6"/>')
text(32, 40, "Legend", size=13, weight="bold", anchor="start", fill="#1f6aa5")
add(f'<rect x="32" y="47" width="26" height="13" rx="4" fill="{STEP_FILL}"/>'); text(66, 58, "computational step", size=11, anchor="start")
add(f'<rect x="32" y="65" width="26" height="13" rx="4" fill="{PLANNED}" stroke="#4a90b8" stroke-dasharray="4 3"/>'); text(66, 76, "planned, not in code", size=11, anchor="start")
add(f'<rect x="32" y="83" width="26" height="13" fill="{ART_FILL}" stroke="{ART_STROKE}"/>'); text(66, 94, "artifact (files)", size=11, anchor="start")
add(f'<path d="M32,107 h26" stroke="{REPORT}" stroke-width="1.5" stroke-dasharray="4 3"/>'); text(66, 111, "read by report R.1", size=11, anchor="start")
badge(45, 128, "X.n", STEP_ID); badge(85, 128, "X.n", ART_ID); text(108, 132, "step / artifact ID", size=11, anchor="start")
text(32, 148, "monospace = script", size=10, anchor="start", fill="#35617d", family="Menlo, Consolas, monospace")

# ---------------------------------------------------------------- top row: QC and host index
step(250, 40, 170, 62, "Q.1", ["Quality check", "(FastQC, MultiQC)"], "00")
artifact(455, 46, 80, 50, "Q.1", ["html", "zip"], "QC reports\n+ MultiQC")
artifact(640, 46, 80, 50, "I.2", ["FASTA"], None)
text(680, 36, "Host genome", size=12, weight="bold")
step(760, 40, 170, 62, "I.2", ["Indexing", "(bowtie2-build)"], "build_bowtie2_index.sh")
artifact(970, 46, 80, 50, "I.3", ["bt2"], None)
text(1012, 36, "Bowtie2 index", size=12, weight="bold")
arrow([(420, 71), (452, 71)])
arrow([(722, 71), (757, 71)])
arrow([(930, 71), (967, 71)])

# ---------------------------------------------------------------- experimental data + preprocessing
# Rule: a step's real output (passed to the next step) and its logs/settings
# (read only by the report) are separate artifacts.
panel(20, 182, 170, 200, "Experimental data")
artifact(55, 228, 90, 56, "I.1", ["FASTQ", "*_1/_2.fq.gz"], None)
text(105, 300, "paired-end reads", size=11)
panel(205, 182, 875, 200, "Preprocessing")
step(250, 222, 190, 66, "P.1", ["Adapter removal", "(AdapterRemoval)"], "00-01")
artifact(470, 228, 80, 54, "P.1.1", ["FASTQ"], "clean\nreads")
artifact(300, 312, 90, 48, "P.1.2", ["settings"], "trimming\nstats")
step(620, 222, 190, 66, "P.2", ["Host DNA filtering", "(Bowtie2)"], "03")
artifact(860, 228, 80, 54, "P.2.1", ["FASTQ"], "microbial\nreads")
artifact(670, 312, 90, 48, "P.2.2", ["log"], "host\nfraction")
arrow([(150, 256), (247, 256)])
arrow([(105, 225), (105, 170), (230, 170), (230, 71), (247, 71)])
arrow([(440, 255), (467, 255)])
arrow([(345, 288), (345, 309)])
arrow([(564, 255), (617, 255)])
arrow([(1012, 98), (1012, 140), (715, 140), (715, 219)])
arrow([(810, 255), (857, 255)])
arrow([(715, 288), (715, 309)])
# microbial reads (P.2.1) go down to the bus line feeding C.1, C.3, A.1, A.3
add('<path d="M900,284 L900,410" fill="none" stroke="#1d2b36" stroke-width="1.4"/>')

# Everything below the preprocessing panel is drawn 90 px lower.
DY = 90
add(f'<g transform="translate(0,{DY})">')

# ---------------------------------------------------------------- community-level profiling
panel(20, 330, 440, 915, "Community-level profiling")
add('<ellipse cx="72" cy="380" rx="38" ry="9" fill="#ffd43b" stroke="#b08900"/>')
add('<path d="M34,380 v30 a38,9 0 0 0 76,0 v-30" fill="#ffd43b" stroke="#b08900"/>')
badge(40, 372, "I.4", ART_ID)
text(72, 438, "MetaPhlAn DB", size=11, weight="bold")
text(72, 452, "(vJun23)", size=11)
step(150, 366, 260, 62, "C.1", ["Taxonomic profiling", "(MetaPhlAn 4.1)"], "04")
arrow([(110, 395), (147, 395)])
artifact(118, 462, 112, 52, "C.1.2", ["read", "alignments"], None)
artifact(270, 462, 115, 52, "C.1.1", ["taxonomic", "profiles"], None)
arrow([(210, 428), (175, 459)])
arrow([(350, 428), (327, 459)])
step(150, 544, 260, 70, "C.2", ["Profile merging,", "abundance and diversity"], "04 + calculate_diversity.R")
arrow([(327, 516), (327, 541)])
artifact(110, 646, 115, 56, "C.2.1", ["merged", "abundance"], None)
artifact(270, 646, 115, 56, "C.2.2", ["alpha/beta", "diversity"], None)
arrow([(210, 614), (167, 643)])
arrow([(350, 614), (327, 643)])
step(150, 740, 260, 70, "C.3", ["Functional profiling", "(HUMAnN 3.9)"], "04b-1")
arrow([(167, 702), (167, 722), (200, 722), (200, 737)])
artifact(195, 840, 110, 56, "C.3", ["genes", "pathways"], "functional\nprofiles\n(per sample)")
arrow([(280, 810), (280, 837)])
step(150, 928, 260, 70, "C.4", ["Community-level aggregation", "(HUMAnN join, renorm, regroup)"], "04b-last")
arrow([(250, 896), (250, 925)])
artifact(205, 1030, 100, 56, "C.4", ["CPM, KO,", "EC, MetaCyc"], "merged\nfunctional\ntables")
arrow([(280, 998), (280, 1027)])
# one bus line carries the microbial reads to C.1, C.3, A.1 and A.3
add('<path d="M280,320 L1060,320" fill="none" stroke="#1d2b36" stroke-width="1.4"/>')
arrow([(280, 320), (280, 363)])
arrow([(430, 320), (430, 775), (413, 775)])

# ---------------------------------------------------------------- assembly-specific profiling
panel(475, 330, 605, 915, "")
text(640, 350, "Assembly-specific profiling", size=14, weight="bold", anchor="start", fill="#1f6aa5")
step(500, 372, 210, 62, "A.1", ["Contig assembly", "(metaSPAdes)"], "05")
artifact(780, 376, 85, 54, "A.1", ["FASTA"], "contigs")
arrow([(710, 403), (777, 403)])
arrow([(605, 320), (605, 369)])
step(500, 470, 210, 62, "A.2", ["Contig filtering", "(> 500 bp)"], "06")
arrow([(822, 432), (822, 455), (690, 455), (690, 467)])
artifact(560, 566, 90, 54, "A.2", ["FASTA"], "filtered\ncontigs")
arrow([(605, 532), (605, 563)])
step(820, 470, 210, 62, "A.8", ["Quality check", "(QUAST)"], "planned", planned=True)
arrow([(710, 501), (817, 501)], dashed=True)
artifact(1000, 560, 70, 50, "A.8", ["reports"], None, planned=True)
arrow([(925, 532), (925, 585), (997, 585)], dashed=True)
step(820, 640, 210, 70, "A.3", ["Read mapping", "(Bowtie2 + samtools)"], "07, 08-09")
arrow([(1060, 320), (1060, 625), (1030, 625), (1030, 637)])
arrow([(650, 612), (760, 612), (760, 675), (817, 675)])
artifact(860, 742, 100, 54, "A.3.1", ["sorted BAM", "+ index"], None)
artifact(990, 742, 70, 50, "A.3.2", ["mapping", "log"], None)
arrow([(910, 710), (910, 739)])
arrow([(1015, 710), (1015, 739)])
step(500, 690, 210, 62, "A.4", ["Coverage and binning", "(MetaBAT2)"], "10-11")
arrow([(580, 620), (580, 687)])
arrow([(860, 770), (750, 770), (750, 721), (713, 721)])
artifact(520, 785, 85, 54, "A.4.1", ["bins", "FASTA"], None)
artifact(640, 785, 80, 54, "A.4.2", ["contig", "depth"], None)
arrow([(562, 752), (562, 782)])
arrow([(680, 752), (680, 782)])
step(820, 860, 210, 62, "A.5", ["Quality check", "(CheckM2)"], "12-13")
arrow([(590, 841), (590, 858), (790, 858), (790, 891), (817, 891)])
artifact(875, 955, 85, 54, "A.5", ["tsv"], "completeness,\ncontamination")
arrow([(925, 922), (925, 952)])
step(500, 950, 210, 62, "A.6", ["MAG filtering and", "metadata preparation"], "14")
arrow([(560, 841), (560, 947)])
arrow([(872, 982), (713, 982)])
arrow([(605, 1012), (605, 1037)])
artifact(560, 1040, 90, 54, "A.6", ["FASTA", "metadata"], "kept MAGs")
step(500, 1125, 210, 62, "A.7", ["MAG annotation", "(TORMES)"], "15")
arrow([(605, 1094), (605, 1122)])
artifact(780, 1131, 90, 54, "A.7", ["GFF, html", "Kraken2 txt"], "annotated MAGs,\nannotation reports")
arrow([(710, 1156), (777, 1156)])

# ---------------------------------------------------------------- joint step: both branches
step(420, 1271, 260, 66, "R.1", ["Integrated report", "(both branches)"], "16")
artifact(745, 1277, 85, 54, "R.1", ["HTML"], "interactive\nreport")
arrow([(680, 1304), (742, 1304)])
# R.1 reads the files of almost every step (scripts/report/make_report.py).
# Each dot is an artifact read by the report; the dots join two collector
# lines (drawn below, outside this group): one in the gap between the two
# branches, one on the right edge.
LB, RB = 467, 1090
feed([(397, 486), (LB, 486)])                               # C.1.1 taxonomic profiles
feed([(397, 670), (LB, 670)])                               # C.2.2 diversity
feed([(130, 704), (130, 714), (LB, 714)])                   # C.2.1 merged abundance
feed([(395, 865), (LB, 865)])                               # C.3 per-sample functional profiles
feed([(380, 1050), (LB, 1050)])                             # C.4 merged functional tables
feed([(556, 585), (LB, 585)])                               # A.2 filtered contigs
feed([(516, 812), (LB, 812)])                               # A.4.1 bins
feed([(556, 1068), (LB, 1068)])                             # A.6 kept MAGs + metadata
feed([(925, 403), (RB, 403)])                               # A.1 contigs
feed([(1064, 767), (RB, 767)])                              # A.3.2 mapping log
feed([(724, 812), (RB, 812)])                               # A.4.2 contig depth
feed([(1058, 975), (RB, 975)])                              # A.5 CheckM2 table
feed([(990, 1155), (RB, 1155)])                             # A.7 TORMES

# ---------------------------------------------------------------- downstream
add('<rect x="20" y="1365" width="1060" height="40" rx="8" fill="#fbe3d6"/>')
text(550, 1390, "Downstream analyses", size=14, weight="bold", fill="#7a3b1d")
add("</g>")

# collector lines into R.1 (absolute coordinates) and the feeds above the group
feed([(LB, 372), (LB, 1268 + DY)], end=True)                       # left collector
feed([(RB, 10), (RB, 1257 + DY), (650, 1257 + DY), (650, 1268 + DY)], end=True)  # right collector
feed([(500, 36), (500, 10), (RB, 10)])                      # Q.1 QC reports
feed([(360, 364), (360, 372), (LB, 372)])                   # P.1.2 trimming settings
feed([(842, 336), (RB, 336)])                               # P.2.2 host-removal log

add("</svg>")
open(sys.argv[1], "w").write("\n".join(out))
