from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import _repro  # noqa: F401  -- deterministic PDF output; see figures/src/_repro.py

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, roc_curve

parser = argparse.ArgumentParser(
    description='Regenerate the COMAVI ISDS-v1 manuscript figures.'
)
parser.add_argument(
    '--analysis-dir',
    type=Path,
    required=True,
    help='Directory containing the regenerated ISDS-v1 analysis outputs.',
)
parser.add_argument(
    '--out-dir',
    type=Path,
    required=True,
    help='Directory for regenerated figure files.',
)
args = parser.parse_args()

AN = args.analysis_dir.expanduser().resolve()
OUT = args.out_dir.expanduser().resolve()

if not AN.is_dir():
    raise SystemExit(f'Analysis directory not found: {AN}')

OUT.mkdir(parents=True, exist_ok=True)

NAVY = '#18364A'
BLUE = '#3B82B8'
TEAL = '#2A9D8F'
ORANGE = '#E67E22'
PURPLE = '#7251B5'
RED = '#C44E52'
GREEN = '#5A9A66'
GOLD = '#C49A35'
LIGHT = '#F4F7F9'
MID = '#D7E0E6'
DARK = '#263238'
GRAY = '#6B7280'

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 10.5,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'xtick.labelsize': 9.5,
    'ytick.labelsize': 9.5,
    'figure.dpi': 170,
    'savefig.dpi': 300,
})


def panel_label(ax, label):
    ax.text(-0.10, 1.06, label, transform=ax.transAxes, fontsize=15, fontweight='bold', color=NAVY, va='top')


def save(fig, name):
    fig.savefig(OUT / name, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def box(ax, xy, wh, title, body, fc, ec, fontsize=9.6):
    x, y = xy; w, h = wh
    p = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.012,rounding_size=0.02', facecolor=fc, edgecolor=ec, linewidth=1.4)
    ax.add_patch(p)
    ax.text(x+w/2, y+h*0.67, title, ha='center', va='center', fontsize=10.2, fontweight='bold', color=DARK)
    ax.text(x+w/2, y+h*0.31, body, ha='center', va='center', fontsize=fontsize, color=DARK, linespacing=1.25)
    return p


def arrow(ax, x1, y1, x2, y2, color=GRAY):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=14,linewidth=1.3,color=color))


# Figure 1: unified workflow
# Layout follows the author's hand-drawn replacement: four labelled columns
# (Input, Structural Evidence, COMAVI Outputs, Practical Use) rather than the
# earlier three-stage flow, which never said what the two outputs are FOR. Drawn
# in code rather than pasted so it renders at the same resolution as the other
# figures and stays inside the embedded-figure gate.
fig, ax = plt.subplots(figsize=(13.4, 4.9))
ax.set_xlim(0, 1); ax.set_ylim(0.075, 1.0); ax.axis('off')
BUS_X = 0.478

# ALL FOUR evidence rows feed BOTH outputs. Verified against the tier-ablation
# record, not asserted: classify_mechanism_at() in apply_concordance_v5.py reads
# comavi_tier (it selects contact-driven vs burial-driven vocabulary), and
# blanking the tier changes 7 of 57 mechanism labels. An earlier version of this
# figure routed the tier to the priority score alone; that was wrong, and the
# assertion below is what stops it recurring.
_abl = json.loads((Path(__file__).resolve().parents[2] / 'reference_outputs'
                   / 'COMAVI_audit_summary.json').read_text())['tier_ablation']
assert _abl['n_labels_changed'] > 0, (
    'tier-ablation record says the tier changes no mechanism label; if that is now '
    'true this figure must be rewired to route the tier to the priority score only')

COLHEAD_Y = 0.955
for cx, head in [(0.085, 'Input'), (0.335, 'Structural Evidence'),
                 (0.615, 'COMAVI Outputs'), (0.875, 'Practical Use')]:
    ax.text(cx, COLHEAD_Y, head, ha='center', va='center',
            fontsize=12.5, fontweight='bold', color=GRAY)

EV_X, EV_W, EV_H = 0.205, 0.245, 0.175
ev = [(0.735, 'Stability in monomer', '#DCE9F7', BLUE),
      (0.525, 'Stability in multimer', '#DFF0E6', GREEN),
      (0.315, 'Binding energy by partner', '#FDECD9', ORANGE),
      (0.105, 'Structural-context tier', '#E7E1F5', PURPLE)]

# Every box first, every connector second. Drawing an arrow before the box it
# points into lets the box's rounded padding cover the arrowhead, which is how
# an earlier version of this figure shipped with headless connectors.
box(ax, (0.005, 0.545), (0.16, 0.20), 'Missense variant', '', LIGHT, GRAY, 9.0)
box(ax, (0.005, 0.175), (0.16, 0.24), 'Isolated / assembled\nstructural models', '', LIGHT, GRAY, 9.0)
for y, t, fc, ec in ev:
    box(ax, (EV_X, y), (EV_W, EV_H), t, '', fc, ec, 9.0)
box(ax, (0.505, 0.545), (0.225, 0.245), 'Priority score (ISDS)',
    'Ranks strength of modeled\nstructural-disruption evidence', '#F7DEDE', '#C06C6C', 8.6)
box(ax, (0.505, 0.155), (0.225, 0.245), 'Mechanism profile',
    'Assigns per-variant\nmechanisms and context', '#D9EBEA', TEAL, 8.6)
box(ax, (0.775, 0.545), (0.220, 0.245), 'Prioritization',
    'Choose variants for\nfurther evaluation', '#F7DEDE', '#C06C6C', 8.6)
box(ax, (0.775, 0.155), (0.220, 0.245), 'Experimental design',
    'Test stability, assembly\nand/or interaction', '#D9EBEA', TEAL, 8.6)

def link(x1, y1, x2, y2, color=GRAY):
    ax.plot([x1, x2], [y1, y2], color=color, lw=1.1, solid_capstyle='round', zorder=5)

def tip(x1, y1, x2, y2, color=GRAY):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                                 mutation_scale=15, linewidth=1.3, color=color, zorder=6))

for y, _, _, _ in ev:
    link(0.168, y + EV_H / 2, EV_X - 0.004, y + EV_H / 2)
    link(EV_X + EV_W + 0.004, y + EV_H / 2, BUS_X, y + EV_H / 2)
BUS_TOP, BUS_BOT = ev[0][0] + EV_H / 2, ev[-1][0] + EV_H / 2
link(BUS_X, BUS_BOT, BUS_X, BUS_TOP)
tip(BUS_X, (BUS_TOP + BUS_BOT) / 2, 0.499, 0.668)
tip(BUS_X, (BUS_TOP + BUS_BOT) / 2, 0.499, 0.278)
tip(0.734, 0.668, 0.770, 0.668, '#C06C6C')
tip(0.734, 0.278, 0.770, 0.278, TEAL)
save(fig,'figure1_unified_comavi_workflow.png')

# Figure 2: population map
fig, ax = plt.subplots(figsize=(11.5, 5.6))
ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
ax.text(0.02,0.94,'The resource supports two complementary evaluation tasks',fontsize=16,fontweight='bold',color=NAVY,va='top')
box(ax,(0.04,0.35),(0.20,0.35),'61-variant resource','14 protein systems\n49 interaction variants\n12 BRCT variants',LIGHT,NAVY,9.1)
arrow(ax,0.24,0.52,0.32,0.68)
arrow(ax,0.24,0.52,0.32,0.34)
box(ax,(0.32,0.55),(0.29,0.27),'57 mechanism-gradeable','Whole-variant pattern score\nThree energetic axes\nInteraction + BRCT systems','#EAF3F8',BLUE,9.0)
box(ax,(0.32,0.20),(0.29,0.27),'4 ungraded resource cases','2 mechanisms missing\nrequired context\n2 mechanism-uncommitted\nreference cases','#F5F5F5',DARK,9.0)
arrow(ax,0.61,0.69,0.66,0.69)
_pop=json.load(open(AN/'ISDS_v1_summary.json'))['population']
box(ax,(0.66,0.55),(0.31,0.27),f"{_pop['n']} structural-prioritization",f"Tier-carrying interaction variants\n{_pop['positive']} modeled structural mechanisms\n{_pop['negative']} with no committed lesion",'#E9F6F2',TEAL,9.0)
ax.text(0.34,0.12,'Other analysis subsets are defined once and referenced consistently: 15 direct-energy comparators, 44 measured destabilizers, and 47 AlphaMissense complete cases.',fontsize=9.3,color=GRAY)
save(fig,'figure2_population_map.png')

# Figure 3: mechanism localization performance
summary=json.load(open(AN/'ISDS_v1_summary.json'))
fig, axes=plt.subplots(1,2,figsize=(11.8,5.2),gridspec_kw={'width_ratios':[1.2,1]})
ax=axes[0]
labels=['Monomer fold','Complex context','Binding','All energetic axes',
        'Structural-context tier']
# Read the Figure 2 anchors from the analysis summary, which computes them from
# the canonical. They were duplicated literals here until v7.7, so a ledger
# correction left this panel plotting a stale partition.
#
# PRIMARY convention (the 57 mechanism-gradeable variants) — the same population
# panel b plots and the manuscript caption quotes. The summary also carries the
# all-row convention under 'all_row'; mixing the two is what desynchronised this
# panel from its caption (19/26 plotted against 19/25 quoted).
_mlv=summary['mechanism_localization_final_values']
assert _mlv['convention'].startswith('primary'), _mlv['convention']
def _nv(key):
    frac=_mlv[key].split(' = ')[0]
    a,b=frac.split('/')
    return frac, float(a)/float(b)
nums,vals=zip(*[_nv(k) for k in ('monomer_agreement','complex_context_agreement',
                                 'binding_agreement','all_energetic_agreement',
                                 'tier_agreement')])
nums,vals=list(nums),list(vals)
y=np.arange(len(labels))
bars=ax.barh(y,vals,color=[BLUE,TEAL,ORANGE,NAVY,PURPLE],height=.58)
ax.set_yticks(y,labels); ax.invert_yaxis(); ax.set_xlim(0,1.03); ax.set_xlabel('Direction-aware agreement')
ax.set_title('Physical axes remain separate',loc='left',fontweight='bold',color=NAVY)
ax.grid(axis='x',alpha=.2)
for bar,val,n in zip(bars,vals,nums): ax.text(val+.015,bar.get_y()+bar.get_height()/2,f'{n}  ({val:.3f})',va='center',fontsize=9.5,fontweight='bold')
panel_label(ax, 'A')
ax=axes[1]
_pnums,vals=zip(*[_nv(k) for k in ('whole_variant_score','interaction_subset_score',
                                   'brct_subset_score')])
_pnums,vals=list(_pnums),list(vals)
labels=['All gradeable variants','Interaction subset','BRCT fold subset']
colors=[NAVY,TEAL,PURPLE]
y=np.arange(3)
bars=ax.barh(y,vals,color=colors,height=.58)
ax.set_yticks(y,labels); ax.invert_yaxis(); ax.set_xlim(0,1.03); ax.set_xlabel('Whole-variant pattern score')
ax.set_title('Mechanism localization across systems',loc='left',fontweight='bold',color=NAVY)
ax.grid(axis='x',alpha=.2)
for bar,val,n in zip(bars,vals,_pnums): ax.text(val+.015,bar.get_y()+bar.get_height()/2,f'{n}  ({val:.3f})',va='center',fontsize=9.5,fontweight='bold')
# Cluster (whole-system) bootstrap CI, read from the stress-test output rather
# than pinned in this source. Test E resamples systems as clusters and is the
# honest interval; test D resamples variants and understates it, so the row is
# selected explicitly.
_REPO=Path(__file__).resolve().parents[2]
# Reads the COMMITTED stress outputs, not verification_output/. That directory
# is untracked -- git ls-files returns nothing for it -- so this figure was
# displaying a cluster CI from a file no one else has, it would abort on a
# clean clone, and the copy on this machine had been stale for two days and
# carried pre-correction values. reference_outputs/stress_tests/ is what ships.
_stress=_REPO/'reference_outputs'/'stress_tests'/'comavi_stress_tests.csv'
assert _stress.is_file(), f'{_stress} missing — run verification/stress_tests.py --out-dir reference_outputs/stress_tests first'
_sb=pd.read_csv(_stress)
_row=_sb[_sb['test'].eq('cluster bootstrap 95% CI (MC)')]
assert len(_row)==1, f'expected one cluster-bootstrap MC row, found {len(_row)}'
_lo,_hi=[float(x) for x in str(_row.iloc[0]['null_mean']).strip('[]').split(',')]
ax.text(0.02,-0.22,f'System-cluster 95% CI for the primary score: {_lo:.3f}-{_hi:.3f}',transform=ax.transAxes,fontsize=9.2,color=GRAY)
panel_label(ax, 'B')
fig.tight_layout(w_pad=3.0)
save(fig,'figure3_mechanism_localization.png')

# Figure 4: ISDS definition and performance
pv=pd.read_csv(AN/'ISDS_v1_per_variant.csv')
metrics=pd.read_csv(AN/'ISDS_v1_primary_metrics.csv')
topk=pd.read_csv(AN/'ISDS_v1_top_k.csv')
boot=pd.read_csv(AN/'ISDS_v1_system_cluster_bootstrap_summary.csv')
# PLOS Computational Biology requires one file per numbered figure. This block
# used to emit a single four-panel composite whose panels A and C were cited as
# S4 Fig while B and D were cited as Fig 6 -- one image serving two figures,
# which cannot ship. It now writes two two-panel files:
#
#   figure6_priority_recovery.png    Fig 6   score distribution + top-k recovery
#   figureS4_priority_definition.png S4 Fig  energy transformation + component ROC
#
# The former composite also carried an internal suptitle duplicating its
# caption, the defect flagged on the workflow schematic; it is not reproduced.
fig_s4,ax_s4=plt.subplots(1,2,figsize=(11.8,4.9))
fig_f6,ax_f6=plt.subplots(1,2,figsize=(12.6,4.9))
axes={(0,0):ax_s4[0],(1,0):ax_s4[1],(0,1):ax_f6[0],(1,1):ax_f6[1]}
# a transform
ax=axes[0,0]
R=np.linspace(0,8,400); E=R/(1+R)
ax.plot(R,E,color=TEAL,linewidth=2.5)
for r in [0.5,1,2,4]: ax.scatter([r],[r/(1+r)],color=NAVY,zorder=3); ax.annotate(f'R={r:g}\nE={r/(1+r):.2f}',(r,r/(1+r)),xytext=(6,8),textcoords='offset points',fontsize=8.5)
ax.set_xlim(0,8); ax.set_ylim(0,1); ax.set_xlabel('Strongest axis-normalized energetic magnitude, R'); ax.set_ylabel('Energetic component, E = R/(1+R)')
ax.set_title('Soft saturation preserves high-end ranking',loc='left',fontweight='bold',color=NAVY); ax.grid(alpha=.2); panel_label(ax, 'A')
# Fig 6 panel A: score distribution
# b score distribution
ax=axes[0,1]
neg=pv[~pv.structural_ground_truth.astype(bool)].isds_v1
pos=pv[pv.structural_ground_truth.astype(bool)].isds_v1
bins=np.linspace(0,1,11)
_pp=summary['population']
ax.hist(neg,bins=bins,alpha=.75,label=f"No modeled lesion (n={_pp['negative']})",color=BLUE)
ax.hist(pos,bins=bins,alpha=.78,label=f"Structural mechanism (n={_pp['positive']})",color=ORANGE)
ax.set_xlabel('ISDS-v1'); ax.set_ylabel('Variants'); ax.set_title('Higher scores enrich structural mechanisms',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False,fontsize=8.8); panel_label(ax, 'A')
# c ROC/PR
ax=axes[1,0]
y=pv.structural_ground_truth.astype(int).to_numpy()
for col,label,color in [('isds_v1','ISDS-v1',TEAL),('isds_energy_component','Energy component',ORANGE),('isds_context_component','Context component',PURPLE)]:
    fpr,tpr,_=roc_curve(y,pv[col]); auc=float(metrics.loc[metrics.score.eq(col),'roc_auc'].iloc[0]); ax.plot(fpr,tpr,label=f'{label}: AUC {auc:.3f}',color=color,linewidth=2.2)
ax.plot([0,1],[0,1],'--',color=MID)
ax.set_xlabel('False-positive rate'); ax.set_ylabel('True-positive rate'); ax.set_title('ISDS-v1 combines complementary evidence',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False,fontsize=8.8); ax.grid(alpha=.2); panel_label(ax, 'B')
# d top k
ax=axes[1,1]
t=topk[topk.score.eq('isds_v1')]
ax.plot(t.k,t.precision_at_k,'o-',color=TEAL,linewidth=2.2,label='Precision among top k')
ax.plot(t.k,t.recovery_at_k,'s-',color=ORANGE,linewidth=2.2,label=f"Fraction of {_pp['positive']} mechanisms recovered")
# Panel B x axis: 'rank depth', not 'experimental budget'. The budget wording
# read as a claim about how many experiments a lab can afford, which this panel
# says nothing about, and the caption now says rank depth. Keep the two in step.
ax.set_xticks(t.k); ax.set_ylim(0,1.05); ax.set_xlabel('Rank depth, k top-ranked variants'); ax.set_ylabel('Fraction'); ax.set_title('Top-ranked variants concentrate structural mechanisms',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False,fontsize=8.7); ax.grid(alpha=.2); panel_label(ax, 'B')
for _f in (fig_s4,fig_f6): _f.tight_layout(w_pad=2.6)
save(fig_s4,'figureS4_priority_definition.png')
save(fig_f6,'figure6_priority_recovery.png')

# Figure 5: tier components and evidence states
# S3 Fig is cited as a single-panel screen comparison, so it is written as its
# own file. The evidence-state panel that shared this composite presents the
# same counts as S10 Table and is not cited by any caption; it is kept as a
# repository diagnostic rather than a manuscript figure.
fig_s3,ax_s3=plt.subplots(figsize=(6.6,4.8))
fig,axes=plt.subplots(1,2,figsize=(12.0,5.2))
for _t in (ax_s3,axes[0]):
    pass
ax=ax_s3
labels=['Interface status alone','Tier without interface bonus','Full Tier 1-2']
# Screen performance read from the analysis summary, which computes it from the
# tier comparator. These were literals here (and, independently, in
# analyze_isds_v1.py) until v7.7, so the ledger correction left this panel
# plotting 17 positives against a comparator that had 20.
_cb={r['screen']:r for r in summary['component_binary_baselines']}
_ord=['interface_status_alone','tier_without_interface_bonus','full_tier_1_2']
sens=[_cb[k]['sensitivity'] for k in _ord]; spec=[_cb[k]['specificity'] for k in _ord]
x=np.arange(3); width=.35
ax.bar(x-width/2,sens,width,color=ORANGE,label='Sensitivity')
ax.bar(x+width/2,spec,width,color=BLUE,label='Specificity')
ax.set_xticks(x,labels,rotation=12,ha='right'); ax.set_ylim(0,1.08); ax.set_ylabel('Fraction'); ax.set_title('The full tier trades specificity for sensitivity',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False); ax.grid(axis='y',alpha=.2)
for i,(a,b) in enumerate(zip(sens,spec)):
    ax.text(i-width/2,a+.025,f'{a:.3f}',ha='center',fontsize=8.8); ax.text(i+width/2,b+.025,f'{b:.3f}',ha='center',fontsize=8.8)
fig_s3.tight_layout()
save(fig_s3,'figureS3_context_comparators.png')
ax=axes[1]
states=['Convergent','Energy only','Context only','Neither']
_fs=summary['four_state_agreement']
assert _fs is not None, 'four_state_agreement missing from ISDS_v1_summary.json — re-run analyze_isds_v1.py'
struct=np.array(_fs['structural']); neg=np.array(_fs['no_lesion']); x=np.arange(4)
ax.bar(x-.18,struct,.36,color=ORANGE,label='Structural mechanism')
ax.bar(x+.18,neg,.36,color=BLUE,label='No modeled lesion')
ax.set_xticks(x,['Convergent','Energy\nonly','Context\nonly','Neither']); ax.set_ylabel('Variants'); ax.set_ylim(0,18.2); ax.set_title('Discrete evidence states remain interpretable',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False,fontsize=8.8); ax.grid(axis='y',alpha=.2)
for i,(a,b) in enumerate(zip(struct,neg)):
    if a: ax.text(i-.18,a+.25,str(a),ha='center',fontweight='bold',color=ORANGE)
    if b: ax.text(i+.18,b+.25,str(b),ha='center',fontweight='bold',color=BLUE)
ax.text(.02,-.19,'States use the 2.5 kcal/mol mechanism-call threshold; they are not defined by an ISDS cutoff.',transform=ax.transAxes,fontsize=8.8,color=GRAY)
panel_label(ax, 'B')
fig.tight_layout(w_pad=3)
save(fig,'figure5_context_components_and_states.png')

# Figure 6: threshold tradeoff (panel a/b); existing calibration can remain separate supplementary/main depending layout
fig,axes=plt.subplots(1,2,figsize=(12,5.0))
# Panel A was hardcoded through v7.8 and had gone stale at every point
# (recovery .760/.720/.640/.600/.540 against the committed .714/.679/.607/
# .554/.500). Read the generator output. Panel B's counts come from
# build_measured_effect_tables.py and are asserted against it below.
_top = pd.read_csv(_REPO / 'reference_outputs' / 'COMAVI_threshold_operating_points.csv')
_top = _top.set_index('threshold').loc[['t10', 't15', 't20', 't25', 'tSAP']]
thresholds=['1.0','1.5','2.0','2.5','2.9/2.9/3.5']
recovery=[float(v) for v in _top.sensitivity]; rejection=[float(v) for v in _top.specificity]
assert len(recovery) == 5 and len(rejection) == 5
ax=axes[0]
x=np.arange(5)
ax.plot(x,recovery,'o-',color=ORANGE,linewidth=2.4,label='Mechanism recovery')
ax.plot(x,rejection,'s-',color=BLUE,linewidth=2.4,label='Correct rejection')
ax.set_xticks(x,thresholds); ax.set_ylim(.35,.9); ax.set_xlabel('Decision threshold (kcal/mol)'); ax.set_ylabel('Fraction'); ax.set_title('Threshold choice changes the error tradeoff',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False); ax.grid(alpha=.2); panel_label(ax, 'A')
ax=axes[1]
_met = json.loads((_REPO / 'reference_outputs' / 'COMAVI_measured_effect_tables.json').read_text())
measured=[int(r['cell'].split('/')[0]) for r in _met['measured_effects_recovered']['by_threshold']]
_mdenom=int(_met['measured_effects_recovered']['by_threshold'][0]['cell'].split('/')[1])
assert len(measured) == 5, measured
ax.bar(x,measured,color=[GREEN,GREEN,GOLD,NAVY,PURPLE])
ax.set_xticks(x,thresholds); ax.set_ylim(0,_mdenom); ax.set_xlabel('Decision threshold (kcal/mol)'); ax.set_ylabel(f'Measured destabilizing effects recovered (of {_mdenom})'); ax.set_title('Higher thresholds miss more measured effects',loc='left',fontweight='bold',color=NAVY); ax.grid(axis='y',alpha=.2)
for i,v in enumerate(measured): ax.text(i,v+.8,f'{v}/{_mdenom}',ha='center',fontweight='bold',fontsize=9)
panel_label(ax, 'B')
fig.tight_layout(w_pad=3)
save(fig,'figure6_threshold_tradeoff.png')

# Figure 7: AlphaMissense, ISDS, mechanism
# Two rows, not one. At 14.2 in wide on a 6.5 in text column every label was
# downscaled by 0.46, so 10.5 pt type reached the page at under 5 pt. Narrowing
# the figure and raising the type is the only thing that makes it legible; a
# wider figure with the same font makes it worse.
clin=pd.read_csv(AN/'ISDS_v1_alphamissense_common_set.csv')
with plt.rc_context({'font.size':12.5,'axes.titlesize':13.0,'axes.labelsize':12.5,
                     'xtick.labelsize':11.5,'ytick.labelsize':11.5}):
    fig=plt.figure(figsize=(9.4,7.6))
    gs=fig.add_gridspec(2,2,height_ratios=[1.0,0.82],hspace=0.40,wspace=0.28)
    axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,:])
    colors=np.where(clin.clinical_y.eq(1),ORANGE,BLUE)
    for ax,ycol,ylab,ttl,lab in [
            (axA,'isds_v1','ISDS-v1','Pathogenicity and priority differ','A'),
            (axB,'isds_energy_component','ISDS energetic component','Energetic evidence is one component','B')]:
        ax.scatter(clin['AM pathogenicity'],clin[ycol],c=colors,edgecolor='white',linewidth=.5,s=52)
        ax.axvspan(.34,.564,color=MID,alpha=.45)
        ax.set_xlabel('AlphaMissense pathogenicity score'); ax.set_ylabel(ylab)
        ax.set_title(ttl,loc='left',fontweight='bold',color=NAVY)
        ax.grid(alpha=.15); panel_label(ax,lab)
    labels=['AlphaMissense','ISDS-v1','Energy','Context']; vals=[.902778,.851010,.835859,.768939]
    cols=[NAVY,TEAL,ORANGE,PURPLE]
    bars=axC.barh(np.arange(4),vals,color=cols)
    axC.set_yticks(np.arange(4),labels); axC.invert_yaxis(); axC.set_xlim(.5,1.0)
    axC.set_xlabel('Pathogenicity AUC (descriptive)')
    axC.set_title('COMAVI is not a replacement pathogenicity model',loc='left',fontweight='bold',color=NAVY)
    axC.grid(axis='x',alpha=.2)
    for b,v in zip(bars,vals):
        axC.text(v+.007,b.get_y()+b.get_height()/2,f'{v:.3f}',va='center',fontsize=11.5,fontweight='bold')
    panel_label(axC,'C')
    fig.text(.012,.012,'Orange: pathogenic or pathogenic gain-of-function; blue: benign. '
             'The grey band marks the AlphaMissense ambiguous range.',fontsize=10.5,color=GRAY)
    fig.subplots_adjust(left=0.105,right=0.985,top=0.945,bottom=0.095)
save(fig,'figure7_alphamissense_isds_mechanism.png')

# Alt-text file
alt = {
'figure1_unified_comavi_workflow.png':'Workflow diagram showing a missense variant and structure feeding four native COMAVI calculations. The outputs branch into ISDS-v1 for prioritization and a signed multi-axis mechanism profile for assay selection. Pathogenicity remains separate.',
'figure2_population_map.png':'Population flow from a 61-variant resource to 57 mechanism-gradeable variants and 47 tier-carrying structural-prioritization variants, with four ungraded scope or commitment cases shown separately.',
'figure3_mechanism_localization.png':'Two horizontal bar charts showing direction-aware agreement for monomer, complex-context, binding, and all energetic axes, and whole-variant mechanism-pattern scores for the full, interaction, and BRCT populations.',
'figure6_priority_recovery.png':'Two panels showing ISDS-v1 score distributions for structural-mechanism and no-lesion variants, and precision and recovery as the experimental budget increases.',
'figureS4_priority_definition.png':'Two panels showing the soft-saturating energetic transformation and structural-prioritization ROC curves for ISDS-v1 and its energy and context components.',
'figureS3_context_comparators.png':'Sensitivity and specificity of three structural-context screens: interface status alone, the tier without its interface bonus, and the full Tier 1-2 screen.',
'figure5_context_components_and_states.png':'Two panels comparing sensitivity and specificity of interface, no-interface tier, and full tier screens, and counts of structural mechanisms and negatives across convergent, energy-only, context-only, and neither states.',
'figure6_threshold_tradeoff.png':'Two panels showing mechanism recovery falling and correct rejection rising across five decision thresholds, and the corresponding decline in measured destabilizing effects recovered.',
'figure7_alphamissense_isds_mechanism.png':'Two scatterplots compare AlphaMissense score with ISDS and its energy component; a horizontal bar chart shows descriptive pathogenicity AUCs for AlphaMissense, ISDS, energy, and context.'
}
(OUT/'ALT_TEXT.json').write_text(json.dumps(alt,indent=2)+'\n')
print('wrote figures to',OUT)
