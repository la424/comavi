from __future__ import annotations

import argparse
import json
from pathlib import Path

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


def box(ax, xy, wh, title, body, fc, ec, fontsize=9.2):
    x, y = xy; w, h = wh
    p = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.012,rounding_size=0.02', facecolor=fc, edgecolor=ec, linewidth=1.4)
    ax.add_patch(p)
    ax.text(x+w/2, y+h*0.67, title, ha='center', va='center', fontsize=10.2, fontweight='bold', color=DARK)
    ax.text(x+w/2, y+h*0.31, body, ha='center', va='center', fontsize=fontsize, color=DARK, linespacing=1.25)
    return p


def arrow(ax, x1, y1, x2, y2, color=GRAY):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=14,linewidth=1.3,color=color))


# Figure 1: unified workflow
fig, ax = plt.subplots(figsize=(12.0, 6.8))
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
ax.text(0.02, 0.96, 'COMAVI returns a priority score and a mechanism profile from the same structural calculations', fontsize=16, fontweight='bold', color=NAVY, va='top')
ax.text(0.02, 0.90, 'ISDS-v1 ranks candidates for structural follow-up; the decomposed profile identifies the proposed physical lesion and assay.', fontsize=10.8, color=GRAY, va='top')
box(ax, (0.02,0.36),(0.12,0.25),'INPUT','Missense variant\n+ structure',LIGHT,NAVY)
arrow(ax,0.14,0.485,0.20,0.485)
# calculation boxes
ys=[0.68,0.49,0.30,0.11]
colors=[('#EAF3F8',BLUE),('#E9F6F2',TEAL),('#FFF1E8',ORANGE),('#F1EDFA',PURPLE)]
titles=['Monomer-fold ΔΔG','Complex-context ΔΔG','Binding-interface ΔΔG','Structural-context tier']
bodies=['FoldX on isolated subunit','FoldX on assembled coordinates','FoldX interaction energy by partner','Chemistry, contacts, interface, burial, confidence']
for y,(fc,ec),t,b in zip(ys,colors,titles,bodies):
    box(ax,(0.20,y),(0.23,0.14),t,b,fc,ec,8.5)
# arrows to two outputs
for y in [0.75,0.56,0.37]: arrow(ax,0.43,y,0.52,0.66,BLUE)
arrow(ax,0.43,0.18,0.52,0.32,PURPLE)
box(ax,(0.52,0.54),(0.20,0.25),'ISDS-v1','Cohort-independent\nstructural-disruption\npriority score', '#EAF6F5', TEAL, 9.0)
box(ax,(0.52,0.18),(0.20,0.25),'Mechanism profile','Signed values and calls for\nmonomer fold, complex context,\nand binding', '#EDF3FA', BLUE, 8.7)
arrow(ax,0.72,0.66,0.79,0.66,TEAL)
arrow(ax,0.72,0.30,0.79,0.30,BLUE)
box(ax,(0.79,0.54),(0.19,0.25),'PRIORITIZE','Rank variants for\nstructural follow-up', '#EEF8F1', GREEN, 9.2)
box(ax,(0.79,0.18),(0.19,0.25),'LOCALIZE','Select stability, assembly,\nor interaction experiments', '#FFF7E7', GOLD, 9.0)
ax.text(0.53,0.08,'Pathogenicity evidence remains separate. ISDS-v1 is not a probability, and no binary ISDS cutoff is validated.',fontsize=9.4,color=RED,fontweight='bold')
save(fig,'figure1_unified_comavi_workflow.png')

# Figure 2: population map
fig, ax = plt.subplots(figsize=(11.5, 5.6))
ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
ax.text(0.02,0.94,'The resource supports two complementary evaluation tasks',fontsize=16,fontweight='bold',color=NAVY,va='top')
box(ax,(0.04,0.35),(0.20,0.35),'61-variant resource','14 protein systems\n49 interaction variants\n12 BRCT variants',LIGHT,NAVY,9.1)
arrow(ax,0.24,0.52,0.34,0.68)
arrow(ax,0.24,0.52,0.34,0.34)
box(ax,(0.34,0.55),(0.25,0.27),'57 mechanism-gradeable','Whole-variant mechanism-pattern score\nThree energetic axes\nInteraction + BRCT systems','#EAF3F8',BLUE,8.8)
box(ax,(0.34,0.20),(0.25,0.27),'4 ungraded resource cases','2 known mechanisms missing required context\n2 mechanism-uncommitted reference cases','#F5F5F5',GRAY,8.7)
arrow(ax,0.59,0.69,0.70,0.69)
box(ax,(0.70,0.55),(0.26,0.27),'47 structural-prioritization','Tier-carrying interaction variants\n17 modeled structural mechanisms\n30 variants with no committed modeled lesion','#E9F6F2',TEAL,8.8)
ax.text(0.34,0.12,'Other analysis subsets are defined once and referenced consistently: 15 direct-energy comparators, 44 measured destabilizers, and 47 AlphaMissense complete cases.',fontsize=9.3,color=GRAY)
save(fig,'figure2_population_map.png')

# Figure 3: mechanism localization performance
summary=json.load(open(AN/'ISDS_v1_summary.json'))
fig, axes=plt.subplots(1,2,figsize=(11.8,5.2),gridspec_kw={'width_ratios':[1.2,1]})
ax=axes[0]
labels=['Monomer fold','Complex context','Binding','All energetic axes']
vals=[21/27,20/26,24/32,65/85]
nums=['21/27','20/26','24/32','65/85']
y=np.arange(len(labels))
bars=ax.barh(y,vals,color=[BLUE,TEAL,ORANGE,NAVY],height=.58)
ax.set_yticks(y,labels); ax.invert_yaxis(); ax.set_xlim(0,1.03); ax.set_xlabel('Direction-aware agreement')
ax.set_title('Physical axes remain separate',loc='left',fontweight='bold',color=NAVY)
ax.grid(axis='x',alpha=.2)
for bar,val,n in zip(bars,vals,nums): ax.text(val+.015,bar.get_y()+bar.get_height()/2,f'{n}  ({val:.3f})',va='center',fontsize=9.5,fontweight='bold')
panel_label(ax,'a')
ax=axes[1]
vals=[41/57,34/47,7/10]
labels=['All gradeable variants','Interaction subset','BRCT fold subset']
colors=[NAVY,TEAL,PURPLE]
y=np.arange(3)
bars=ax.barh(y,vals,color=colors,height=.58)
ax.set_yticks(y,labels); ax.invert_yaxis(); ax.set_xlim(0,1.03); ax.set_xlabel('Whole-variant pattern score')
ax.set_title('Mechanism localization across systems',loc='left',fontweight='bold',color=NAVY)
ax.grid(axis='x',alpha=.2)
for bar,val,n in zip(bars,vals,['41/57','34/47','7/10']): ax.text(val+.015,bar.get_y()+bar.get_height()/2,f'{n}  ({val:.3f})',va='center',fontsize=9.5,fontweight='bold')
ax.text(0.02,-0.22,'System-cluster 95% CI for the primary score: 0.625-0.821',transform=ax.transAxes,fontsize=9.2,color=GRAY)
panel_label(ax,'b')
fig.tight_layout(w_pad=3.0)
save(fig,'figure3_mechanism_localization.png')

# Figure 4: ISDS definition and performance
pv=pd.read_csv(AN/'ISDS_v1_per_variant.csv')
metrics=pd.read_csv(AN/'ISDS_v1_primary_metrics.csv')
topk=pd.read_csv(AN/'ISDS_v1_top_k.csv')
boot=pd.read_csv(AN/'ISDS_v1_system_cluster_bootstrap_summary.csv')
fig,axes=plt.subplots(2,2,figsize=(12.2,9.2))
# a transform
ax=axes[0,0]
R=np.linspace(0,8,400); E=R/(1+R)
ax.plot(R,E,color=TEAL,linewidth=2.5)
for r in [0.5,1,2,4]: ax.scatter([r],[r/(1+r)],color=NAVY,zorder=3); ax.annotate(f'R={r:g}\nE={r/(1+r):.2f}',(r,r/(1+r)),xytext=(6,8),textcoords='offset points',fontsize=8.5)
ax.set_xlim(0,8); ax.set_ylim(0,1); ax.set_xlabel('Strongest axis-normalized energetic magnitude, R'); ax.set_ylabel('Energetic component, E = R/(1+R)')
ax.set_title('Soft saturation preserves high-end ranking',loc='left',fontweight='bold',color=NAVY); ax.grid(alpha=.2); panel_label(ax,'a')
# b score distribution
ax=axes[0,1]
neg=pv[~pv.structural_ground_truth.astype(bool)].isds_v1
pos=pv[pv.structural_ground_truth.astype(bool)].isds_v1
bins=np.linspace(0,1,11)
ax.hist(neg,bins=bins,alpha=.75,label='No modeled lesion (n=30)',color=BLUE)
ax.hist(pos,bins=bins,alpha=.78,label='Structural mechanism (n=17)',color=ORANGE)
ax.set_xlabel('ISDS-v1'); ax.set_ylabel('Variants'); ax.set_title('Higher scores enrich modeled structural mechanisms',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False,fontsize=8.8); panel_label(ax,'b')
# c ROC/PR
ax=axes[1,0]
y=pv.structural_ground_truth.astype(int).to_numpy()
for col,label,color in [('isds_v1','ISDS-v1',TEAL),('isds_energy_component','Energy component',ORANGE),('isds_context_component','Context component',PURPLE)]:
    fpr,tpr,_=roc_curve(y,pv[col]); auc=float(metrics.loc[metrics.score.eq(col),'roc_auc'].iloc[0]); ax.plot(fpr,tpr,label=f'{label}: AUC {auc:.3f}',color=color,linewidth=2.2)
ax.plot([0,1],[0,1],'--',color=MID)
ax.set_xlabel('False-positive rate'); ax.set_ylabel('True-positive rate'); ax.set_title('ISDS-v1 combines complementary evidence',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False,fontsize=8.8); ax.grid(alpha=.2); panel_label(ax,'c')
# d top k
ax=axes[1,1]
t=topk[topk.score.eq('isds_v1')]
ax.plot(t.k,t.precision_at_k,'o-',color=TEAL,linewidth=2.2,label='Precision among top k')
ax.plot(t.k,t.recovery_at_k,'s-',color=ORANGE,linewidth=2.2,label='Fraction of 17 mechanisms recovered')
ax.set_xticks(t.k); ax.set_ylim(0,1.05); ax.set_xlabel('Experimental budget, k variants'); ax.set_ylabel('Fraction'); ax.set_title('Top-ranked variants concentrate structural mechanisms',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False,fontsize=8.7); ax.grid(alpha=.2); panel_label(ax,'d')
fig.suptitle('ISDS-v1 is a fixed prioritization index, not a probability or validated binary classifier',fontsize=15.5,fontweight='bold',color=NAVY,y=.995)
fig.tight_layout(rect=[0,0,1,.97],h_pad=3,w_pad=2.5)
save(fig,'figure4_isds_definition_performance.png')

# Figure 5: tier components and evidence states
fig,axes=plt.subplots(1,2,figsize=(12.0,5.2))
ax=axes[0]
labels=['Interface status alone','Tier without interface bonus','Full Tier 1-2']
sens=[15/17,14/17,1.0]; spec=[23/30,18/30,17/30]
x=np.arange(3); width=.35
ax.bar(x-width/2,sens,width,color=ORANGE,label='Sensitivity')
ax.bar(x+width/2,spec,width,color=BLUE,label='Specificity')
ax.set_xticks(x,labels,rotation=12,ha='right'); ax.set_ylim(0,1.08); ax.set_ylabel('Fraction'); ax.set_title('The full tier trades specificity for sensitivity',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False); ax.grid(axis='y',alpha=.2)
for i,(a,b) in enumerate(zip(sens,spec)):
    ax.text(i-width/2,a+.025,f'{a:.3f}',ha='center',fontsize=8.8); ax.text(i+width/2,b+.025,f'{b:.3f}',ha='center',fontsize=8.8)
panel_label(ax,'a')
ax=axes[1]
states=['Convergent','Energy only','Context only','Neither']
struct=np.array([15,0,2,0]); neg=np.array([3,3,10,14]); x=np.arange(4)
ax.bar(x-.18,struct,.36,color=ORANGE,label='Structural mechanism')
ax.bar(x+.18,neg,.36,color=BLUE,label='No modeled lesion')
ax.set_xticks(x,['Convergent','Energy\nonly','Context\nonly','Neither']); ax.set_ylabel('Variants'); ax.set_ylim(0,16.5); ax.set_title('Discrete evidence states remain interpretable',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False,fontsize=8.8); ax.grid(axis='y',alpha=.2)
for i,(a,b) in enumerate(zip(struct,neg)):
    if a: ax.text(i-.18,a+.25,str(a),ha='center',fontweight='bold',color=ORANGE)
    if b: ax.text(i+.18,b+.25,str(b),ha='center',fontweight='bold',color=BLUE)
ax.text(.02,-.19,'States use the 2.5 kcal/mol mechanism-call threshold; they are not defined by an ISDS cutoff.',transform=ax.transAxes,fontsize=8.8,color=GRAY)
panel_label(ax,'b')
fig.tight_layout(w_pad=3)
save(fig,'figure5_context_components_and_states.png')

# Figure 6: threshold tradeoff (panel a/b); existing calibration can remain separate supplementary/main depending layout
fig,axes=plt.subplots(1,2,figsize=(12,5.0))
thresholds=['1.0','1.5','2.0','2.5','2.9/2.9/3.5']
recovery=[.760,.720,.640,.600,.540]; rejection=[.438,.594,.719,.812,.844]
ax=axes[0]
x=np.arange(5)
ax.plot(x,recovery,'o-',color=ORANGE,linewidth=2.4,label='Mechanism recovery')
ax.plot(x,rejection,'s-',color=BLUE,linewidth=2.4,label='Correct rejection')
ax.set_xticks(x,thresholds); ax.set_ylim(.35,.9); ax.set_xlabel('Decision threshold (kcal/mol)'); ax.set_ylabel('Fraction'); ax.set_title('Threshold choice changes the error tradeoff',loc='left',fontweight='bold',color=NAVY); ax.legend(frameon=False); ax.grid(alpha=.2); panel_label(ax,'a')
ax=axes[1]
measured=[35,26,19,13,12]
ax.bar(x,measured,color=[GREEN,GREEN,GOLD,NAVY,PURPLE])
ax.set_xticks(x,thresholds); ax.set_ylim(0,44); ax.set_xlabel('Decision threshold (kcal/mol)'); ax.set_ylabel('Measured destabilizing effects recovered (of 44)'); ax.set_title('Higher thresholds miss more measured effects',loc='left',fontweight='bold',color=NAVY); ax.grid(axis='y',alpha=.2)
for i,v in enumerate(measured): ax.text(i,v+.8,f'{v}/44',ha='center',fontweight='bold',fontsize=9)
panel_label(ax,'b')
fig.tight_layout(w_pad=3)
save(fig,'figure6_threshold_tradeoff.png')

# Figure 7: AlphaMissense, ISDS, mechanism
clin=pd.read_csv(AN/'ISDS_v1_alphamissense_common_set.csv')
fig,axes=plt.subplots(1,3,figsize=(14.2,4.7))
colors=np.where(clin.clinical_y.eq(1),ORANGE,BLUE)
ax=axes[0]; ax.scatter(clin['AM pathogenicity'],clin['isds_v1'],c=colors,edgecolor='white',linewidth=.5,s=48); ax.axvspan(.34,.564,color=MID,alpha=.45); ax.set_xlabel('AlphaMissense pathogenicity score'); ax.set_ylabel('ISDS-v1'); ax.set_title('Pathogenicity and structural priority differ',loc='left',fontweight='bold',color=NAVY); ax.grid(alpha=.15); panel_label(ax,'a')
ax=axes[1]; ax.scatter(clin['AM pathogenicity'],clin['isds_energy_component'],c=colors,edgecolor='white',linewidth=.5,s=48); ax.axvspan(.34,.564,color=MID,alpha=.45); ax.set_xlabel('AlphaMissense pathogenicity score'); ax.set_ylabel('ISDS energetic component'); ax.set_title('Energetic evidence is one component',loc='left',fontweight='bold',color=NAVY); ax.grid(alpha=.15); panel_label(ax,'b')
ax=axes[2]
labels=['AlphaMissense','ISDS-v1','Energy','Context']; vals=[.902778,.851010,.835859,.768939]; cols=[NAVY,TEAL,ORANGE,PURPLE]
bars=ax.barh(np.arange(4),vals,color=cols); ax.set_yticks(np.arange(4),labels); ax.invert_yaxis(); ax.set_xlim(.5,1.0); ax.set_xlabel('Pathogenicity AUC (descriptive)'); ax.set_title('COMAVI is not a replacement pathogenicity model',loc='left',fontweight='bold',color=NAVY); ax.grid(axis='x',alpha=.2)
for b,v in zip(bars,vals): ax.text(v+.008,b.get_y()+b.get_height()/2,f'{v:.3f}',va='center',fontsize=9,fontweight='bold')
panel_label(ax,'c')
fig.text(.02,.01,'Orange: pathogenic/pathogenic gain-of-function; blue: benign. The gray band marks the AlphaMissense ambiguous range.',fontsize=9,color=GRAY)
fig.tight_layout(rect=[0,.035,1,1],w_pad=2.4)
save(fig,'figure7_alphamissense_isds_mechanism.png')

# Alt-text file
alt = {
'figure1_unified_comavi_workflow.png':'Workflow diagram showing a missense variant and structure feeding four native COMAVI calculations. The outputs branch into ISDS-v1 for prioritization and a signed multi-axis mechanism profile for assay selection. Pathogenicity remains separate.',
'figure2_population_map.png':'Population flow from a 61-variant resource to 57 mechanism-gradeable variants and 47 tier-carrying structural-prioritization variants, with four ungraded scope or commitment cases shown separately.',
'figure3_mechanism_localization.png':'Two horizontal bar charts showing direction-aware agreement for monomer, complex-context, binding, and all energetic axes, and whole-variant mechanism-pattern scores for the full, interaction, and BRCT populations.',
'figure4_isds_definition_performance.png':'Four panels showing the soft-saturating energetic transformation, ISDS distributions by structural-mechanism class, ROC curves for ISDS and its components, and top-k precision and recovery.',
'figure5_context_components_and_states.png':'Two panels comparing sensitivity and specificity of interface, no-interface tier, and full tier screens, and counts of structural mechanisms and negatives across convergent, energy-only, context-only, and neither states.',
'figure6_threshold_tradeoff.png':'Two panels showing mechanism recovery falling and correct rejection rising across five decision thresholds, and the corresponding decline in measured destabilizing effects recovered.',
'figure7_alphamissense_isds_mechanism.png':'Two scatterplots compare AlphaMissense score with ISDS and its energy component; a horizontal bar chart shows descriptive pathogenicity AUCs for AlphaMissense, ISDS, energy, and context.'
}
(OUT/'ALT_TEXT.json').write_text(json.dumps(alt,indent=2)+'\n')
print('wrote figures to',OUT)
