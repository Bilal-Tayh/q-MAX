#!/usr/bin/python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import sys

fig, ax = plt.subplots(1,1)
bar_cycle = cycler(color=list('bgrcmyw')) + cycler(hatch=['--', '++', 'xX', '\\\\', '**', '..', 'OO'])
styles = bar_cycle()
ax.set_ylabel('CPU Throughput [MPPS]', fontsize=20)
ax.yaxis.grid(True)

hash = "xxh"
st = ""
df_c = pd.read_csv('../results/hh_caida.'+hash+'.raw_res', skipinitialspace=True)
df_univ = pd.read_csv('../results/hh_univ1.'+hash+'.raw_res', skipinitialspace=True)
df_c18 = pd.read_csv('../results/hh_caida18.'+hash+'.raw_res', skipinitialspace=True)

df_org = pd.concat([df_c, df_univ, df_c18])
df_org = df_org.replace(np.nan, 0)
df_org['dataset'] = df_org.apply(lambda row: 'caida16' if row['dataset'] == 'caida' else row['dataset'], axis=1)

if st=='qurey':
	ops=df_org['size']
else:
	ops=df_org['insertions']

df_org['throughput'] = ops / 1000000 / df_org['duration']
del df_org['insertions']
del df_org['duration']

sizes=map(lambda s : 10**int(s), sys.argv[1:])
for s in sizes:
        df_org = df_org[df_org['size'] == s]
	del df_org['size']

	df_org['type'] = df_org.apply(
	    lambda row: r'$q$-MAX,$\gamma=$' + str(row['gamma']) if row['type'] == 'AmortizedQMax'+st
	    else 'Heap' if row['type'] == 'Heap'+st else 'SkipList' if row['type'] == 'SkipList'+st else '', axis=1)
	df_org = df_org[df_org['type'] != '']

	j = df_org.groupby(['type', 'dataset'])['throughput'].mean().unstack(0)
	jstd = df_org.groupby(['type', 'dataset'])['throughput'].sem().unstack(0)
	traces = ["caida16", "univ1", "caida18"]

	w=-0.6
	for c in ["Heap", "SkipList", r'$q$-MAX,$\gamma=$0.05', r'$q$-MAX,$\gamma=$0.1', r'$q$-MAX,$\gamma=$0.25']:
	    print c
	    print j[c]
	    y = []
	    yerr = []
	    for t in traces:
	        y.append(j[c][t])
	        yerr.append(jstd[c][t])
	    ax.bar([x+w for x in [1,3,5]],y, yerr=yerr, error_kw=dict(ecolor='r', lw=3, capsize=7, capthick=1), width=0.2, alpha=1, edgecolor='k', linewidth=1, label=c, **next(styles))
	    w+=0.2

	ax.set_xticks([1,3,5])
	ax.set_xticklabels([x.upper() for x in traces], fontsize=16)
	plt.tick_params(axis='y', which='major', labelsize=16)
#	ax.legend(prop={'size': 16}, loc='lower center', framealpha=1.0, bbox_to_anchor=(0.5,1), ncol=2)
	fig.tight_layout()
#	plt.show()
	plt.ylim(ymin=0, ymax=25)
        plt.savefig('../figures/graph6-q='+str(sizes)+'.pdf', bbox_inches='tight')

