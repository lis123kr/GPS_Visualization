from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
from mcolors import mcols

filename = '..\\N_백신주 ORF62.xlsx'
basepair = ['A', 'G', 'C', 'T']
duma_position = 'Du_position'
mcolor = [ 'Blues', 'Orgs', 'Grns', 'Reds' ]
dcolor = [ [None, 'BuOg', 'BuGn', 'BuRd' ], ['OgBu', None, 'OgGn', 'OgRd'], 
		['GnBu', 'GnOg', None, 'GnRd'], ['RdBu', 'RdOg', 'RdGn', None] ]


def get_major_minor(Passage):	
	from numpy import zeros
	from pandas import DataFrame
	"""	
		- major : 염기 A, G, C, T 중 가장 큰 값
		- minor : 염기 중 두번째로 큰 값, argsort 후 인덱스 2의 값을 추출
	"""
	ranked_df = Passage.values.argsort(axis=1)		
	arr = zeros((len(Passage),1))
	for i, d in enumerate(ranked_df[:,2]):
		arr[i] = Passage[i:i+1][basepair[d]]
	major = DataFrame(Passage.max(axis=1))
	minor = DataFrame(arr, index=major.index)		
	return major, minor

def preprocessing(bp):
	from pandas import DataFrame
	from numpy import divide

	bp['sum'] = bp[basepair].sum(axis=1)
	bp = bp[ bp['sum'] >= 35.0 ]

	_, bp['minor'] = get_major_minor(bp[basepair])

	argx = -bp[basepair]
	bps = argx.values.argsort(axis=1)
	bp['Mj_seq'] = DataFrame(bps[:,0], index=bp.index)
	bp['Mn_seq'] = DataFrame(bps[:,1], index=bp.index)

	bp['maf'] = divide(bp[['minor']], bp[['sum']]) * 100

	return bp[[duma_position] + ['maf', 'Mj_seq', 'Mn_seq']]

if __name__ == '__main__':
	import matplotlib.patches as mpatches
	from pandas import ExcelFile, DataFrame, merge
	from numpy import linspace, logical_and, logical_or, array, abs
	xls = ExcelFile(filename)
	sheet = [ preprocessing(xls.parse(name)) for name in xls.sheet_names ]
	nsheet = len(sheet)

	fig, ax = plt.subplots(3, 1, sharex=True)
	fig.set_size_inches(12, 7, forward=True)
	plt.subplots_adjust(left=0.16, bottom=0.05, right=0.98, hspace=0.15)

	bpdiff = []
	bpdiff_trans = [] # major transition

	colors = ['blue', 'orange', 'lightgreen', 'red']
	spatches = [ mpatches.Patch(color=colors[idx], label=basepair[idx]) for idx in range(4) ]

	fig.legend(handles=spatches, labels=basepair, bbox_to_anchor=(0.33, .90, 0.65, .102), loc=3,
           ncol=4, mode="expand", borderaxespad=0.)

	for i in range(nsheet):
		for j in range(nsheet):
			if i >= j:
				continue

			merged = merge(sheet[i], sheet[j], how='outer', on=duma_position, left_index=True, right_index=True)

			# 공통 
			cond = abs(merged['maf_x'] - merged['maf_y']) >= 5.0
			cond1 = merged['Mj_seq_x'] == merged['Mj_seq_x']
			cond2 = merged['Mj_seq_x'] != merged['Mj_seq_x']
			bpdiff.append(merged[ logical_and(cond, cond1).values ])
			bpdiff_trans.append( merged[ logical_and(cond, cond2).values ])

	# from matplotlib import cm
	seq_size = 100
	x = np.linspace(0.0, 1.0, seq_size)
	ylabel = ['low-middle', 'low-high', 'middle-high']
	for i in range(len(bpdiff)):
		for a in range(4):

			cond1 = bpdiff[i]['Mj_seq_x'] == a

			for b in range(4):
				if a == b:
					continue
				# major a -> minor b				

				# 1. minor가 같은 것들 처리 : 같은 색상
				cond2 = logical_and(bpdiff[i]['Mn_seq_x'] == b, bpdiff[i]['Mn_seq_y'] == b)

				d = bpdiff[i][ logical_and(cond1, cond2).values ]

				y_ = [ linspace(d_.maf_x, d_.maf_y, seq_size) for d_ in d[['maf_x', 'maf_y']].itertuples() ]

				x_ = array(y_)
				for j in range(len(d)):
					x_[j].fill( d[duma_position].values[j] )
				c = [ x for _ in range(len(x_))]

				ax[i].scatter(x_, y_, c=c, cmap=mcols.cmap(mcolor[b]), s=10, linewidth=0.0) # 
				
				# # 2. minor가 바뀌는 것들 처리 : colors diverging
				# 함수로 각각 처리해야 할 듯

				# cond3 = logical_and(bpdiff[i]['Mn_seq_x'] == b, bpdiff[i]['Mn_seq_y'] != b)
				# d = bpdiff[i][ logical_and(cond1, cond3).values ]

				# y_ = [ linspace(d_.maf_x, d_.maf_y, seq_size) for d_ in d[['maf_x', 'maf_y']].itertuples() ]				

				# x_ = array(y_)
				# for j in range(len(d)):
				# 	x_[j].fill( d[duma_position].values[j] )
				# c = [ x for _ in range(len(x_))]

				# ax[i].scatter(x_, y_, c=c, cmap=mcols.cmap(dcolor[b][]))
		ax[i].set_ylim([0, 50])
		ax[i].set_ylabel(ylabel[i])







	
	



	# c_ = x
	# cmap = 'viridis'
	# y_ = np.linspace(0, 10, 100)
	# ax.scatter(x, y_, c=c_, cmap=cmap, s=300, linewidth=0.0)

	plt.show()