from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

filename = 'YC01_ORF_616263.xlsx'
basepair = ['A', 'G', 'C', 'T']
duma_position = 'Dm_position'

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
	bp['Mn_seq'] = DataFrame(bps[:,1], index=[bp.index])

	bp['maf'] = divide(bp[['minor']], bp[['sum']]) * 100

	return bp[[duma_position] + basepair + ['maf', 'Mn_seq']]

if __name__ == '__main__':
	from pandas import ExcelFile, DataFrame, merge

	# xls = ExcelFile(filename)
	# sheet_ = [ preprocessing(xls.parse(name)) for name in xls.sheet_names ]

	x = np.linspace(0.0, 1.0, 100)
	from matplotlib import cm
	# (x)
	print(cm.get_cmap('viridis'))

	fig, ax = plt.subplots()

	c_ = x
	cmap = 'viridis'
	y_ = np.arange(0, 10, 100)
	ax.scatter(x, y_, c=c_, cmap=cmap, s=300, linewidth=0.0)
