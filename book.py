"""
	Created by ilseop, Lee on 2018. 8. 28..
"""

from vizlogging import infolog, errorlog
from numpy import logical_and, logical_or
class Book:

	def __init__(self, filename, params):
		self.duma_position = params['duma_position']
		self.duma_Genome = params['duma_Genome']
		self.duma_Repeat = params['duma_Repeat']
		self.duma_ORF = params['duma_ORF']

		self.basepair = ['A', 'G', 'C', 'T']
		self.dumas_col = [ self.duma_position, self.duma_Genome, self.duma_Repeat, self.duma_ORF ]
		self.constraints = 35.0

		self.filename = filename

		self.sheet = []
		self.dumas_info = { self.duma_Genome : [], self.duma_Repeat: [], self.duma_ORF : []}
		self.range_dumas = { self.duma_Genome : dict(), self.duma_Repeat : dict(), self.duma_ORF : dict()}

		self.minpos, self.maxpos = 1e9, 0
		self.nsheet = 0
		self.xls = None

		self.dcolor = [ [None, 'BuOg', 'BuGn', 'BuRd' ], 
				 ['OgBu', None, 'OgGn', 'OgRd'], 
				 ['GnBu', 'GnOg', None, 'GnRd'], 
				 ['RdBu', 'RdOg', 'RdGn', None] ]
		self.colors = ['blue', 'orange', 'lightgreen', 'red']

	def preprocessing(self):
		infolog("Reading files...")
		self.readfiles_()
		infolog("Done...")

		infolog("preprocessing data...")
		self.preprocessing_()
		self.extract_dumas_info_()
		self.extract_dumas_range_()
		self.extract_incdec()
		infolog("Done...")

	def readfiles_(self):
		from pandas import ExcelFile
		self.xls = ExcelFile(self.filename)
		self.nsheet = len(self.xls.sheet_names)
		

	def get_major_minor(self, Passage):	
		from numpy import zeros
		from pandas import DataFrame
		"""	
			- major : 염기 A, G, C, T 중 가장 큰 값
			- minor : 염기 중 두번째로 큰 값, argsort 후 인덱스 2의 값을 추출
		"""
		ranked_df = Passage.values.argsort(axis=1)		
		arr = zeros((len(Passage),1))
		for i, d in enumerate(ranked_df[:,2]):
			arr[i] = Passage[i:i+1][self.basepair[d]]
		major = DataFrame(Passage.max(axis=1))
		minor = DataFrame(arr, index=major.index)
		return major, minor

	def preprocessing_(self):
		from pandas import DataFrame
		from numpy import divide

		assert self.xls is not None, "read file first"
		
		for name in self.xls.sheet_names:
			bp = self.xls.parse(name)

			bp['sum'] = bp[self.basepair].sum(axis=1)
			bp = bp[ bp['sum'] >= self.constraints ]

			_, bp['minor'] = self.get_major_minor(bp[self.basepair])

			argx = -bp[self.basepair]
			bps = argx.values.argsort(axis=1)
			bp['Mj_seq'] = DataFrame(bps[:,0], index=bp.index)
			bp['Mn_seq'] = DataFrame(bps[:,1], index=bp.index)
			bp['maf'] = divide(bp[['minor']], bp[['sum']]) * 100

			bp['minor_'+name] = bp['minor']
			bp['Mj_seq_'+name] = bp['Mj_seq']
			bp['Mn_seq_'+name] = bp['Mn_seq']
			bp['maf_'+name] = bp['maf']
			self.minpos = min(self.minpos, bp[self.duma_position].values[0]) # if self.minpos else bp[self.duma_position].values[0]
			self.maxpos = max(self.maxpos, bp[self.duma_position].values[-1]) # if self.maxpos else bp[self.duma_position].values[-1]

			self.sheet.append(bp[ self.dumas_col + ['maf', 'minor', 'Mj_seq', 'Mn_seq', 'maf_'+name, 'minor_'+name, 'Mj_seq_'+name, 'Mn_seq_'+name] ])
			# minor 컬럼 추가..
			# 180725 : GenomeSt, RepeatRe, ORF 의 hspace를 위해 dumas_col 추가
			# 변경 전 : return bp[[duma_position] + ['maf', 'minor', 'Mj_seq', 'Mn_seq']]
			# return bp[ self.dumas_col + ['maf', 'minor', 'Mj_seq', 'Mn_seq'] ]
		
	def extract_dumas_info_(self):
		"""
			Get dumas info list... in duma_ORF, duma_Genome, duma_Repeat
		"""
		from pandas import unique

		assert len(self.sheet) > 0, "sheet list is empty"

		for i in unique( self.sheet[0][self.duma_Genome]):
			if str(i) != 'nan': 
				self.dumas_info[ self.duma_Genome ].append(i)

		for i in unique( self.sheet[0][self.duma_Repeat]):
			if str(i) != 'nan': 
				self.dumas_info[ self.duma_Repeat ].append(i)

		for i in unique( self.sheet[0][self.duma_ORF]):
			# ORF doesn't have '/'...
			if str(i) != 'nan' and str(i).find('/') == -1: 
				self.dumas_info[ self.duma_ORF].append(i)
	
	def extract_dumas_range_(self):
		"""
			get range of duma_position of each dumas_info
		"""
		for c in self.range_dumas.keys():
			for d in self.dumas_info[c]:
				d_ = self.sheet[0][ self.sheet[0][c] == d][self.duma_position].values
				self.range_dumas[c][d] = {'min': d_[0], 'max' : d_[-1] }

	def extract_incdec(self):
		"""
			Get change of data over the passage...
		"""
		from pandas import merge
		strp = ['low', 'middle', 'high']
		self.datalabel = []
		self.bpdiff, self.bpdiff_trans = [], []
		for i in range(self.nsheet):
			for j in range(self.nsheet):
				if i >= j: continue
				self.datalabel.append(strp[i]+'-'+strp[j])
				merged = merge(self.sheet[i], self.sheet[j], how='outer', on=self.dumas_col, left_index=True, right_index=True)

				# common
				cond = abs(merged['maf_x'] - merged['maf_y']) >= 5.0

				cond1 = merged['Mj_seq_x'] == merged['Mj_seq_y']
				cond2 = merged['Mj_seq_x'] != merged['Mj_seq_y']
				self.bpdiff.append(merged[ logical_and(cond, cond1).values ])

				# Actually, it doesn't show this data in current version.
				self.bpdiff_trans.append( merged[ logical_and(cond, cond2).values ])

	def get_bpdiff_cond(self, index, cond):
		return self.bpdiff[index][cond.values]

	def get_bpdiff(self, index, major, minor, base_=None, variation_=None):
		
		cond1 = self.bpdiff[index]['Mj_seq_x'] == major if major else True

		cond2 = self.bpdiff[index]['Mn_seq_y'] == minor if minor else True

		cond3 = self.bpdiff[index]['maf_x'] >= base_ if base_ else True

		cond4 = abs(self.bpdiff[index]['maf_y'] - self.bpdiff[index]['maf_x']) >= variation_ if variation_ else True

		return self.get_bpdiff_cond(index, logical_and( logical_and(cond1, cond2), logical_and(cond3, cond4) ) )

	# def get_nsheet(self):
	# 	return self.nsheet

	# def get_datalabel(self):
	# 	return self.datalabel

	# def get_bpdiff(self):
	# 	return self.bpdiff

	def get_dumas_info(self):
		return self.dumas_info

	def get_dumas_range(self):
		return self.range_dumas

	def get_bpdiff_cols(self, index, columns):
		return self.bpdiff[index][columns]

	# filtering된 major...처리 필요
	def get_parallel_data(self, bpdiff_, xmin, xmax, base_, variation_,  major, minor=None):
		def lambda_(x):
			if x['minor_x'] != 0 and x['Mn_seq_x'] != x['Mn_seq_y']:
				return self.basepair[int(x['Mn_seq_x'])].lower() + self.basepair[int(x['Mn_seq_y'])].lower()
			else:
				return self.basepair[int(x['Mn_seq_y'])].lower()
		def lambda_list(x):
			return [ x['maf_x'], x['maf_y'] ]
		def lambda_style(x):
			from mcolors import mcols
			if len(x['ct']) == 1:
				return self.colors[self.basepair.index(x['ct'].upper())]
			else:
				return mcols.cmap(self.dcolor[self.basepair.index(x['ct'][0].upper())][self.basepair.index(x['ct'][1].upper())])

		tmp = bpdiff_[ bpdiff_['Mj_seq_x'] == major ]
		data, label, style = [], [], []

		# for i in range(1, self.nsheet):
		cond0 = logical_and(tmp['maf_x']>= base_, abs(tmp['maf_y'] - tmp['maf_x']) >= variation_)
		cond1 = logical_and(tmp[self.duma_position] >= int(xmin), tmp[self.duma_position] <= int(xmax))
		cond2 = tmp['Mn_seq_y'] == minor if minor!=None else True
		tmp = tmp[ logical_and(cond0, logical_and(cond2,cond1)).values ]

		if len(tmp) > 0:
			tmp['ct'] = tmp.apply(lambda x : lambda_(x), axis=1)

			data = tmp.apply(lambda x: lambda_list(x), axis=1).values.tolist()
			label = tmp['ct'].values.tolist()
			style = tmp[['ct']].apply(lambda x: lambda_style(x), axis=1).values.tolist()

		return data, label, style

	def get_parallel_data_all(self, bpdiff_, xmin, xmax, base_, variation_, major, minor):
		pass
