"""
	Created by ilseop, Lee on 2018. 8. 31..
"""
import matplotlib.pyplot as plt
from numpy import linspace
from mcolors import mcols

class ParallelCoordinates:
	# all filter는 하나의 parallel로... 만들기
	def __init__(self, xmin, xmax, base_, variation_, size=1):
		self.fig, self.ax = plt.subplots(4,4, sharex=True, sharey=True)
		self.fig.set_size_inches(10, 10, forward=True)
		self.basepair = ['A', 'G', 'C', 'T']

		self.minorcolor = { 'a' : 'Blues', 'g' : 'Orgs' , 'c' : 'Grns', 't' : 'Reds'}

		self.dcolor = [ [None, 'BuOg', 'BuGn', 'BuRd' ], 
				 ['OgBu', None, 'OgGn', 'OgRd'], 
				 ['GnBu', 'GnOg', None, 'GnRd'], 
				 ['RdBu', 'RdOg', 'RdGn', None] ]
		self.fig.suptitle("Variations of Major(left) / Minor(top)\nfiltered by DumasPosition {} ~ {}, base : {}, variation : {}".format(int(xmin), int(xmax), int(base_), int(variation_)),  fontsize=20)
		
		#self.bpdiff['category'] = self.bpdiff.apply(lambda x : self.category(x), axis=1)

		def ticklabel(x, pos):
			if x== 0 :
				return 'low'
			elif x==5:
				return 'middle'
			else:
				return 'high'
		def ticklabel_size1(x, pos):
			if x== 0 :
				return 'low'
			elif x==5:
				return 'high'
		from matplotlib.ticker import FuncFormatter
		
		formatter = FuncFormatter(ticklabel_size1) if size == 1 else FuncFormatter(ticklabel)
		
		for a in range(4):
			for b in range(4):

				# self.ax[a, b].xaxis.set_visible(False)
				if b==0:
					self.ax[a, b].set_ylabel(self.basepair[a], fontdict={"fontsize":20})
				if a==0:
					self.ax[a, b].set_title(self.basepair[b].lower(), fontdict={'fontsize':20})	
				
				if a==b:
					self.ax[a, b].set_yticks([])
					continue
				self.ax[a, b].set_xlim([0, 5])
				self.ax[a, b].xaxis.set_major_formatter(formatter)
				# xlim unvisible
				if size == 2:
					self.ax[a, b].set_xlim([0, 10])
					# self.ax[a,b].xaxis.set_major_locator(plt.MaxNLocator(3))
					self.ax[a, b].set_xticks([0,5,10])

					par = [ self.ax[a, b].twinx() for i in range(2) ]
					par[0].spines['right'].set_position(('axes', 0.5))
					for p in par:
						p.set_ylim([0, 50])
				else:
					self.ax[a, b].twinx().set_ylim([0, 50])
					self.ax[a, b].set_xticks([0,5])

	# to book.py
	def category(self, x):
		mj = self.basepair[ int(x['Mj_seq_'+self.book.xls.sheet_names[0]]) ]
		mn,ischange = [], False
		for i in range(len(self.book.xls.sheet_names)):
			mn.append(self.basepair[ int(x['Mn_seq_'+self.book.xls.sheet_names[i]])].lower())
	
			if i != 0 and mn[i] != mn[i-1] and x['minor_'+self.book.xls.sheet_names[i-1]] != 0:
				ischange = True
		if ischange:
			mn_ = ''.join( i for i in mn)
		else:
			mn_ = mn[0]
		return mj+'/'+mn_


	def parallel_coordinates_in_axes(self, ax, data, labels, style, major, size=1):
		# print(data, labels, style, major)
		seq_size = 300
		x = [0, 5] if size == 1 else [5, 10]

		for d, s, l in zip(data, style, labels):
			# print(major, "  :  ", d, l)
			if len(l) == 1:
				ax.plot(x, d, color=s, label=major+'/'+l)
			else:
				y = linspace(d[0], d[1], seq_size)
				ax.scatter( linspace(x[0], x[1], seq_size), y, c=linspace(0,1,seq_size), cmap=s, label=major+'/'+l, s=1.5)
			# ax.legend()
		# print()
		return ax

class ParallelCoordinates_all:
	def __init__(self, xmin, xmax, base_, variation_, size=2, plotsize=4):
		self.fig, self.ax = plt.subplots(1, plotsize, sharey=True)
		self.fig.set_size_inches(10, 10, forward=True)
		self.basepair = ['A', 'G', 'C', 'T']

		self.minorcolor = { 'a' : 'Blues', 'g' : 'Orgs' , 'c' : 'Grns', 't' : 'Reds'}

		self.dcolor = [ [None, 'BuOg', 'BuGn', 'BuRd' ], 
				 ['OgBu', None, 'OgGn', 'OgRd'], 
				 ['GnBu', 'GnOg', None, 'GnRd'], 
				 ['RdBu', 'RdOg', 'RdGn', None] ]
		self.fig.suptitle("filtered by DumasPosition {} ~ {}, base : {}, variation : {}".format(int(xmin), int(xmax), int(base_), int(variation_)),  fontsize=20)

		def ticklabel(x, pos):
			if x== 0 :
				return 'low'
			elif x==5:
				return 'middle'
			else:
				return 'high'

		from matplotlib.ticker import FuncFormatter
		formatter = FuncFormatter(ticklabel)

		if plotsize == 1:
			self.ax.xaxis.set_major_formatter(formatter)
			self.ax.set_xlim([0, 10])

			self.ax.set_xticks([0, 5, 10])

			par = [ self.ax.twinx() for i in range(2) ]
			par[0].spines['right'].set_position(('axes', 0.5))
			for p in par:
				p.set_ylim([0, 50])
		else:
			for a in range(4):
				self.ax[a].xaxis.set_major_formatter(formatter)
				self.ax[a].set_title(self.basepair[a], fontdict={"fontsize":20})

				self.ax[a].set_xlim([0, 10])

				self.ax[a].set_xticks([0, 5, 10])

				par = [ self.ax[a].twinx() for i in range(2) ]
				par[0].spines['right'].set_position(('axes', 0.5))
				for p in par:
					p.set_ylim([0, 50])


	def parallel_coordinates_in_axes(self, ax, data, labels, style, major, size=1):
		# print(data, labels, style, major)
		seq_size = 300
		x = [0, 5] if size == 1 else [5, 10]
		for d, s, l in zip(data, style, labels):
			# print(major, "  :  ", d, l)
			if len(l) == 1:
				ax.plot(x, d, color=s, label=major+'/'+l)
			else:
				y = linspace(d[0], d[1], seq_size)
				ax.scatter( linspace(x[0], x[1], seq_size), y, c=linspace(0,1,seq_size), cmap=s, label=major+'/'+l, s=1.5)
			ax.legend()
		# print()
		return ax

		# # print(data, labels, style, major)
		# seq_size = 300
		
		# for d, s, l in zip(data, style, labels):
		# 	print(major, "  :  ", d, l)
		# 	for dd, ss, ll in zip(d,s,l):
		# 		if len(ll) == 1:
		# 			ax.plot([0,5], dd, color=ss, label=major+'/'+ll)
		# 		else:
		# 			y = linspace(dd[0], dd[1], seq_size)
		# 			ax.scatter( linspace(5, 10, seq_size), y, c=linspace(0,1,seq_size), cmap=ss, label=major+'/'+ll, s=1.5)
		# 	# ax.legend()
		# print()
		# return ax		