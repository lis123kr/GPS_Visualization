from numpy import tile
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Slider, CheckButtons

filename = 'N_백신주 통합.xlsx'

class Visualization:

	duma_col = ['Genome strucure', 'Repeat region', 'ORF', 'Du_position', 'Du_seq']

	basepair = ["A", "G", "C", "T"]

	colors = ['blue', 'orange', 'green', 'red']

	def __init__(self, title, sheet):
		from pandas import DataFrame
		from numpy import divide
		bp = sheet[ self.duma_col+['seq', 'pos']+self.basepair ]
		bp = bp[ bp['seq'] != '-' ]
		
		bp['sum'] = bp[self.basepair].sum(axis=1)
		bp = bp[ bp['sum'] >= 35.0 ]
		
		major, minor = self.get_major_minor(bp[self.basepair])
		bp['minor'] = minor

		argx = -bp[self.basepair]
		bps = argx.values.argsort(axis=1)
		bp['Mn_seq'] = DataFrame(bps[:,1], index=[bp.index])

		bp['maf'] = divide(bp[['minor']], bp[['sum']]) * 100

		self.genomes = bp['Genome strucure'].dropna().unique().tolist()
		self.repeats = bp['Repeat region'].dropna().unique().tolist()

		self.ge_visible = [True for _ in range(len(self.genomes))]
		self.re_visible = [True for _ in range(len(self.repeats))]
		self.bplist = [bp[ bp['Mn_seq'] == i ][self.duma_col + ['pos', 'maf']] for i in range(len(self.basepair))]

		self.title = title
		self.sct = []

		del bp, bps, major, minor, argx

	def make_plot(self):
		self.fig, self.ax = plt.figure(self.title), plt.axes()

		self.ax.set_title(self.title)
		self.fig.set_size_inches(10, 6, forward=True)
		self.fig.subplots_adjust(left=0.20, bottom=0.45, right=0.95)
		self.ax.set_xlabel("position")
		self.ax.set_ylabel('Minor(%)')
		self.ax.set_ylim([0.0, 50.0])		

		spatches = []

		for idx, t in enumerate(self.bplist):
			sc = self.ax.scatter(t['pos'], t['maf'], s=4, c=self.colors[idx], label=self.basepair[idx])
			self.sct.append(sc)
			spatches.append(mpatches.Patch(color=self.colors[idx], label=self.basepair[idx]))
			# del sc
		
		self.ax.legend(handles=spatches, bbox_to_anchor=(0.65, 1.01, 0.35, 0.03), loc=3,
	           ncol=4, mode="expand", borderaxespad=0.)

		marker = self.fig.add_axes([0.65, 0.95, 0.30, 0.03], facecolor='lightgoldenrodyellow')
		self.msize = Slider(marker, 'marker size', 0.0, 30.0, valinit=4/3, color='lightblue')

		mnmin = self.fig.add_axes([0.22, 0.3, 0.70, 0.03], facecolor='lightgoldenrodyellow')
		mnmax = self.fig.add_axes([0.22, 0.35, 0.70, 0.03], facecolor='lightgoldenrodyellow')

		self.smax = Slider(mnmax, 'y_max', 0.1, 50.0, valinit=50.0, color='lightblue')
		self.smin = Slider(mnmin, 'y_min', 0.0, 50.0, valinit=0.0, color='lightblue')

		self.smax.on_changed(self.update)
		self.smin.on_changed(self.update)
		self.msize.on_changed(self.markersize)

		# rectangles = [ mpatches.Rectangle((1,1),1,1, facecolor=colors[i]) for i in range(0,4) ] #dot.get_facecolor() for dot in sct ]
		self.labels = [ dot.get_label() for dot in self.sct]
		visibility = [dot.get_visible() for dot in self.sct]

		# Genome Structure
		self.text2 = self.fig.text(0.025, 0.955, 'Genome Structure', fontsize=9)
		cbax2 = self.fig.add_axes([0.025, 0.8, 0.1, 0.15], facecolor='#fafafa')
		self.cb2 = CheckButtons(cbax2, self.genomes, [True for _ in range(len(self.genomes))])
		self.cb2.on_clicked(self.ge_clicked)
		
		# seq
		self.text = self.fig.text(0.025, 0.705, 'base type', fontsize=9)
		cbax = self.fig.add_axes([0.025, 0.55, 0.1, 0.15], facecolor='#fafafa')
		self.cb = CheckButtons(cbax, self.labels, visibility)
		self.cb.on_clicked(self.seq_clicked)

		# Repeat Region
		self.text3 = self.fig.text(0.025, 0.455, 'Repeat Region', fontsize=9)
		cbax3 = self.fig.add_axes([0.025, 0.3, 0.1, 0.15], facecolor='#fafafa')
		self.cb3 = CheckButtons(cbax3, self.repeats, [True for _ in range(len(self.repeats))])
		self.cb3.on_clicked(self.re_clicked)

		for i, r in enumerate(self.cb.rectangles):
			r.set_facecolor(self.colors[i])
			self.cb.labels[i].set_color(self.colors[i])
			r.set_edgecolor('k')
			r.set_alpha(0.65)

		self.annot = self.ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
	                    bbox=dict(boxstyle="round", fc="w"),
	                    arrowprops=dict(arrowstyle="->"))

		self.fig.canvas.mpl_connect("motion_notify_event", self.hover)

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

	def update(self, val):
	    if self.smax.val <= self.smin.val:
	    	return
	    self.ax.axes.set_ylim([self.smin.val, self.smax.val])
	    # fig.canvas.draw_idle()

	def markersize(self, val):
		val = 3 * val
		for sc in self.sct:
			sc.set_sizes([val,])

	def seq_clicked(self, label):
		index = self.labels.index(label)
		self.sct[index].set_visible(not self.sct[index].get_visible())
		# visiable[index] = not visiable[index]
		plt.draw()

	def ge_clicked(self, label):
		index = self.genomes.index(label)
		for i in range(4):
			Npts = len(self.sct[i].get_offsets())
			fc = self.sct[i].get_facecolors()
			if len(fc) == 1:
				fc = tile(fc, Npts).reshape(Npts, -1)
			# bp를 A,G,C,T별로 분리
			ind = self.bplist[i]['Genome strucure'].isin([label])
			if len(fc[ind]) is not 0:
				fc[ind, -1] = 0.08 if self.ge_visible[index] else 1.0
				self.sct[i].set_facecolors(fc)
		self.ge_visible[index] = not self.ge_visible[index]
		del fc, ind
		plt.draw()

	def re_clicked(self, label):
		index = self.repeats.index(label)
		for i in range(4):
			Npts = len(self.sct[i].get_offsets())
			fc = self.sct[i].get_facecolors()
			if len(fc) == 1:
				fc = tile(fc, Npts).reshape(Npts, -1)
			ind = self.bplist[i]['Repeat region'].isin([label])
			if len(fc[ind]) is not 0:
				fc[ind, -1] = 0.08 if self.re_visible[index] else 1.0			
				self.sct[i].set_facecolors(fc)
		self.re_visible[index] = not self.re_visible[index]
		del fc, ind
		plt.draw()

	# instance of scatter, index
	def update_annot(self, sc, ind):
		pos = sc.get_offsets()[ind['ind'][0]]
		seq = sc.get_label()
		self.annot.xy = pos
		text = "{}, {}".format(seq, pos[0])
		self.annot.set_text(text)
		self.annot.get_bbox_patch().set_facecolor(self.colors[self.labels.index(seq)])
		self.annot.get_bbox_patch().set_alpha(0.4)

	def hover(self, event):
		vis = self.annot.get_visible()
		if event.inaxes == self.ax:
			for sc in self.sct:
				cont, ind = sc.contains(event)
				if cont:
					self.update_annot(sc, ind)
					self.annot.set_visible(True)
					self.fig.canvas.draw_idle()
				elif vis:
					self.annot.set_visible(False)
					self.fig.canvas.draw_idle()

if __name__ == '__main__':
	from pandas import ExcelFile

	xls = ExcelFile(filename)

	v = []

	for name in xls.sheet_names:
		sheet = xls.parse( name )
		v.append(Visualization(name, sheet))
	for s in v:
		s.make_plot()

	plt.show()



