from matplotlib import colors
import matplotlib.pyplot as plt
from numpy import where, linspace
from mcolors import mcols
from Slider import ZoomSlider

# filename = 'N_pOB strain_기본.xlsx'
filename = 'N_백신주 통합.xlsx'
basepair = ['A', 'G', 'C', 'T']
duma_position, duma_Genome, duma_Repeat, duma_ORF = 'Du_position', 'Genome strucure', 'Repeat region', 'ORF'
dumas_col = [ duma_position, duma_Genome, duma_Repeat, duma_ORF ]
dumas = { duma_Genome : [], duma_Repeat: [], duma_ORF : []}
range_dumas = { duma_Genome : dict(), duma_Repeat : dict(), duma_ORF : dict()}
min_, max_ = 0,0
seq_size, markersize = 100, 10
x = linspace(0.0, 1.0, seq_size)

mcolor = [ 'Blues', 'Orgs', 'Grns', 'Reds' ]
dcolor = [ [None, 'BuOg', 'BuGn', 'BuRd' ], 
		 ['OgBu', None, 'OgGn', 'OgRd'], 
		 ['GnBu', 'GnOg', None, 'GnRd'], 
		 ['RdBu', 'RdOg', 'RdGn', None] ]
colors = ['blue', 'orange', 'lightgreen', 'red']
labels = ['A', 'G', 'C', 'T']
annot, sct, leg = [], [], None

def update_xlim(val):
	if zslider.valleft > zslider.valright:
		return
	for x in ax:
		x.axes.set_xlim([zslider.valleft, zslider.valright])
	fig.canvas.draw_idle()	

# major 기준으로 filtering
def onpick(event):
	patch = event.artist
	vis = True if patch.get_alpha()==0.3 else False
	label = patch.get_label()
	color = 'grey'
	if label == 'A' or label == 'G' or label == 'C' or label == 'T':
		nmajor = basepair.index( label )
		color = colors[nmajor]
		for sc in sct:
			for mj, mn in zip(sc[nmajor]['major'], sc[nmajor]['minor']):
				mj.set_visible( vis )
				mn.set_visible( vis )
	else:
		global range_dumas
		for sc in range_dumas[label]['spans']:
			sc.set_visible(vis)
	if vis:
		patch.set_alpha(1.0)
	else:
		patch.set_facecolor(color)
		patch.set_alpha(0.3)
	fig.canvas.draw_idle()

def click_(event):
	if event.xdata is None:
		return
	pos = int(event.xdata)
	# print(pos)

	for i in range(len(bpdiff)):
		annot[i].set_visible(False)
		bp_ = bpdiff[i][ bpdiff[i][duma_position] == pos ]
		if len(bp_) > 0:
			# print(5)
			for sc_ in sct[i]:
				for s_ in sc_['major']:
					if s_ == None:
						continue
					# xy data -> xy (position in scatter)로 변환하는 방법
					if s_.get_offsets().__contains__(pos):
						ind = where(s_.get_offsets() == pos)
						print(ind)
						pos_ = s_.get_offsets()[ind[0]]						
						pos_[0][1] = bp_['maf_y']
						# print(pos_)
						annot[i].xy = pos_[0]
						text = "{}/{}\n{}".format(labels[bp_['Mj_seq_y'].values[0]], labels[bp_['Mn_seq_y'].values[0]].lower(), round(float(pos_[0][1]), 3))
						annot[i].set_text(text)
						annot[i].get_bbox_patch().set_facecolor(colors[labels.index(s_.get_label())])
						annot[i].get_bbox_patch().set_alpha(0.4)
						annot[i].set_visible(True)
	fig.canvas.draw_idle()

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

	global min_, max_
	min_ = min(min_, bp[duma_position].values[0]) if min_ else bp[duma_position].values[0]
	max_ = max(max_, bp[duma_position].values[-1]) if max_ else bp[duma_position].values[-1]

	# minor 컬럼 추가..
	# 180725 : GenomeSt, RepeatRe, ORF 의 hspace를 위해 dumas_col 추가
	# 변경 전 : return bp[[duma_position] + ['maf', 'minor', 'Mj_seq', 'Mn_seq']]
	return bp[ dumas_col + ['maf', 'minor', 'Mj_seq', 'Mn_seq'] ]

def get_uppery(y):
	y_ = [ i[-1] for i in y ]
	return y_

def get_dumas(sheet):
	from pandas import unique
	dumas = { duma_Genome : [], duma_Repeat: [], duma_ORF : []}
	for i in unique( sheet[duma_Genome] ):
		if str(i) != 'nan':
			dumas[duma_Genome].append(i)

	for i in unique( sheet[duma_Repeat] ):
		if str(i) != 'nan':
			dumas[duma_Repeat].append(i)

	for i in unique( sheet[duma_ORF] ):
		if str(i) == 'nan':
			continue
		if str(i).find('/') == -1:
			dumas[duma_ORF].append(i)
	return dumas

def get_dumas_range(sheet):
	global dumas
	col = [ duma_Genome, duma_Repeat, duma_ORF]
	for c in col:
		range_dumas[c]['spans'] = []
		for d in dumas[c]:
			d_ = sheet[ sheet[c] == d][duma_position].values
			range_dumas[c][d] = {'min': d_[0], 'max' : d_[-1] }
	return range_dumas

def draw_vspans(ax_):
	global range_dumas
	from numpy import random
	random.seed(1234)
	ymax = 50
	for c in range_dumas.keys():
		if c == duma_Genome:
			ymax = 0.7
		elif c == duma_Repeat:
			ymax = 0.85
		else:		
			ymax = 1		
		for d in range_dumas[c].keys():
			if d=='spans':
				continue
			r, g, b = random.uniform(0, 1, 3)
			range_dumas[c]['spans'].append(ax_.axvspan(range_dumas[c][d]['min'], range_dumas[c][d]['max'], 0, ymax, color=(r,g,b), alpha=0.2))
def button_click(label):
	print(label)
	pass

def draw_minor_transition(ax, sc, d, a, b):
	# transition from i to b
	from numpy import logical_and
	for i in range(4):
		if i == b:
			continue
		d_ = d[ d['Mn_seq_x'] == i ]

		y_ = [ linspace(p_.maf_x, p_.maf_y, seq_size) for p_ in d_[['maf_x', 'maf_y']].itertuples() ]
		x_ = array(y_)
		for j in range(len(d_)):
			x_[j].fill( d_[duma_position].values[j] )
		c = [ x for _ in range(len(x_))]
		
		sc[a]['minor'].append(ax.scatter(x_, y_, c=c, cmap=mcols.cmap(dcolor[i][b]), s=markersize))

		# draw major scatter
		y_ = get_uppery(y_)
		x_ = [ i[0] for i in x_]
		sc[a]['major'].append(ax.scatter(x_, y_, color=colors[a], label=basepair[a], s=markersize-1))


if __name__ == '__main__':
	import matplotlib.patches as mpatches
	from matplotlib.widgets import CheckButtons
	from pandas import ExcelFile, DataFrame, merge
	from numpy import linspace, logical_and, logical_or, array, abs
	xls = ExcelFile(filename)
	sheet = [ preprocessing(xls.parse(name)) for name in xls.sheet_names ]

	dumas = get_dumas(sheet[0])	
	nsheet = len(sheet)
	range_dumas = get_dumas_range(sheet[0])

	fig, ax = plt.subplots(3, 1, sharex=True)
	fig.set_size_inches(12, 7, forward=True)
	plt.subplots_adjust(left=0.05, bottom=0.15, right=0.98, hspace=0.15)

	bpdiff = []
	bpdiff_trans = [] # major transition

	strp = ['low', 'middle', 'high']
	ylabel = []#['low-middle', 'low-high', 'middle-high']
	ax = ax[::-1]

	spatches = [ mpatches.Patch(color=colors[idx], label=basepair[idx]) for idx in range(4) ]

	leg = fig.legend(handles=spatches, labels=basepair, bbox_to_anchor=(0.53, .90, 0.45, .102), loc=3,
           ncol=4, mode="expand", borderaxespad=0.)

	for i in range(nsheet):
		for j in range(nsheet):
			if i >= j:
				continue
			ylabel.append(strp[i]+'-'+strp[j])
			merged = merge(sheet[i], sheet[j], how='outer', on=dumas_col, left_index=True, right_index=True)

			# 공통 
			cond = abs(merged['maf_x'] - merged['maf_y']) >= 5.0
			cond1 = merged['Mj_seq_x'] == merged['Mj_seq_y']
			cond2 = merged['Mj_seq_x'] != merged['Mj_seq_y']
			bpdiff.append(merged[ logical_and(cond, cond1).values ])
			bpdiff_trans.append( merged[ logical_and(cond, cond2).values ])

	# 임의지정을 통한 교환
	ylabel[1], ylabel[2] = ylabel[2], ylabel[1]
	bpdiff[1], bpdiff[2] = bpdiff[2], bpdiff[1]
	bpdiff_trans[1], bpdiff_trans[2] = bpdiff_trans[2], bpdiff_trans[1]

	# major가 같은 것 들만
	for i in range(len(bpdiff)):
		sc = [ {'major' : [], 'minor': [] }, {'major' : [], 'minor': [] }, 
		       {'major' : [], 'minor': [] }, {'major' : [], 'minor': [] } ]
		
		for a in range(4):

			cond1 = bpdiff[i]['Mj_seq_x'] == a

			# major a -> minor b
			for b in range(4):
				if a == b:
					# sc[a]['major'].append(None)
					# sc[a]['minor'].append(None)
					continue

				# 1. minor가 같은 것들 처리 : 같은 색상
				cond2 = bpdiff[i]['Mn_seq_y'] == b
				# x는 
				cond3 = logical_or(bpdiff[i]['minor_x'] == 0, bpdiff[i]['Mn_seq_x'] == b)				
				cond3 = logical_and(cond2, cond3)

				d = bpdiff[i][ logical_and(cond1, cond3).values ]
				if len(d) > 0:
					y_ = [ linspace(d_.maf_x, d_.maf_y, seq_size) for d_ in d[['maf_x', 'maf_y']].itertuples() ]					
					x_ = array(y_)
					for j in range(len(d)):
						x_[j].fill( d[duma_position].values[j] )
					c = [ x for _ in range(len(x_))]
					sc[a]['minor'].append(ax[i].scatter(x_, y_, c=c, cmap=mcols.cmap(mcolor[b]), s=markersize, linewidth=0.0))

					# draw major scatter
					y_ = get_uppery(y_)
					x_ = [ ii[0] for ii in x_]
					sc[a]['major'].append(ax[i].scatter(x_, y_, color=colors[a], label=basepair[a], s=markersize-1))
				
				# # 2. minor가 바뀌는 것들 처리 : colors diverging
				# 각 basetype -> b type, target이 b, 즉, b로 변한것
				cond2 = logical_and(bpdiff[i]['Mn_seq_x'] != b, bpdiff[i]['minor_x'] > 0)				
				cond3 = logical_and(cond2, bpdiff[i]['Mn_seq_y'] == b)
				d = bpdiff[i][ logical_and(cond1, cond3).values ]
				if len(d) > 0:
					draw_minor_transition(ax[i], sc, d, a, b)

			# end for b
		#end for a
		sct.append(sc)

		draw_vspans(ax[i])

		ax[i].set_ylim([0, 50])
		ax[i].set_ylabel('Variation of MAF(%)')
		ax[i].set_title(ylabel[i], loc='left', fontdict={'fontsize':13, 'verticalalignment':'top', 'color':'black', 'backgroundcolor':'#FEFEFE'})

	# ZoomSlider of Dumas position
	zaxes = plt.axes([0.08, 0.07, 0.90, 0.03], facecolor='lightgoldenrodyellow')
	zslider = ZoomSlider(zaxes, 'Dumas Position', valmin=min_, valmax=max_, valleft=min_, valright=max_, color='lightblue', valstep=1.0)
	zslider.on_changed(update_xlim)

	labels = [ key for key in dumas.keys()]
	spatches2 = [ mpatches.Patch(color='grey', label=labels[idx]) for idx in range(3) ]
	leg2 = fig.legend(handles=spatches2, labels=labels, bbox_to_anchor=(0.05, .90, 0.44, .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)

	# lef picker
	leg.get_frame().set_alpha(0.4)
	leg2.get_frame().set_alpha(0.4)
	for legpatch in leg.get_patches():
		legpatch.set_picker(7)
	for legpatch in leg2.get_patches():
		legpatch.set_picker(7)


	# legend pick event
	fig.canvas.mpl_connect('pick_event', onpick)

	from matplotlib.backend_bases import MouseEvent
	for lb, legpatch in zip(labels,leg2.get_patches()):
		legpatch.pick(MouseEvent(lb, fig.canvas, 0, 0))		

	for x in ax:
		annot.append(x.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->")))

	fig.canvas.mpl_connect("button_press_event", click_)
	plt.show()

	# 추가할 것
	# Repeat, Genomes, ORF text 태깅..
	# legend text 클릭하라는 안내글..