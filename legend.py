
class Legend:

	def __init__(self, fig, config):
		self.fig = fig

		self.config = config

		self.leg = self.fig.legend(handles=config['spatches'], labels=config['label'], bbox_to_anchor=config['bbox'], loc=config['loc'],
			ncol=config['ncol'], mode=config['mode'], borderaxespad=config['borderaxespad'])

		for legpatch in self.leg.get_patches():
			legpatch.set_picker(7)

		self.onpick = None

	def pick(self):
		from matplotlib.backend_bases import MouseEvent
		for lb, legpatch in zip(self.config['label'],self.leg.get_patches()):
			legpatch.pick(MouseEvent(lb, self.fig.canvas, 0, 0))

	def get_legend(self):
		return self.leg