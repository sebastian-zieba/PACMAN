class AncillaryData:
	def __init__(self, params):
		self.path = params['path']
		self.output = params['direct_image_output']#flag specifying whether coordinates are output to a file
		self.diagnostics = params['direct_image_diagnostics']#makes diagnostic plot if true

		# selects the boundaries of the region where the 2d Gaussian is fit,
		# these were selected by eye
		self.rmin = params['di_rmin']
		self.rmax = params['di_rmax']
		self.cmin = params['di_cmin']
		self.cmax = params['di_cmax']
		self.FILTER = params['FILTER']