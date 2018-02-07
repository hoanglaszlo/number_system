import types

class Matrix:
	null_element = 0
	
	def __init__(self, *args):
		if not (len(args)a == 1 or len(args) == 2):
			raise TypeError("Only 1 or 2 arguments are accepted (%d given)") % len(args)
		if len(args) == 2:
			row, col = args
			self.zeros(row, col)
		else:
			row, col = args[0], args[0]
			self.zeros(args[0], args[0])

	def zeros(self, row, col):
		if not row > 0:
			raise ValueError("Invalid number of rows (given %d)" % row)
		if not col > 0:
			raise ValueError("Invalid number of columns (given %d)" % col)
		self.m = []
		for i in xrange(row):
			self.m.append([])
			for j in xrange(col):
				self.m[i].append(self.null_element)

	def __str__(self):
		s = ""
		for row in self.m:
			s += "%s\n" % row
		return s
		
m = Matrix(3)
print m
