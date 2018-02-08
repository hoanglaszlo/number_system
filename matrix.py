import types
import sys

class Matrix_Error(Exception):
	pass
	
class Square_Error(Matrix_Error):
	def __init__(self, func):
		self.func = func

	def __str__(self):
		return "%s only defined for square matricies." % self.func

class Matrix:
	def __init__(self, *args):
		if not (len(args) == 1 or len(args) == 2):
			raise TypeError("Only 1 or 2 arguments are accepted (%d given)") % len(args)
		if len(args) == 2:
			row, col = args
			self.zeros(row, col)
		else:
			if isinstance(args[0], types.IntType):
				row, col = args[0], args[0]
				self.zeros(args[0], args[0])
			else:
				self.matrix = args[0]

	def __str__(self):
		s = ""
		for row in self.matrix:
			s += "%s\n" % row
		return s

	def __getitem__(self, (row, col)):
		return self.matrix[row][col]

	def __setitem__(self, (row, col), value):
		self.matrix[row][col] = value
		
	@classmethod
	def makeMatrix(cls, rows):
		m = len(rows)
		n = len(rows[0])

		if any([len(row) != n for row in rows[1:]]):
			raise MatrixError, "Not a valid matrix."
		mat = Matrix(rows)

		return mat

	@classmethod
	def readConsole(cls):
		print 'Enter matrix row by row. Type "q" to quit.'
		rows = []
		while True:
			line = sys.stdin.readline().strip()
			if line=='q': break
			row = [int(x) for x in line.split()]
			rows.append(row)
		return cls.makeMatrix(rows)

	@classmethod
	def readFile(cls, fname):
		rows = []
		for line in open(fname).readlines():
			row = [int(x) for x in line.split()]
			rows.append(row)
		return cls.makeMatrix(rows)

	def rows(self):
		return len(self.matrix)

	def cols(self):
		return len(self.matrix[0])
		
	def zeros(self, row, col):
		if not row > 0:
			raise ValueError("Invalid number of rows (given %d)" % row)
		if not col > 0:
			raise ValueError("Invalid number of columns (given %d)" % col)
		self.matrix = []
		for i in xrange(row):
			self.matrix.append([])
			for j in xrange(col):
				self.matrix[i].append(0)


class NumberSystem:
	def __init__(self, matrix, digitSet):
		self.matrix = matrix
		self.digitSet = digitSet

	def hashFunction(self, U, G):
		return sum((U[i] % G[i,i])*prod(G[j,j] for j in range(i)) for i in range(rank(U)))

	def isCongruent(self, elementOne, elementTwo):
		G, U, V = m.smith_form()
		return self.hashFunction(U*elementOne, G) == self.hashFunction(U*elementTwo, G)

	def isCompleteResidues(self):
		for i in self.digitSet:
			if any(self.isCongruent(i,j) for j in self.digitSet-{i}):
				return false
		return true

	def isExpansive(self):
		eigens = self.matrix.eigenvalues()
		return all((abs(i)>1) for i in eigens)

	def unitCondition(self):
		tmp = matrix.identity(rank(self.matrix)) - self.matrix
		return abs(tmp.determinant()) != 1

	def check(self):
		if self.isExpansive() and self.unitCondition() and self.isCompleteResidues():
			return true
		return false

mat = Matrix.readFile("input.txt")
print mat
digitSet = {0,1,2,3,4,5,6,7,8,9}

#numsys = NumberSystem(mat, digitSet)
#numsys.check()
