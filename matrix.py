import types
import sys
import random

class MatrixError(Exception):
	pass

class MatrixArithmeticError(MatrixError):
	def __init__(self, a, b, operation):
		self.a = a
		self.b = b
		self.operation = operation

	def __str__(self):
		return "Cannot %s a %dx%d and a %dx%d matrix" % (self.operation, self.a.rows(), self.a.cols(), self.b.rows(), self.b.cols())

class MatrixAdditionError(MatrixArithmeticError):
	def __init__(self, a, b):
		MatrixArithmeticError.__init__(self, a, b, "add")

class SquareError(MatrixError):
	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return "%s only defined for square matricies." % self.msg

class Matrix:
	def __init__(self, *args):
		if len(args) == 2:
			if isinstance(args[0], types.IntType) and isinstance(args[1], types.IntType):
				self.zeros(args[0], args[1])
			else:
				raise TypeError("Only two integer arguments are accepted.")
		elif len(args) == 1:
			if isinstance(args[0], types.IntType):
				self.zeros(args[0], args[0])
			elif isinstance(args[0], types.ListType):
				self.matrix = args[0]
			else:
				raise TypeError("Only an integer or a list is accepted for one argument.")
		else:
			raise TypeError("Only 1 or 2 arguments are accepted (%d given).") % len(args)
				
	def __str__(self):
		s = ""
		for row in self.matrix:
			s += "%s\n" % row
		return s

	def __getitem__(self, (row, col)):
		return self.matrix[row][col]

	def __setitem__(self, (row, col), value):
		self.matrix[row][col] = value
		
	def __add__(self, other):
		if not isinstance(other, Matrix):
			raise TypeError("Cannot add a matrix to type %s" % type(other))
		if not (self.cols() == other.cols() and self.rows() == other.rows()):
			raise MatrixAdditionError(self, other)
		r = []
		for row in xrange(self.rows()):
			r.append([])
			for col in xrange(self.cols()):
				r[row].append(self[(row, col)] + other[(row, col)])
		return Matrix(r)

	def __sub__(self, other):
		return self + -other

	def __neg__(self):
		return -1 * self

	def __mul__(self, other):
		if self.isScalarElement(other):
			return self.scalarMultiply(other)
		if not isinstance(other, Matrix):
			raise TypeError("Cannot multiply matrix and type %s" % type(other))
		#if other.is_row_vector():
		#	raise MatrixMultiplication_Error(self, other)
		#return self.matrix_multiply(other)

	def __rmul__(self, other):
		if not self.isScalarElement(other):
			raise TypeError("Cannot right-multiply by %s" % type(other))
		return self.scalarMultiply(other)

	def scalarMultiply(self, scalar):
		rows = []
		for row in self.matrix:				
			rows.append(map(lambda x: x * scalar, row))
		return Matrix(rows)

	def isScalarElement(self, x):
		return isinstance(x, types.IntType) or isinstance(x, types.FloatType) or isinstance(x, types.ComplexType)

	def isRowVector(self):
		return self.rows() == 1 and self.cols() > 1

	def isColumnVector(self):
		return self.cols() == 1 and self.rows() > 1
		
	def rows(self):
		return len(self.matrix)

	def cols(self):
		return len(self.matrix[0])
		
	def isSquare(self):
		return self.rows() == self.cols()
		
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
				
	def determinant(self):
		n = self.rows()
		if (n > 2):
			sign, t, sum = 1, 0, 0
			while t < n:
				d = {}
				t1 = 1
				while t1 < n:
					m = 0
					d[t1] = []
					while m < n:
						if (m == t):
							u = 0
						else:
							d[t1].append(self.matrix[t1][m])
						m += 1
					t1 += 1
				l1 = Matrix([d[x] for x in d])
				sum += sign * (self.matrix[0][t]) * (l1.determinant())
				sign *= (-1)
				t += 1
			return sum
		else:
			return (self.matrix[0][0] * self.matrix[1][1] - self.matrix[0][1] * self.matrix[1][0])

	@classmethod
	def makeRandom(cls, m, n, low=0, high=10):
		rows = []
		for x in range(m):
			rows.append([random.randrange(low, high) for i in range(n)])
		return Matrix(rows)

	@classmethod
	def readConsole(cls):
		print 'Enter matrix row by row. Type "q" to quit.'
		rows = []
		while True:
			line = sys.stdin.readline().strip()
			if line=='q': break
			row = [int(number) for number in line.split()]
			rows.append(row)
		return Matrix(rows)

	@classmethod
	def readFile(cls, fname):
		rows = []
		for line in open(fname).readlines():
			row = [int(number) for number in line.split()]
			rows.append(row)
		return Matrix(rows)
		
	@classmethod
	def identity(cls, rank):
		matrix = Matrix(rank)
		for index in xrange(matrix.rows()):
			matrix[(index, index)] = 1
		return matrix


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
		n = self.matrix.rows()
		tmp = self.matrix.identity(n) - self.matrix
		return abs(tmp.determinant()) != 1

	def check(self):
		if self.isExpansive() and self.unitCondition() and self.isCompleteResidues():
			return true
		return false

mat = Matrix([[6,1,1],[4,-2,5],[2,8,7]])
id = Matrix.identity(3)
digitSet = {0,1,2,3,4,5,6,7,8,9}

numsys = NumberSystem(mat, digitSet)
print numsys.unitCondition()
