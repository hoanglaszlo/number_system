import types
import sys
import random
import operator
import unittest
from sympy import *

class MatrixError(Exception):
	pass

class MatrixArithmeticError(MatrixError):
	def __init__(self, elementLeft, elementRight, operation):
		self.elementLeft = elementLeft
		self.elementRight = elementRight
		self.operation = operation

	def __str__(self):
		return "Cannot %s a %dx%d and a %dx%d matrix" % (self.operation, self.elementLeft.rows(), self.elementLeft.cols(), self.elementRight.rows(), self.elementRight.cols())

class MatrixAdditionError(MatrixArithmeticError):
	def __init__(self, elementLeft, elementRight):
		MatrixArithmeticError.__init__(self, elementLeft, elementRight, "add")

class MatrixMultiplicationError(MatrixArithmeticError):
	def __init__(self, elementLeft, elementRight):
		MatrixArithmeticError.__init__(self, elementLeft, elementRight, "multiply")

class SquareError(MatrixError):
	def __init__(self, function):
		self.function = function

	def __str__(self):
		return "The %s function is only defined for square matricies." % self.function

class DeterminantError(SquareError):
	def __init__(self):
		SquareError.__init__(self, "determinant")

class InverseError(SquareError):
	def __init__(self):
		SquareError.__init__(self, "inverse")

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
			raise TypeError("Cannot add a matrix and a %s" % type(other))
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
		if self.is_scalar_element(other):
			return self.scalar_multiply(other)
		if not isinstance(other, Matrix):
			raise TypeError("Cannot multiply matrix and type %s" % type(other))
		return self.matrix_multiply(other)

	def __rmul__(self, other):
		if not self.is_scalar_element(other):
			raise TypeError("Cannot right-multiply by %s" % type(other))
		return self.scalar_multiply(other)

	def __eq__(self, other):
		if not isinstance(other, Matrix):
			raise TypeError("Cannot equal a matrix and a %s" % type(other))
		return all(self.row(i) == other.row(i) for i in xrange(self.rows()))

	def scalar_multiply(self, scalar):
		rows = []
		for row in self.matrix:				
			rows.append(map(lambda x: x * scalar, row))
		return Matrix(rows)

	def matrix_multiply(self, other):
		r = []
		if not isinstance(other, Matrix):
			raise TypeError("Cannot multiply matrix and type %s" % type(other))
		if not self.cols() == other.rows():
			raise MatrixMultiplicationError(self, other)
		for row in xrange(self.rows()):
			r.append([])
			for col in xrange(other.cols()):
				r[row].append(self.vector_inner_product(self.row(row), other.col(col)))
		if len(r) == 1 and len(r[0]) == 1:
			return r[0][0]
		else:
			return Matrix(r)

	def vector_inner_product(self, a, b):
		if not isinstance(a, types.ListType):
			raise TypeError("Only two lists are accepted.")
		if not isinstance(b, types.ListType):
			raise TypeError("Only two lists are accepted.")
		return reduce(operator.add, map(operator.mul, a, b))

	def is_scalar_element(self, x):
		return isinstance(x, types.IntType) or isinstance(x, types.FloatType) or isinstance(x, types.ComplexType)

	def is_row_vector(self):
		return self.rows() == 1 and self.cols() > 1

	def is_column_vector(self):
		return self.cols() == 1 and self.rows() > 1
	
	def row(self, i):
		return self.matrix[i]
		
	def col(self, j):
		r = []
		for row in self.matrix:
			r.append(row[j])
		return r
		
	def rows(self):
		return len(self.matrix)

	def cols(self):
		return len(self.matrix[0])
		
	def is_square(self):
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
		
	def smith_form(self):
		m = Matrix(self.matrix)
		U, G, V = Solver(m).smith_form()
		return U, G, V

	@classmethod
	def make_random(cls, m, n, low=0, high=10):
		rows = []
		for x in range(m):
			rows.append([random.randrange(low, high) for i in range(n)])
		return Matrix(rows)

	@classmethod
	def read_console(cls):
		print 'Enter matrix row by row. Type "q" to quit.'
		rows = []
		while True:
			line = sys.stdin.readline().strip()
			if line=='q': break
			row = [int(number) for number in line.split()]
			rows.append(row)
		return Matrix(rows)

	@classmethod
	def read_file(cls, fname):
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
		if not matrix.is_square():
			raise DeterminantError()
		if matrix.determinant() == 0:
			raise InverseError()
		#self.lattice = lattice
		self.matrix = matrix
		self.digitSet = digitSet

	def hash_function(self, U, G):
		s = self.find_in_diagonal(G, 1)
		return sum((U[i] % G[i,i]) * prod(G[j,j] for j in range(s + 1, i)) for i in range(s + 1, rank(U)))

	def is_congruent(self, elementOne, elementTwo):
		U, G, V = self.matrix.smith_form()
		return self.hash_function(U*elementOne, G) == self.hash_function(U*elementTwo, G)

	def find_in_diagonal(self, G, number):
		for index in xrange(G.rows()):
			if G[(index, index)] != number:
				return index 
		return G.rows()

	def is_complete_residues_system(self):
		for i in self.digitSet:
			if any(self.is_congruent(i,j) for j in self.digitSet - {i}):
				return False
		return True

	def is_expansive(self):
		eigens = self.matrix.eigenvalues()
		return all((abs(i)>1) for i in eigens)

	def unit_condition(self):
		n = self.matrix.rows()
		tmp = self.matrix.identity(n) - self.matrix
		return abs(tmp.determinant()) != 1

	def check(self):
		if self.is_expansive() and self.is_complete_residues_system():
			if self.unit_condition():
				return True
			else:
				print "It is okay, but... unit_condition failed"
		return False

class Solver:
	def __init__(self, matrix):
		self.matrix = matrix
		
	def leftmult2(self, m, i0, i1, a, b, c, d):
		for j in range(self.matrix.cols()):
			x, y = m[(i0,j)], m[(i1,j)]
			m[(i0,j)] = a * x + b * y
			m[(i1,j)] = c * x + d * y
	 
	def rightmult2(self, m, j0, j1, a, b, c, d):
		for i in range(self.matrix.rows()):
			x, y = m[(i,j0)], m[(i,j1)]
			m[(i,j0)] = a * x + c * y
			m[(i,j1)] = b * x + d * y
	 
	def smith_form(self, domain=ZZ):
		
		s = Matrix.identity(self.matrix.rows())
		t = Matrix.identity(self.matrix.cols())
		last_j = -1
		for i in range(self.matrix.rows()):
			for j in range(last_j+1, self.matrix.cols()):
				if any(i != 0 for i in self.matrix.col(j)):
					break
			else:
				break
			if self.matrix[(i,j)] == 0:
				for ii in range(self.matrix.rows()):
					if self.matrix[ii][j] != 0:
						break
				self.leftmult2(self.matrix, i, ii, 0, 1, 1, 0)
				self.rightmult2(s, i, ii, 0, 1, 1, 0)
			self.rightmult2(self.matrix, j, i, 0, 1, 1, 0)
			self.leftmult2(t, j, i, 0, 1, 1, 0)
			j = i
			upd = True
			while upd:
				upd = False
				for ii in range(i+1, self.matrix.rows()):
					if self.matrix[(ii,j)] == 0:
						continue
					upd = True
					if domain.rem(self.matrix[ii, j], self.matrix[i, j]) != 0:
						coef1, coef2, g = domain.gcdex(self.matrix[i,j], self.matrix[ii, j])
						coef3 = domain.quo(self.matrix[ii, j], g)
						coef4 = domain.quo(self.matrix[i, j], g)
						self.leftmult2(self.matrix,i, ii, coef1, coef2, -coef3, coef4)
						self.rightmult2(s, i, ii, coef4, -coef2, coef3, coef1)
					coef5 = domain.quo(self.matrix[ii, j], self.matrix[i, j])
					self.leftmult2(self.matrix, i, ii, 1, 0, -coef5, 1)
					self.rightmult2(s, i, ii, 1, 0, coef5, 1)
				for jj in range(j+1, self.matrix.cols()):
					if self.matrix[i, jj] == 0:
						continue
					upd = True
					if domain.rem(self.matrix[i, jj], self.matrix[i, j]) != 0:
						coef1, coef2, g = domain.gcdex(self.matrix[i,j], self.matrix[i, jj])
						coef3 = domain.quo(self.matrix[i, jj], g)
						coef4 = domain.quo(self.matrix[i, j], g)
						self.rightmult2(self.matrix, j, jj, coef1, -coef3, coef2, coef4)
						self.leftmult2(t, j, jj, coef4, coef3, -coef2, coef1)
					coef5 = domain.quo(self.matrix[i, jj], self.matrix[i, j])
					self.rightmult2(self.matrix, j, jj, 1, -coef5, 0, 1)
					self.leftmult2(t, j, jj, 1, coef5, 0, 1)
			last_j = j
		for i1 in range(min(self.matrix.rows(), self.matrix.cols())):
			for i0 in reversed(range(i1)):
				coef1, coef2, g = domain.gcdex(self.matrix[i0, i0], self.matrix[i1,i1])
				if g == 0:
					continue
				coef3 = domain.quo(self.matrix[i1, i1], g)
				coef4 = domain.quo(self.matrix[i0, i0], g)
				self.leftmult2(self.matrix, i0, i1, 1, coef2, coef3, coef2*coef3-1)
				self.rightmult2(s, i0, i1, 1-coef2*coef3, coef2, coef3, -1)
				self.rightmult2(self.matrix, i0, i1, coef1, 1-coef1*coef4, 1, -coef4)
				self.leftmult2(t, i0, i1, coef4, 1-coef1*coef4, 1, -coef1)
		return (s, self.matrix, t)
		
class MatrixTests(unittest.TestCase):
	def test_add_1(self):
		m1 = Matrix([[1, 2, 3], [4, 5, 6]])
		m2 = Matrix([[7, 8, 9], [10, 11, 12]])		
		m3 = m1 + m2
		self.assertTrue(m3 == Matrix([[8, 10, 12], [14,16,18]]))
		
	def test_add_2(self):
		m1 = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
		m2 = Matrix([[7, 8, 9], [10, 11, 12], [13, 14, 15]])		
		m3 = m1 + m2
		self.assertTrue(m3 == Matrix([[8, 10, 12], [14,16,18], [20, 22, 24]]))
		
	def test_add_3(self):
		m1 = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
		m2 = Matrix.identity(3)		
		m3 = m1 + m2
		self.assertTrue(m3 == Matrix([[2, 2, 3], [4, 6, 6], [7, 8, 10]]))

	def test_sub(self):
		m1 = Matrix([[1, 2, 3], [4, 5, 6]])
		m2 = Matrix([[7, 8, 9], [10, 11, 12]])		
		m3 = m2 - m1
		self.assertTrue(m3 == Matrix([[6, 6, 6], [6, 6, 6]]))

	def test_mul(self):
		m1 = Matrix([[1, 2, 3], [4, 5, 6]])
		m2 = Matrix([[7, 8], [10, 11], [12, 13]])
		id = Matrix.identity(3)
		self.assertTrue(m1 * m2 == Matrix([[63, 69], [150, 165]]))
		self.assertTrue(m2 * m1 == Matrix([[39, 54, 69], [54, 75, 96], [64, 89, 114]]))
		self.assertTrue(m1 * id == m1)
		self.assertTrue(id * m2 == m2)

	def test_det(self):
		m1 = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
		m2 = Matrix([[-1, -1], [1, -1]])
		m3 = Matrix([[6,1,1],[4,-2,5],[2,8,7]])
		id = Matrix.identity(3)
		self.assertTrue(m1.determinant() == 0)
		self.assertTrue(m2.determinant() == 2)
		self.assertTrue(m3.determinant() == -306)
		self.assertTrue(id.determinant() == 1)

class NumberSystemTests(unittest.TestCase):
	def test_unit_condition(self):
		mat1 = Matrix([[6,1,1],[4,-2,5],[2,8,7]])
		digitSet1 = {0,1,2,3,4,5,6,7,8,9}

		mat2 = Matrix([[-1,-1],[1,-1]])
		digitSet2 = {(0,0),(1,0)}
		
		self.assertTrue(NumberSystem(mat1, digitSet1).unit_condition() == True)
		self.assertTrue(NumberSystem(mat2, digitSet2).unit_condition() == True)
		
	def test_find_in_diagonal(self):
		mat1 = Matrix([[6,1,1],[4,-2,5],[2,8,7]])
		mat2 = Matrix([[1,1,1],[4,-2,5],[2,8,7]])
		mat3 = Matrix([[1,1,1],[4,1,5],[2,8,7]])
		mat4 = Matrix([[1,1,1],[4,1,5],[2,8,1]])
		digitSet = {0,1,2,3,4,5,6,7,8,9}
		self.assertTrue(NumberSystem(mat1, digitSet).find_in_diagonal(mat1, 1) == 0)
		self.assertTrue(NumberSystem(mat1, digitSet).find_in_diagonal(mat2, 1) == 1)
		self.assertTrue(NumberSystem(mat1, digitSet).find_in_diagonal(mat3, 1) == 2)
		self.assertTrue(NumberSystem(mat1, digitSet).find_in_diagonal(mat4, 1) == 3)

test = False

if test:
	if __name__ == "__main__":
		unittest.main()
else:
	mat = Matrix([[2,4,4],[-6,6,12],[10,-4,-16]])
	digitSet = {0,1,2,3,4,5,6,7,8,9}
