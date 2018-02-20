import types
import sys
import random
import operator
import unittest
from sympy import *
from sympy import Matrix as Mat
import copy
import math

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

class FormError(MatrixError):
	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return "%s" % self.msg
		
		
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
		if  isinstance(other, Matrix):
			return self.matrix_multiply(other)
		elif isinstance(other, tuple):
			return self.tuple_multiply(other)
		else:
			raise TypeError("Cannot multiply matrix and type %s" % type(other))

	def __rmul__(self, other):
		if not self.is_scalar_element(other):
			raise TypeError("Cannot right-multiply by %s" % type(other))
		return self.scalar_multiply(other)

	def __eq__(self, other):
		if not isinstance(other, Matrix):
			raise TypeError("Cannot equal a matrix and a %s" % type(other))
		return all(self.row(i) == other.row(i) for i in xrange(self.rows()))

	def scalar_multiply(self, scalar):
		if not self.is_scalar_element(scalar):
			raise TypeError("Cannot multiply matrix and type %s" % type(other))
		rows = []
		for row in self.matrix:				
			rows.append(map(lambda x: x * scalar, row))
		return Matrix(rows)

	def matrix_multiply(self, other):
		if not isinstance(other, Matrix):
			raise TypeError("Cannot multiply matrix and type %s" % type(other))
		if not self.cols() == other.rows():
			raise MatrixMultiplicationError(self, other)
		rows = []
		for row in xrange(self.rows()):
			rows.append([])
			for col in xrange(other.cols()):
				rows[row].append(self.vector_inner_product(self.row(row), other.col(col)))
		if len(rows) == 1 and len(rows[0]) == 1:
			return rows[0][0]
		else:
			return Matrix(rows)

	def tuple_multiply(self, other):
		if not isinstance(other, tuple):
			raise TypeError("Cannot multiply matrix and type %s" % type(other))
		if not self.cols() == len(other):
			raise MatrixMultiplicationError(self, other)
		rows = []
		for row in xrange(self.rows()):
			rows.append([])
			rows[row].append(self.vector_inner_product(self.row(row), other))
		if len(rows) == 1 and len(rows[0]) == 1:
			return rows[0][0]
		else:
			return Matrix(rows)

	def vector_inner_product(self, elementLeft, elementRight):
		if not isinstance(elementLeft, types.ListType):
			raise TypeError("Only two lists are accepted.")
		if not isinstance(elementRight, types.ListType) and not isinstance(elementRight, types.TupleType):
			raise TypeError("Only two lists are accepted.")
		return reduce(operator.add, map(operator.mul, elementLeft, elementRight))

	def is_scalar_element(self, element):
		return isinstance(element, types.IntType) or isinstance(element, types.FloatType) or isinstance(element, types.ComplexType)

	def is_row_vector(self):
		return self.rows() == 1 and self.cols() > 1

	def is_column_vector(self):
		return self.cols() == 1 and self.rows() > 1
	
	def is_vector(self):
		return self.is_row_vector() or self.is_column_vector()
	
	def row(self, index):
		return self.matrix[index]
		
	def col(self, index):
		col = []
		for row in self.matrix:
			col.append(row[index])
		return col
		
	def rows(self):
		return len(self.matrix)

	def cols(self):
		return len(self.matrix[0])
		
	def make_list(self):
		return [number for sublist in self.matrix for number in sublist]
		
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
				
	def ones(self, row, col):
		if not row > 0:
			raise ValueError("Invalid number of rows (given %d)" % row)
		if not col > 0:
			raise ValueError("Invalid number of columns (given %d)" % col)
		self.matrix = []
		for i in xrange(row):
			self.matrix.append([])
			for j in xrange(col):
				self.matrix[i].append(1)
		
	def determinant(self):
		if not self.is_square():
			raise SquareError("determinant")
		if self.rows() == 1:
			return self.matrix[0][0]
		i = 0
		sum = 0
		for j in xrange(self.rows()):
			if self.matrix[i][j] == 0:
				continue
			value = (-1) ** (i + j) * self.matrix[i][j] * self._A_ij(i, j).determinant()
			sum += value
		return sum
	
	def inverse(self):
		if not self.is_square():
			raise SquareError("inverse")
		if self.rows() == 1:
			return Matrix([[operator.truediv(1,self.matrix[0][0])]])
		d = self.determinant()
		if abs(d) < 10**-4:
			raise Exception('Matrix is not invertible')
		return operator.truediv(1,d) * self.adjugate()
	
	def cut(self, left = 0, right = None, top = 0, bottom = None):
		if right is None:
			right = self.cols()
		if bottom is None:
			bottom = self.rows()
		if not (left >= 0 and left < self.cols()):
			raise ValueError("left out of bounds")
		if not (right > 0 and right <= self.cols()):
			raise ValueError("right out of bounds")
		if not (top >= 0 and top < self.rows()):
			raise ValueError("top out of bounds")
		if not (bottom > 0 and bottom <= self.rows()):
			raise ValueError("bottom out of bounds'")
		if not (left < right):
			raise ValueError("left must be smaller than right")
		if not (top < bottom): 
			raise ValueError("top must be smaller than bottom")
		width = right - left
		height = bottom - top
		flat_values = self.make_list()
		rows = []
		for row in xrange(height):
			newrow = []
			for col in xrange(width):
				value = flat_values[self.cols() * top + left + self.cols() * row + col]
				newrow.append(value)
			rows.append(newrow)
		return Matrix(rows)	

	def _A_ij(self, i, j):
		if not (i >= 0 and i < self.rows()):
			raise ValueError("i out of bounds")
		if not (j >= 0 and j < self.cols()):
			raise ValueError("j out of bounds")
		if i == 0:
			m1 = self.cut(top = 1)
		elif i == self.rows() - 1:
			m1 = self.cut(bottom = self.rows() - 1)
		else:
			tm1 = self.cut(bottom = i)
			tm2 = self.cut(top = i + 1)
			m1 = stackv(tm1, tm2)
		if j == 0:
			m2 = m1.cut(left = 1)
		elif j == m1.cols() - 1:
			m2 = m1.cut(right = m1.cols() - 1)
		else:
			tm1 = m1.cut(right = j)
			tm2 = m1.cut(left = j + 1)
			m2 = stackh(tm1, tm2)
		return m2		

	def adjugate(self):
		if not self.is_square():
			raise SquareError("adjugate")
		rows = []
		for i in xrange(self.rows()):
			row = []
			for j in xrange(self.rows()):
				value = (-1) ** (i + j) * self._A_ij(j, i).determinant()
				row.append(value)
			rows.append(row)
		return Matrix(rows)		
			
	def smith_form(self):
		m = copy.deepcopy(self)
		U, G, V = Solver(m).smith_form()
		return U, G, V

	def eigenvalues(self):
		return Mat(self.matrix).eigenvals()
		
	def norm(self, type=2):
		if type == 1:
			return self._norm1()
		elif type == 2:
			return self._norm2()
		elif type == 'inf':
			return self._norm_inf()
		elif type == 'fro':
			return self._norm_fro()
		else:
			raise Exception('Illegal norm type')

	def _norm1(self):
		max = -1
		for j in xrange(self.cols()):
			value = sum(tuple(map(abs, self.col(j))))
			if value > max:
				max = value
		return max

	def _norm2(self):
		if not (self.is_row_vector() or self.is_column_vector()):
			#return math.sqrt(max((self.transpose() * self).eigenvalues()))
			raise FormError("Form not accepted.")
		elif self.is_row_vector():
			return math.sqrt(sum(tuple(map(lambda x: abs(x ** 2), self.row(0)))))
		elif self.is_column_vector():
			return math.sqrt(sum(tuple(map(lambda x: abs(x ** 2), self.col(0)))))

	def _norm_inf(self):
		max = -1
		for i in xrange(self.rows()):
			value = sum(tuple(map(abs, self.row(i))))
			if value > max:
				max = value
		return max

	def _norm_fro(self):
		sum = 0
		for i in xrange(self.rows()):
			for j in xrange(self.cols()):
				value = self.matrix[i][j]
				sum += abs(value ** 2)
		return math.sqrt(sum)
		
	def transpose(self):
		return Matrix([self.col(i) for i in xrange(self.cols())])
		
	def power(self, c):
		m = copy.deepcopy(self)
		for i in xrange(1,c):
			m = m * self
		return m
		
	@classmethod
	def from_tuple(cls, tup):
		if not isinstance(tup, tuple):
			raise TypeError("Cannot convert into matrix from type %s" % type(tup))
		rows = []
		for element in tup:
			rows.append([element])
		return Matrix(rows)
		
	@classmethod
	def make_random(cls, m, n, low=0, high=10):
		rows = []
		for x in xrange(m):
			rows.append([random.randrange(low, high) for i in xrange(n)])
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
	def __init__(self, matrix, digitSet, lattice={}):
		if not matrix.is_square():
			raise SquareError("determinant")
		if matrix.determinant() == 0:
			raise SquareError("inverse")
		self.lattice = lattice
		self.matrix = matrix
		self.digitSet = digitSet

	def find_in_diagonal(self, G, number):
		for index in xrange(G.rows()):
			if G[(index, index)] != number:
				return index 
		return G.rows()

	def hash_function(self, U, G):
		s = self.find_in_diagonal(G, 1)
		return sum((U[(i,0)] % G[(i,i)]) * prod(G[(j,j)] for j in range(s, i)) for i in range(s, U.rows()))

	def is_congruent(self, elementOne, elementTwo):
		U, G, V = self.matrix.smith_form()
		if not (abs(U.determinant()) == 1 and abs(V.determinant()) == 1):
			raise FormError("Smith normal form error")
		return self.hash_function(U * elementOne, G) == self.hash_function(U * elementTwo, G)

	def find_congruent(self, element):
		for i in self.digitSet:
			if self.is_congruent(i, element):
				return i
				
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
		
	def phi(self, element, n = 1, save = False):
		digSet = []
		M_inv = self.matrix.inverse()
		for i in range(n):
			d = self.find_congruent(element)
			digSet.append(d)
			if self.matrix.is_scalar_element(element):
				k = M_inv * (element - d)
				element = int(k[(0,0)])
			elif isinstance(element, Matrix):
				k = M_inv * (element - Matrix.from_tuple(d))
				element = k
			elif isinstance(element, tuple):
				k = M_inv * tuple(map(lambda x, y: x - y, element, d))
				element = k
		if save:
			return digSet
		else:
			return element
		
	def find_c(self):
		c = 1
		while not self.matrix.power(c).inverse().norm("inf") < 1:
			c += 1
		return c
		
	def find_gamma(self):
		c = self.find_c()
		norm = self.matrix.power(c).inverse().norm("inf")
		gamma = operator.truediv(1,1-norm)
		return gamma
	
class Solver:
	def __init__(self, matrix):
		self.matrix =  matrix
		
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
	def setUp(self):
		self.v1 = Matrix([[1, 2, 3]])
		self.v2 = Matrix([[4, 5, 6]])
		self.v3 = Matrix([[1], [2], [3]])
		
		self.c1 = Matrix([[10]])
		self.c2 = Matrix([[0]])
		
		self.m1 = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
		self.m2 = Matrix([[4, 1, -7, 2], [-1, 9, 6, 3]])
		self.m3 = Matrix([[8, -3, 1], [4, -6, 2], [7, 3, 5], [-2, -5, 1]])
		self.m4 = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
		self.m5 = Matrix([[6, 1, 1], [4,-2, 5], [2, 8, 7]])
		self.m6 = Matrix([[1, 3, 3], [1, 4, 3], [1, 3, 4]])
		self.m7 = Matrix([[1, 2, 3], [0, 1, 4], [5, 6, 0]])
		
	def testAdd(self):
		self.assertEqual(self.v1 + self.v2, Matrix([[5, 7, 9]]))
		self.assertEqual(self.m1 + self.m4, Matrix([[2, 2, 3], [4, 6, 6], [7, 8, 10]]))
		self.assertEqual(self.m1 + self.m5, Matrix([[7, 1, 1], [4, -1, 5], [2, 8, 8]]))
		self.assertEqual(self.m4 + self.m5, Matrix([[7, 3, 4], [8, 3, 11], [9, 16, 16]]))
		self.assertEqual(self.c1 + self.c2, Matrix([[10]]))
		with self.assertRaises(MatrixAdditionError):
			self.m2 + self.m3
		with self.assertRaises(MatrixAdditionError):
			self.v1 + self.v3

	def testSub(self):
		self.assertEqual(self.v2 - self.v1, Matrix([[3, 3, 3]]))
		self.assertEqual(self.m4 - self.m1, Matrix([[0, 2, 3], [4, 4, 6], [7, 8, 8]]))
		self.assertEqual(self.m5 - self.m1, Matrix([[5, 1, 1], [4,-3, 5], [2, 8, 6]]))
		self.assertEqual(self.m5 - self.m4, Matrix([[5,-1,-2], [0,-7,-1], [-5, 0,-2]]))
		with self.assertRaises(MatrixAdditionError):
			self.m2 - self.m3
		with self.assertRaises(MatrixAdditionError):
			self.v1 - self.v3
		
	def testMul(self):
		self.assertEqual(self.v1 * self.v3, 14)
		self.assertEqual(self.v2 * self.v3, 32)
		self.assertEqual(self.m1 * self.m4, Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
		self.assertEqual(self.m1 * self.m5, Matrix([[6, 1, 1], [4,-2, 5], [2, 8, 7]]))
		self.assertEqual(self.m4 * self.m5, Matrix([[20, 21, 32], [56, 42, 71], [92, 63, 110]]))
		self.assertEqual(self.m3 * self.m1, Matrix([[8, -3, 1], [4, -6, 2], [7, 3, 5], [-2, -5, 1]]))
		self.assertEqual(self.m3 * self.m5, Matrix([[38, 22, 0], [4, 32,-12], [64, 41, 57], [-30, 16, -20]]))		
		self.assertEqual(self.m4 * self.m1, Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
		self.assertEqual(self.m5 * self.m1, Matrix([[6, 1, 1], [4,-2, 5], [2, 8, 7]]))
		with self.assertRaises(MatrixMultiplicationError):
			self.m1 * self.m3
		with self.assertRaises(MatrixMultiplicationError):
			self.v1 * self.v2

	def testDet(self):
		self.assertEqual(self.m1.determinant(), 1)
		self.assertEqual(self.m4.determinant(), 0)
		self.assertEqual(self.m5.determinant(), -306)
		self.assertEqual(self.c1.determinant(), 10)
		self.assertEqual(self.c2.determinant(), 0)
		with self.assertRaises(SquareError):
			self.m2.determinant()
		with self.assertRaises(SquareError):
			self.m3.determinant()

	def testInv(self):
		self.assertEqual(self.c1.inverse(), Matrix([[0.1]]))
		self.assertEqual(self.m6.inverse(), Matrix([[7, -3, -3], [-1, 1, 0], [-1, 0, 1]]))
		self.assertEqual(self.m7.inverse(), Matrix([[-24, 18, 5], [20, -15, -4], [-5, 4, 1]]))
		with self.assertRaises(SquareError):
			self.m2.inverse()
		with self.assertRaises(SquareError):
			self.m3.inverse()
		with self.assertRaises(Exception):
			self.m4.inverse()
		
	def testNorm(self):
		self.assertAlmostEqual(self.m1.norm(1), 1)
		#self.assertAlmostEqual(self.m1.norm(2), 1)
		self.assertAlmostEqual(self.m1.norm("inf"), 1)
		self.assertAlmostEqual(self.m1.norm("fro"), 1.732, 3)
		
		self.assertAlmostEqual(self.m2.norm(1), 13)
		#self.assertAlmostEqual(self.m2.norm(2), 11.8579)
		self.assertAlmostEqual(self.m2.norm("inf"), 19)
		self.assertAlmostEqual(self.m2.norm("fro"), 14.0357, 3)
		
		self.assertAlmostEqual(self.m3.norm(1), 21)
		#self.assertAlmostEqual(self.m3.norm(2), 12.5191)
		self.assertAlmostEqual(self.m3.norm("inf"), 15)
		self.assertAlmostEqual(self.m3.norm("fro"), 15.5885, 3)
		
		self.assertAlmostEqual(self.m5.norm(1), 13)
		#self.assertAlmostEqual(self.m5.norm(2), 11.7476)
		self.assertAlmostEqual(self.m5.norm("inf"), 17)
		self.assertAlmostEqual(self.m5.norm("fro"), 14.1421, 3)
			
class NumberSystemTests(unittest.TestCase):
	def setUp(self):
		self.mat1 = Matrix([[6,1,1],[4,-2,5],[2,8,7]])
		self.mat2 = Matrix([[1,1,1],[4,-2,5],[2,8,7]])
		self.mat3 = Matrix([[1,1,1],[4,1,5],[2,8,7]])
		self.mat4 = Matrix([[1,1,1],[4,1,5],[2,8,1]])
		
		self.numsys1 = NumberSystem(Matrix([[10]]), {0,1,2,3,4,5,6,7,8,9})
		self.numsys2 = NumberSystem(Matrix([[-1,-1],[1,-1]]), {(0,0),(1,0)})
		self.numsys3 = NumberSystem(Matrix([[10]]), {0,1,2,3,4,5,6,7,8,18})
		self.numsys4 = NumberSystem(Matrix([[10]]), {0,1,2,3,4,5,6,7,8,39})
		self.numsys5 = NumberSystem(Matrix([[-1,-1],[1,-1]]), {(1,1),(1,0)})
		
	def testFindInDiagonal(self):
		self.assertEqual(self.numsys1.find_in_diagonal(self.mat1, 1), 0)
		self.assertEqual(self.numsys1.find_in_diagonal(self.mat2, 1), 1)
		self.assertEqual(self.numsys1.find_in_diagonal(self.mat3, 1), 2)
		self.assertEqual(self.numsys1.find_in_diagonal(self.mat4, 1), 3)

	def testIsCongruent(self):		
		self.assertTrue(self.numsys1.is_congruent(10,0))
		self.assertTrue(self.numsys1.is_congruent(14,4))
		self.assertTrue(self.numsys1.is_congruent(6,36))
		self.assertFalse(self.numsys1.is_congruent(13,12))
		self.assertFalse(self.numsys1.is_congruent(66,43))
		
		self.assertFalse(self.numsys2.is_congruent((0,0),(1,0)))

	def testFindCongruent(self):
		self.assertEqual(self.numsys1.find_congruent(13), 3)
		self.assertEqual(self.numsys1.find_congruent(15), 5)
		self.assertEqual(self.numsys1.find_congruent(64), 4)
		self.assertEqual(self.numsys1.find_congruent(8486), 6)
		
		self.assertEqual(self.numsys2.find_congruent((0,1)), (1,0))
		self.assertEqual(self.numsys2.find_congruent((1,1)), (0,0))
	
	def testIsCompleteResiduesSystem(self):
		self.assertTrue(self.numsys1.is_complete_residues_system())
		self.assertTrue(self.numsys2.is_complete_residues_system())
		self.assertFalse(self.numsys3.is_complete_residues_system())
		self.assertTrue(self.numsys4.is_complete_residues_system())
		self.assertTrue(self.numsys5.is_complete_residues_system())
	
	def testUnitCondition(self):
		self.assertTrue(self.numsys1.unit_condition())
		self.assertTrue(self.numsys2.unit_condition())
		
	def testIsExpansive(self):
		self.assertTrue(self.numsys1.is_expansive())
		self.assertTrue(self.numsys2.is_expansive())
	
	def testPhi(self):
		self.assertEqual(self.numsys1.phi(123456789, 0), 123456789)
		self.assertEqual(self.numsys1.phi(123456789, 1), 12345678)
		self.assertEqual(self.numsys1.phi(123456789, 2), 1234567)
		self.assertEqual(self.numsys1.phi(123456789, 5), 1234)
		self.assertEqual(self.numsys1.phi(123456789, 8), 1)
		self.assertEqual(self.numsys1.phi(123456789, 9), 0)
	
	def testCheck(self):
		self.assertTrue(self.numsys1.check())
		self.assertTrue(self.numsys2.check())
		self.assertFalse(self.numsys3.check())
		self.assertTrue(self.numsys4.check())
		self.assertTrue(self.numsys5.check())
		
		
def stackh(*matrices):
	matrices = _normalize_args(matrices)
	assert len(matrices) > 0, 'Can\'t stack zero matrices'
	for matrix in matrices:
		assert isinstance(matrix, Matrix), 'Can only stack matrices'
	height = matrices[0].rows()
	for matrix in matrices:
		assert matrix.rows() == height, 'Can\'t horizontally stack matrices with different heights'
	values = []
	for row in range(0, height):
		newrow = []
		for matrix in matrices:
			newrow += matrix.row(row)
		values.append(newrow)
	return Matrix(values)

def stackv(*matrices):
	matrices = _normalize_args(matrices)
	assert len(matrices) > 0, 'Can\'t stack zero matrices'
	for matrix in matrices:
		assert isinstance(matrix, Matrix), 'Can only stack matrices'
	width = matrices[0].cols()
	for matrix in matrices:
		assert matrix.cols() == width, 'Can\'t vertically stack matrices with different widths'
	values = []
	for matrix in matrices:
		values += matrix.matrix
	return Matrix(values)		

def _normalize_args(matrices):
	if len(matrices) > 0:
		first_elem = matrices[0]
		if isinstance(first_elem, list) or isinstance(first_elem, tuple):
			assert len(matrices) == 1, 'Couldn\'t normalize arguments'
			return first_elem
		return matrices
	return matrices

test = False

if __name__ == "__main__":
	if test:
		unittest.main()
	else:
		print "Hooray"
		numsys1 = NumberSystem(Matrix([[10]]), {0,1,2,3,4,5,6,7,8,9})
		numsys2 = NumberSystem(Matrix([[-1,-1],[1,-1]]), {(0,0),(1,0)})
		print numsys2.phi((4,5), )
