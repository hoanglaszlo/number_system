import types
import sys
import random
import operator
import unittest
from sympy import Matrix as Mat
import copy
import math
import cmath

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
		if self.cols() > 1:
			for row in self.matrix:
				s += "%s\n" % row
		elif self.cols() == 1:
			for row in self.matrix:
				s += "%s, " % row[0]
			s = "(" + s[:-2] + ")"
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
		rows = []
		for row in xrange(self.rows()):
			rows.append([])
			for col in xrange(self.cols()):
				rows[row].append(self[(row, col)] + other[(row, col)])
		return Matrix(rows)

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
		if isinstance(other, int):
			return self[(0,0)] == other
		if not isinstance(other, Matrix):
			raise TypeError("Cannot equal a matrix and a %s" % type(other))
		return all(self.row(i) == other.row(i) for i in xrange(self.rows()))

	def __ne__(self, other):
		return not self.__eq__(other)
		
	def __pow__(self, n):
		result = self.identity(self.rows())
		for i in xrange(abs(n)):
			result *= self
		if  n < 0:
			return result.inverse()
		return result
		
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
		return Matrix(rows)

	def vector_inner_product(self, elementLeft, elementRight):
		if not isinstance(elementLeft, types.ListType):
			raise TypeError("Only two lists are accepted.")
		if not isinstance(elementRight, types.ListType) and not isinstance(elementRight, types.TupleType):
			raise TypeError("Only two lists are accepted.")
		return reduce(operator.add, map(operator.mul, elementLeft, elementRight))
	
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
		return [element for row in self.matrix for element in row]
		
	def is_square(self):
		return self.rows() == self.cols()
		
	def is_scalar_element(self, element):
		return isinstance(element, types.IntType) or isinstance(element, types.FloatType) or isinstance(element, types.ComplexType)

	def is_row_vector(self):
		return self.rows() == 1 and self.cols() > 1

	def is_column_vector(self):
		return self.cols() == 1 and self.rows() > 1
	
	def is_vector(self):
		return self.is_row_vector() or self.is_column_vector()
				
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
		
	def to_tuple(self):
		if not self.cols() == 1:
			raise Form("Only n:1 matrix can be converted.")
		list = []
		for row in self.matrix:
			list.append(row[0])
		return tuple(list)
		
	def to_int(self, func=int):
		rows = []
		for row in self.matrix:
			rows.append([int(func(element)) for element in row])
		return Matrix(rows)
		
	def zeros(self, row, col, init = 0):
		if not row > 0:
			raise ValueError("Invalid number of rows (given %d)" % row)
		if not col > 0:
			raise ValueError("Invalid number of columns (given %d)" % col)
		self.matrix = []
		for i in xrange(row):
			self.matrix.append([])
			for j in xrange(col):
				self.matrix[i].append(init)	
		
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
	def read_tuple(cls, tup):
		if not isinstance(tup, tuple):
			raise TypeError("Cannot convert into matrix from type %s" % type(tup))
		rows = []
		for element in tup:
			rows.append([element])
		return Matrix(rows)

	@classmethod
	def identity(cls, rank):
		matrix = Matrix(rank)
		for index in xrange(matrix.rows()):
			matrix[(index, index)] = 1
		return matrix

class NumberSystem:
	def __init__(self, matrix, digitSet, lattice={}):
		if isinstance(matrix, int):
			matrix = Matrix([[matrix]])
		if isinstance(matrix, Matrix):
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
		for digit in self.digitSet:
			if self.is_congruent(digit, element):
				return digit
				
	def is_complete_residues_system(self):
		for digit in self.digitSet:
			if any(self.is_congruent(digit, other) for other in self.digitSet - {digit}):
				return False
		return True

	def is_expansive(self):
		eigens = self.matrix.eigenvalues()
		return all((abs(value)>1) for value in eigens)

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
			digSet.append(element)
			if self.matrix.is_scalar_element(element):
				k = M_inv * (element - d)
				element = int(k[(0,0)])
			elif isinstance(element, tuple):
				k = M_inv * tuple(map(lambda x, y: x - y, element, d))
				element = k.to_int().to_tuple()
		if save:
			return digSet
		else:
			return element
		
	def find_c(self):
		c = 1
		while not Algorithm(self.matrix ** -c).norm("inf") < 0.001:
			c += 1
		return c
		
	def find_gamma(self):
		c = self.find_c()
		norm = (self.matrix**-c).norm("inf")
		gamma = 1 / (1 - norm)
		return gamma
	
	def calculate_xi_product(self, digit):
		c = self.find_c()
		matrix = Matrix(self.matrix.rows(),1)
		for j in xrange(1,c+1):
			matrix = stackh(matrix, (self.matrix**-j) * digit)
		matrix = matrix.cut(left=1)
		return matrix
		
	def calculate_from_matrices(self, func):
		c = self.find_c()
		xi = Matrix(self.matrix.rows(), c)	
		matrices = []
		for digit in self.digitSet:
			matrices.append(self.calculate_xi_product(digit))
	
		for i in xrange(xi.rows()):
			for j in xrange(xi.cols()):
				xi[(i,j)] = func(matrix[(i,j)] for matrix in matrices)
		return xi
		
	def calculate_box(self, type = "round"):
		c = self.find_c()
		gamma = self.find_gamma()
		eta = self.calculate_from_matrices(min).transpose()
		xi = self.calculate_from_matrices(max).transpose()
		sum_min = sum((Matrix.read_tuple(tuple(col)) for col in eta.matrix), Matrix(self.matrix.rows(), 1))
		sum_max = sum((Matrix.read_tuple(tuple(col)) for col in xi.matrix) , Matrix(self.matrix.rows(), 1))
		if type == "original" or type == 1:
			l = gamma * sum_min
			u = gamma * sum_max
		if type == "round" or type == 2:
			l = ((gamma * sum_min).to_int(math.ceil))
			u = ((gamma * sum_max).to_int(math.floor))
		if type == "periodic" or type == 3:
			l = - ((gamma * sum_min).to_int(math.ceil))
			u = - ((gamma * sum_max).to_int(math.floor))
		return l, u
		
	def calculate_box_phi(self):
		l, u = self.calculate_box()
		print l, u
		if l.rows() == 1:
			for i in range(l[(0,0)], u[(0,0)]+1):
				print i
		if l.rows() == 2:
			n=5
			graph = Matrix(n*2+1,n*2+1)
			for i in range(l[(0,0)], u[(0,0)]+1):
				for j in range(l[(1,0)], u[(1,0)]+1):
					print (i,j),  "->",  self.phi((i,j))
					coord = tuple(map(lambda x: x + n, self.phi((i,j))))
					#if i == u[(1,0)] or i == l[(1,0)]:
					#	graph[(i+n,j+n)] = 1
					#if j == l[(0,0)] or j == u[(0,0)]:
					#	graph[(i+n,j+n)] = 1
					graph[coord] = 1
		return graph

	def classify(self):
		finished = set()
		P = set()
		K = []
		l, u = self.calculate_box()
		print l, u
		for i in range(u[(0,0)], l[(0,0)]+1):
			for j in range(u[(1,0)], l[(1,0)]+1):
				K.append((i,j))
		print K
		for z in K:
			if z not in finished:
				orbit = set()
				while True:
					orbit.add(z)
					finished.add(z)
					z = self.phi(z)
					if z in finished or z not in K:
						break
				if z in orbit:
					P.update(self.get_cycle(z))
		return P
		
	def get_cycle(self, element):
		orbit = {element}
		while self.phi(element) not in orbit:
			orbit.add(self.phi(element))
		return orbit

class Algorithm:
	def __init__(self, matrix):
		self.matrix = matrix
		
	def transpose(self):
		return Matrix([self.matrix.col(i) for i in xrange(self.matrix.cols())])
		
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
		for j in xrange(self.matrix.cols()):
			value = sum(tuple(map(abs, self.matrix.col(j))))
			if value > max:
				max = value
		return max

	def _norm2(self):
		if not (self.matrix.is_row_vector() or self.matrix.is_column_vector()):
			#return math.sqrt(max((self.transpose() * self).eigenvalues()))
			raise FormError("Form not accepted.")
		elif self.matrix.is_row_vector():
			return math.sqrt(sum(tuple(map(lambda x: abs(x ** 2), self.matrix.row(0)))))
		elif self.matrix.is_column_vector():
			return math.sqrt(sum(tuple(map(lambda x: abs(x ** 2), self.matrix.col(0)))))

	def _norm_inf(self):
		max = -1
		for i in xrange(self.matrix.rows()):
			value = sum(tuple(map(abs, self.matrix.row(i))))
			if value > max:
				max = value
		return max

	def _norm_fro(self):
		sum = 0
		for i in xrange(self.matrix.rows()):
			for j in xrange(self.matrix.cols()):
				value = self.matrix[(i,j)]
				sum += abs(value ** 2)
		return math.sqrt(sum)	
		
	def cut(self, left = 0, right = None, top = 0, bottom = None):
		if right is None:
			right = self.matrix.cols()
		if bottom is None:
			bottom = self.matrix.rows()
		if not (left >= 0 and left < self.matrix.cols()):
			raise ValueError("left out of bounds")
		if not (right > 0 and right <= self.matrix.cols()):
			raise ValueError("right out of bounds")
		if not (top >= 0 and top < self.matrix.rows()):
			raise ValueError("top out of bounds")
		if not (bottom > 0 and bottom <= self.matrix.rows()):
			raise ValueError("bottom out of bounds'")
		if not (left < right):
			raise ValueError("left must be smaller than right")
		if not (top < bottom): 
			raise ValueError("top must be smaller than bottom")
		width = right - left
		height = bottom - top
		flat_values = self.matrix.make_list()
		rows = []
		for row in xrange(height):
			newrow = []
			for col in xrange(width):
				value = flat_values[self.matrix.cols() * top + left + self.matrix.cols() * row + col]
				newrow.append(value)
			rows.append(newrow)
		return Matrix(rows)	
		
class Solver:
	def __init__(self, matrix):
		self.matrix = matrix
		self.S = None;
		self.J = None;
		self.T = None;

	def cSwap(self,i,j):
		for k in range(self.J.rows()):
			temp = self.J[(k, i)]
			self.J[(k,i)] = self.J[(k,j)]
			self.J[(k,j)] = temp
		adjustment = Matrix.identity(self.T.rows())
		adjustment[(i,i)] = 0
		adjustment[(j,j)] = 0
		adjustment[(i,j)] = 1
		adjustment[(j,i)] = 1
		self.T = self.T*adjustment

	def cLC(self,I,i,j,a,b,gcd=None):
		if gcd is None or a==1 or a==-1:
			c = 0
			d = 1
		else:
			c = -self.J[(I,j)]/gcd
			d = self.J[(I,i)]/gcd
		temp = []
		for k in range(self.J.rows()):
			temp = self.J[(k,i)]
			self.J[(k,i)] = a*self.J[(k,i)] + b*self.J[(k,j)]
			self.J[(k,j)] = c*temp + d*self.J[(k,j)]
		adjustment = Matrix.identity(self.T.rows())
		adjustment[(i,i)] = a
		if i!=j:
			adjustment[(j,i)] = b 
			adjustment[(i,j)] = c
			adjustment[(j,j)] = d
		assert adjustment.determinant()== 1 or adjustment.determinant()== -1
		self.T = self.T*adjustment

	def rSwap(self,i,j):
		for k in range(self.J.cols()):
			temp = self.J[(i, k)]
			self.J[(i,k)] = self.J[(j,k)]
			self.J[(j,k)] = temp
		adjustment = Matrix.identity(self.S.rows())
		adjustment[(i,j)] = 1
		adjustment[(j,i)] = 1
		adjustment[(i,i)] = 0
		adjustment[(j,j)] = 0
		self.S = adjustment*self.S

	def rLC(self,I,i,j,a,b,gcd=None):
		if (gcd is None or a==1 or a==-1):
			c=0
			d=1
		else:
			c = -self.J[(j,I)]/gcd
			d = self.J[(i,I)]/gcd
		for k in range(self.J.cols()):
			temp = self.J[(i,k)]
			self.J[(i,k)] = a*self.J[(i,k)] + b*self.J[(j,k)]
			self.J[(j,k)] = c*temp + d*self.J[(j,k)]
		adjustment = Matrix.identity(self.S.rows())
		adjustment[(i,i)] = a
		if i!=j:
			adjustment[(i,j)] = b 
			adjustment[(j,i)] = c
			adjustment[(j,j)] = d
		assert adjustment.determinant() == 1 or adjustment.determinant() == -1
		self.S = adjustment*self.S

	def euclid(self,a, b):
		x0 = 1
		x1 = 0
		y0 = 0
		y1 = 1
		while b != 0:
			tempa = a
			tempb = b
			q = tempa / tempb
			a = tempb
			b = tempa % tempb
			tempx0 = x0
			x0 = x1
			x1 = tempx0 - q * x0
			tempy0 = y0
			y0 = y1
			y1 = tempy0 - q * y0
		return [a, x0, y0]

	def smith_form(self):
		self.J = copy.deepcopy(self.matrix)
		self.S = Matrix.identity(self.matrix.rows())
		self.T = Matrix.identity(self.matrix.cols())

		for i in range(min(self.J.rows(),self.J.rows())):
			if self.J[(i,i)] == 0: 
				foundReplacement = False;
				for j in range(i, self.J.rows()):
					if foundReplacement:
						break
					for k in range(i,self.J.cols()):
						if (self.J[(j,k)] != 0):
							foundReplacement = True
							break
				if not foundReplacement:
					break
				else:
					self.rSwap(i,j)
					self.cSwap(i,k)
			gcd=self.J[(i,i)]
			doneIteration = False
			while not doneIteration:
				if (self.J[(i,i)]==1 or self.J[(i,i)]==-1):
					break
				doneIteration = True
				for j in range(i+1, self.J.rows()):
					gcd, x, y = self.euclid(self.J[(i,i)], self.J[(j,i)])
					if (gcd < 0):
						gcd = -gcd
						x = -x
						y = -y
					if gcd == self.J[(i,i)]:
						pass
					elif gcd == self.J[(j,i)]:
						self.rSwap(i,j)
						doneIteration=False
					elif gcd == -self.J[(j,i)]:
						self.rSwap(i,j)
						self.rLC(i,i,i,-1, 0,gcd)
						doneIteration=False
					elif gcd < J[(i,i)] or gcd < -self.J[(i,i)]:
						self.rLC(i, i, j, x, y, gcd)
						doneIteration=False
				for j in range(i+1, self.J.cols()):
					gcd, x, y = self.euclid(self.J[(i,i)], self.J[(i,j)])
					if (gcd < 0):
						gcd = -gcd
						x = -x
						y = -y
					if gcd == self.J[(i,i)]:
						pass
					elif gcd == self.J[(i,j)]:
						self.cSwap(i,j)
						doneIteration=False
					elif gcd == -self.J[(i,j)]: #TODO WORK THIS BLOCK
						self.cSwap(i,j)
						self.cLC(i,i,i,-1, 0,gcd)
						doneIteration=False
					elif gcd < self.J[(i,i)] or gcd < -self.J[(i,i)]:
						self.cLC(i, i, j, x, y, gcd)
						doneIteration=False
			doneZeroing = False
			while not doneZeroing:
				doneZeroing = True
				for j in range(i+1, self.J.rows()):
					if self.J[(j,i)] != 0:
						self.rLC(i,j,i,1,-self.J[(j,i)]/self.J[(i,i)])
						if self.J[(j,i)] != 0:
							doneZeroing = False

				for j in range(i+1, self.J.cols()):
					if self.J[(i,j)] != 0:
						self.cLC(i,j,i,1,-self.J[(i,j)]/self.J[(i,i)])
						if self.J[(i,j)] != 0:
							doneZeroing = False

		for i in range(min(self.J.cols(), self.J.rows())-1):
			if self.J[(i+1,i+1)] == 0:
				return self.S,self.J,self.T
			gcd, x, y = self.euclid(self.J[(i,i)], self.J[(i+1,i+1)])

			if (gcd == self.J[(i+1,i+1)]):
				self.cSwap(i, i+1)
				self.rSwap(i, i+1)
			elif (gcd != self.J[(i,i)]):
				self.rLC(i,i, i+1, 1, 1) 
				self.cLC(i, i, i+1, x, y, gcd)
				self.cLC(i, i+1, i, 1, -self.J[(i,i+1)]/self.J[(i,i)])
				self.rLC(i, i+1, i, 1, -self.J[(i+1,i)]/self.J[(i,i)])
		return self.S,self.J,self.T
				
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
		
		self.numsys1 = NumberSystem(10, {0,1,2,3,4,5,6,7,8,9})
		self.numsys2 = NumberSystem(Matrix([[-1,-1],[1,-1]]), {(0,0),(1,0)})
		self.numsys3 = NumberSystem(10, {0,1,2,3,4,5,6,7,8,18})
		self.numsys4 = NumberSystem(10, {0,1,2,3,4,5,6,7,8,39})
		self.numsys5 = NumberSystem(Matrix([[-1,-1],[1,-1]]), {(1,1),(1,0)})
		self.numsys6 = NumberSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(0,1),(0,-1),(-2,-3)})
		self.numsys7 = NumberSystem(3, {-2,0,2})
		self.numsys8 = NumberSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(2,0),(3,0),(4,0)})
			
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
		
	def testCalculatingBox(self):
		self.assertEqual(self.numsys2.calculate_box(), (Matrix([[0],[0]]),Matrix([[0],[0]])))
		self.assertEqual(self.numsys1.calculate_box(), (Matrix([[0]]),Matrix([[0]])))
		self.assertEqual(self.numsys7.calculate_box(), (Matrix([[-1]]),Matrix([[1]])))
		self.assertEqual(self.numsys6.calculate_box(), (Matrix([[-2],[-1]]),Matrix([[0],[0]])))
		self.assertEqual(self.numsys8.calculate_box(), (Matrix([[0],[-2]]),Matrix([[2],[0]])))
		
class AlgorithmTests(unittest.TestCase):
	def setUp(self):
		self.al = Algorithm(Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
		
	def testNorm(self):
		self.assertAlmostEqual(self.al.norm(1), 1)
		#self.assertAlmostEqual(self.m1.norm(2), 1)
		self.assertAlmostEqual(self.al.norm("inf"), 1)
		self.assertAlmostEqual(self.al.norm("fro"), 1.732, 3)
		
def stackh(*matrices):
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
	
def prod(iterable):
    return reduce(operator.mul, iterable, 1)

test = True

if __name__ == "__main__":
	if test:
		unittest.main()
	else:
		#print "Hooray"
		numsys0 = NumberSystem(Matrix([[10]]), {0,1,2,3,4,5,6,7,8,9}) 						#(0)-(0)
		numsys1 = NumberSystem(Matrix([[-1,-1],[1,-1]]), {(0,0),(1,0)}) 					#(0,0)-(0,0)
		numsys2 = NumberSystem(Matrix([[3]]), {-2,0,2})										#(-1)-(1)
		numsys3 = NumberSystem(Matrix([[0,-3],[1,0]]), {(0,0),(1,0),(-1,1)})				#
		numsys4a = NumberSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(2,0),(3,0),(4,0)})	#(0,-2)-(2,0)
		numsys4b1 = NumberSystem(Matrix([[3,0],[0,3]]), {(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)}) #(-0.5,-0.5)-(0.5,0.5)
		numsys4b2 = NumberSystem(Matrix([[3,0],[0,3]]), {(0,0),(1,0),(2,0),(0,1),(0,2),(1,2),(2,1),(-1,2),(-2,1)}) #(-1,0)-(1,1)
		numsys5 = NumberSystem(Matrix([[1,-2],[1,1]]), {(0,0),(1,0),(-1,0)})				#
		numsys7 = NumberSystem(Matrix([[-3,-1],[1,-3]]), {(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(6,0),(7,0),(8,0),(9,0)})				#?
		numsys8a = NumberSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(0,1),(0,-1),(-6,-5)})	#
		numsys8b = NumberSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(0,1),(0,-1),(-2,-3)})
		print numsys8a.classify()
