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

class Encoder:
	def __init__(self, filename, inverse = False):
		self.sequence = self.read(filename)
		self.hash_table = self.create_hash_table()
		self.encoded = ""
		self.inverse = inverse
		self.code(inverse)
		
	def __str__(self):
		return self.encoded
		
	def read(self, filename):
		with open('Phyllorhiza_punctata.txt', 'r') as myfile:
			self.sequence = myfile.read().replace('\n', '')
		return self.sequence
		
	def create_hash_table(self):
		return {'a': '0', 'c': '1', 'g': '2', 't': '3', " " : ""}
		
	def code(self, inverse):
		if not inverse:
			for base in self.sequence:
				self.encoded += self.hash_table[base]
		else:
			for base in self.sequence:
				self.encoded = self.hash_table[base] + self.encoded
		
class Matrix:
	def __init__(self, *args):
		if len(args) == 2:
			if isinstance(args[0], int) and isinstance(args[1], int):
				self.zeros(args[0], args[1])
			else:
				raise TypeError("Only two integer arguments are accepted.")
		elif len(args) == 1:
			if isinstance(args[0], int):
				self.zeros(args[0], args[0])
			elif isinstance(args[0], list):
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
			return Algorithm(result).inverse()
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
		if not isinstance(elementLeft, list):
			raise TypeError("Only two lists are accepted.")
		if not isinstance(elementRight, list) and not isinstance(elementRight, tuple):
			raise TypeError("Only two lists are accepted.")
		return reduce(operator.add, map(operator.mul, elementLeft, elementRight))
	
	def row(self, index):
		return self.matrix[index]
		
	def row_add(self, indexOne, indexTwo, mul):
		value =  (Matrix([self.matrix[indexOne]]) + mul * Matrix([self.matrix[indexTwo]])).row(0)
		self.matrix[indexOne] = value
		
	def col(self, index):
		col = []
		for row in self.matrix:
			col.append(row[index])
		return col
		
	def col_add(self, indexOne, indexTwo, mul):
		value =  (Matrix([self.col(indexOne)]) + mul * Matrix([self.col(indexTwo)])).row(0)
		self.matrix[indexOne] = value
		
	def rows(self):
		return len(self.matrix)

	def cols(self):
		return len(self.matrix[0])
		
	def make_list(self):
		return [element for row in self.matrix for element in row]
		
	def is_square(self):
		return self.rows() == self.cols()
		
	def is_scalar_element(self, element):
		return isinstance(element, int) or isinstance(element, float) or isinstance(element, complex)

	def is_row_vector(self):
		return self.rows() == 1 and self.cols() > 1

	def is_column_vector(self):
		return self.cols() == 1 and self.rows() > 1
	
	def is_vector(self):
		return self.is_row_vector() or self.is_column_vector()
		
	def is_expansive(self):
		eigens = Algorithm(self).eigenvalues()
		return all((abs(value)>1) for value in eigens)
		
	def to_tuple(self):
		if not self.cols() == 1:
			raise Form("Only n:1 matrix can be converted.")
		list = []
		for row in self.matrix:
			list.append(row[0])
		return tuple(list)
		
	def to_int(self, func = int):
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
	def diagonal(self, diag):
		n = len(diag)
		rows = []
		for i in xrange(n):
			rows.append([])
			for j in xrange(n):
				if i == j:
					rows[i].append(diag[i])	
				else:
					rows[i].append(0)
		return Matrix(rows)
	
	@classmethod
	def companion(self, polynom):
		n = len(polynom) - 1
		rows = []
		for i in xrange(n):
			rows.append([])
			for j in xrange(n):
				if j == n - 1:
					rows[i].append(-polynom[i])	
				elif i == j + 1:
					rows[i].append(1)
				else:
					rows[i].append(0)
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

#class Vector(Matrix):
	
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
		
	def _A_ij(self, i, j):
		if not (i >= 0 and i < self.matrix.rows()):
			raise ValueError("i out of bounds")
		if not (j >= 0 and j < self.matrix.cols()):
			raise ValueError("j out of bounds")
		if i == 0:
			m1 = self.cut(top = 1)
		elif i == self.matrix.rows() - 1:
			m1 = self.cut(bottom = self.matrix.rows() - 1)
		else:
			tm1 = self.cut(bottom = i)
			tm2 = self.cut(top = i + 1)
			m1 = stackv(tm1, tm2)
		if j == 0:
			m2 = Algorithm(m1).cut(left = 1)
		elif j == m1.cols() - 1:
			m2 = Algorithm(m1).cut(right = m1.cols() - 1)
		else:
			tm1 = Algorithm(m1).cut(right = j)
			tm2 = Algorithm(m1).cut(left = j + 1)
			m2 = stackh(tm1, tm2)
		return m2
		
	def determinant(self):
		if not self.matrix.is_square():
			raise DeterminantError()
		if self.matrix.rows() == 1:
			return self.matrix[(0,0)]
		i = 0
		sum = 0
		for j in xrange(self.matrix.rows()):
			if self.matrix[(i,j)] == 0:
				continue
			value = (-1) ** (i + j) * self.matrix[(i,j)] * Algorithm(self._A_ij(i, j)).determinant()
			sum += value
		return sum
		
	def adjugate(self):
		if not self.matrix.is_square():
			raise DeterminantError()
		rows = []
		for i in xrange(self.matrix.rows()):
			row = []
			for j in xrange(self.matrix.rows()):
				value = (-1) ** (i + j) * Algorithm(self._A_ij(j, i)).determinant()
				row.append(value)
			rows.append(row)
		return Matrix(rows)	
		
	def inverse(self):
		if not self.matrix.is_square():
			raise InverseError()
		if self.matrix.rows() == 1:
			return Matrix([[operator.truediv(1,self.matrix[(0,0)])]])
		d = self.determinant()
		if abs(d) < 10**-4:
			raise Exception('Matrix is not invertible')
		return operator.truediv(1,d) * self.adjugate()
	
	def find_in_diagonal(self, number):
		for index in xrange(self.matrix.rows()):
			if self.matrix[(index, index)] != number:
				return index 
		return self.matrix.rows()

	def smith_form(self):
		m = copy.deepcopy(self.matrix)
		U, G, V = Solver(m).smith_form()
		return U, G, V

	def eigenvalues(self):
		return Mat(self.matrix.matrix).eigenvals()
		
class DigitSet:
	def __init__(self, matrix, digitSet = set(), jcan = 0, jsym = 0):
		self.matrix = matrix
		self.size = abs(Algorithm(self.matrix).determinant())
		self.set = set()
		if len(digitSet) != 0:
			self.set = digitSet
		elif jcan:
			self.get_j_canonical(jcan)
		elif jsym:
			self.get_j_symmetric(jsym)
			
	def __str__(self):
		return repr(self.set)
		
	def get_j_canonical(self, j):
		digit = [0] * self.matrix.rows()
		self.set.add(tuple(digit))
		for i in xrange(1, self.size):
			digit[j - 1] = i
			self.set.add(tuple(digit))

	def get_j_symmetric(self, j):
		digit = [0] * self.matrix.rows()
		self.set.add(tuple(digit))
		min = int(math.floor(-self.size/2))
		max = int(math.floor(self.size/2))
		for i in xrange(min, max):
			digit[j - 1] = i
			self.set.add(tuple(digit))
		
class RadixSystem:
	def __init__(self, matrix, digitSet = {}, lattice={}, jcan = 0, jsym = 0):
		if isinstance(matrix, int):
			matrix = Matrix([[matrix]])
		if isinstance(matrix, Matrix):
			if not matrix.is_square():
				raise DeterminantError()
			if Algorithm(matrix).determinant() == 0:
				raise InverseError()
		self.lattice = lattice
		self.matrix = matrix
		if len(digitSet) != 0:
			self.digitSet = digitSet
		elif jcan:
			self.digitSet = DigitSet(self.matrix, jcan = jcan).set
		elif jsym:
			self.digitSet = DigitSet(self.matrix, jsym = jsym).set
		self.smithForm = SmithForm(matrix)
		self.hashTable = self.create_hash_table()
		
	def hash_function(self, z):
		s = Algorithm(self.smithForm.G).find_in_diagonal(1)
		Uz = self.smithForm.U * z
		return sum((Uz[(i,0)] % self.smithForm.G[(i,i)]) * prod(self.smithForm.G[(j,j)] for j in range(s, i)) for i in range(s, Uz.rows()))
	
	def create_hash_table(self):
		hashTable = [None] * len(self.digitSet)
		for digit in self.digitSet:
			hashTable[self.hash_function(digit)] = digit;
		return hashTable

	def is_congruent(self, elementOne, elementTwo):
		return self.hash_function(elementOne) == self.hash_function(elementTwo)
		
	def find_congruent(self, element):
		return self.hashTable[self.hash_function(element)]
				
	def is_complete_residues_system(self):
		for digit in self.digitSet:
			if any(self.is_congruent(digit, other) for other in self.digitSet - {digit}):
				return False
		return True

	def necessary_condition(self):
		if self.matrix.is_expansive() and self.is_complete_residues_system():
			return True
		else:
			return False
		
class NumberSystem(RadixSystem):
	def __init__(self, matrix, digitSet = {}, lattice={}, jcan = 0, jsym = 0):
		RadixSystem.__init__(self, matrix, digitSet, lattice, jcan, jsym)
		self.M_inv = Algorithm(self.matrix).inverse()
		self.C_con = self.find_c()
		self.gamma = self.find_gamma()

	def __str__(self):
		print self.matrix
		print repr(self.digitSet)
		return ""
		
	def unit_condition(self):
		n = self.matrix.rows()
		tmp = self.matrix.identity(n) - self.matrix
		return abs(Algorithm(tmp).determinant()) != 1

	def necessary_condition(self):
		if self.matrix.is_expansive() and self.is_complete_residues_system() and self.unit_condition():
			return True
		else:
			return False
			
	def phi(self, element, n = 1, save = False):
		digSet = []
		for i in range(n):
			d = self.find_congruent(element)
			digSet.append(element)
			if self.matrix.is_scalar_element(element):
				k = self.M_inv * (element - d)
				element = int(k[(0,0)])
			elif isinstance(element, tuple):
				k = self.M_inv * tuple(map(lambda x, y: x - y, element, d))
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
		norm = Algorithm(self.matrix ** -self.C_con).norm("inf")
		gamma = 1 / (1 - norm)
		return gamma
	
	def calculate_fundamental_matrix(self, digit):
		matrix = Matrix(self.matrix.rows(),1)
		for j in xrange(1, self.C_con + 1):
			matrix = stackh(matrix, (self.matrix ** -j) * digit)
		matrix = Algorithm(matrix).cut(left=1)
		return matrix
		
	def calculate_fundamental_minmax(self, func):
		minmax = Matrix(self.matrix.rows(), self.C_con)	
		matrices = []
		for digit in self.digitSet:
			matrices.append(self.calculate_fundamental_matrix(digit))
		for i in xrange(minmax.rows()):
			for j in xrange(minmax.cols()):
				minmax[(i,j)] = func(matrix[(i,j)] for matrix in matrices)
		return minmax
		
	def calculate_fundamental_box_border(self, type = "periodic"):
		min_mat = Algorithm(self.calculate_fundamental_minmax(min)).transpose()
		max_mat = Algorithm(self.calculate_fundamental_minmax(max)).transpose()
		sum_min = sum((Matrix.read_tuple(tuple(col)) for col in min_mat.matrix), Matrix(self.matrix.rows(), 1))
		sum_max = sum((Matrix.read_tuple(tuple(col)) for col in max_mat.matrix), Matrix(self.matrix.rows(), 1))
		if type == "original" or type == 1:
			l = self.gamma * sum_min
			u = self.gamma * sum_max
		if type == "round" or type == 2:
			l = ((self.gamma * sum_min).to_int(math.ceil))
			u = ((self.gamma * sum_max).to_int(math.floor))
		if type == "periodic" or type == 3:
			l = - ((self.gamma * sum_min).to_int(math.ceil))
			u = - ((self.gamma * sum_max).to_int(math.floor))
		return l, u
		
	def classify(self):
		finished = []
		cycles = []
		points = []
		l, u = self.calculate_fundamental_box_border()
		for i in range(u[(0,0)], l[(0,0)]+1):
			for j in range(u[(1,0)], l[(1,0)]+1):
				points.append((i,j))
		for z in points:
			if z not in finished:
				orbit = []
				while z not in finished and z in points:
					orbit.append(z)
					finished.append(z)
					z = self.phi(z)
				if z in orbit:
					cycles.append(self.get_cycle(z))
		return cycles
		
	def get_cycle(self, element):
		cycle = []
		while element not in cycle:
			cycle.append(element)
			element = self.phi(element)
		return cycle
	
class SmithForm:
	def __init__(self, matrix):
		self.matrix = matrix
		self.G = copy.deepcopy(self.matrix)
		self.U = Matrix.identity(self.matrix.rows())
		self.V = Matrix.identity(self.matrix.cols())
		self.smith_form()

	def cSwap(self, i, j):
		for k in range(self.G.rows()):
			temp = self.G[(k, i)]
			self.G[(k,i)] = self.G[(k,j)]
			self.G[(k,j)] = temp
		adjustment = Matrix.identity(self.V.rows())
		adjustment[(i,i)] = 0
		adjustment[(j,j)] = 0
		adjustment[(i,j)] = 1
		adjustment[(j,i)] = 1
		self.V = self.V * adjustment

	def cLC(self, I, i, j, a, b, gcd = None):
		if gcd is None or abs(a) == 1:
			c = 0
			d = 1
		else:
			c = -self.G[(I,j)] / gcd
			d = self.G[(I,i)] / gcd
		temp = []
		for k in range(self.G.rows()):
			temp = self.G[(k,i)]
			self.G[(k,i)] = a * self.G[(k,i)] + b * self.G[(k,j)]
			self.G[(k,j)] = c * temp + d * self.G[(k,j)]
		adjustment = Matrix.identity(self.V.rows())
		adjustment[(i,i)] = a
		if i != j:
			adjustment[(j,i)] = b 
			adjustment[(i,j)] = c
			adjustment[(j,j)] = d
		assert Algorithm(adjustment).determinant() == 1 or Algorithm(adjustment).determinant() == -1
		self.V = self.V * adjustment

	def rSwap(self, i, j):
		for k in range(self.G.cols()):
			temp = self.G[(i, k)]
			self.G[(i,k)] = self.G[(j,k)]
			self.G[(j,k)] = temp
		adjustment = Matrix.identity(self.U.rows())
		adjustment[(i,j)] = 1
		adjustment[(j,i)] = 1
		adjustment[(i,i)] = 0
		adjustment[(j,j)] = 0
		self.U = adjustment * self.U

	def rLC(self, I, i, j, a, b, gcd = None):
		if (gcd is None or abs(a) == 1):
			c = 0
			d = 1
		else:
			c = -self.G[(j,I)] / gcd
			d = self.G[(i,I)] / gcd
		for k in range(self.G.cols()):
			temp = self.G[(i,k)]
			self.G[(i,k)] = a * self.G[(i,k)] + b * self.G[(j,k)]
			self.G[(j,k)] = c * temp + d * self.G[(j,k)]
		adjustment = Matrix.identity(self.U.rows())
		adjustment[(i,i)] = a
		if i != j:
			adjustment[(i,j)] = b 
			adjustment[(j,i)] = c
			adjustment[(j,j)] = d
		assert Algorithm(adjustment).determinant() == 1 or Algorithm(adjustment).determinant() == -1
		self.U = adjustment * self.U

	def euclid(self, a, b):
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
		for i in range(min(self.G.rows(),self.G.rows())):
			if self.G[(i,i)] == 0: 
				foundReplacement = False;
				for j in range(i, self.G.rows()):
					if foundReplacement:
						break
					for k in range(i,self.G.cols()):
						if (self.G[(j,k)] != 0):
							foundReplacement = True
							break
				if not foundReplacement:
					break
				else:
					self.rSwap(i,j)
					self.cSwap(i,k)
			gcd=self.G[(i,i)]
			doneIteration = False
			while not doneIteration:
				if (abs(self.G[(i,i)]) == 1):
					break
				doneIteration = True
				for j in range(i+1, self.G.rows()):
					gcd, x, y = self.euclid(self.G[(i,i)], self.G[(j,i)])
					if (gcd < 0):
						gcd = -gcd
						x = -x
						y = -y
					if gcd == self.G[(i,i)]:
						pass
					elif gcd == self.G[(j,i)]:
						self.rSwap(i,j)
						doneIteration = False
					elif gcd == -self.G[(j,i)]:
						self.rSwap(i,j)
						self.rLC(i, i, i, -1, 0, gcd)
						doneIteration = False
					elif gcd < self.G[(i,i)] or gcd < -self.G[(i,i)]:
						self.rLC(i, i, j, x, y, gcd)
						doneIteration = False
				for j in range(i+1, self.G.cols()):
					gcd, x, y = self.euclid(self.G[(i,i)], self.G[(i,j)])
					if (gcd < 0):
						gcd = -gcd
						x = -x
						y = -y
					if gcd == self.G[(i,i)]:
						pass
					elif gcd == self.G[(i,j)]:
						self.cSwap(i,j)
						doneIteration = False
					elif gcd == -self.G[(i,j)]: #TODO WORK THIS BLOCK
						self.cSwap(i,j)
						self.cLC(i, i, i, -1, 0, gcd)
						doneIteration = False
					elif gcd < self.G[(i,i)] or gcd < -self.G[(i,i)]:
						self.cLC(i, i, j, x, y, gcd)
						doneIteration = False
			doneZeroing = False
			while not doneZeroing:
				doneZeroing = True
				for j in range(i+1, self.G.rows()):
					if self.G[(j,i)] != 0:
						self.rLC(i, j, i, 1, -self.G[(j,i)] / self.G[(i,i)])
						if self.G[(j,i)] != 0:
							doneZeroing = False
				for j in range(i+1, self.G.cols()):
					if self.G[(i,j)] != 0:
						self.cLC(i, j, i, 1, -self.G[(i,j)] / self.G[(i,i)])
						if self.G[(i,j)] != 0:
							doneZeroing = False
		for i in range(min(self.G.cols(), self.G.rows())-1):
			if self.G[(i+1, i+1)] == 0:
				return self.U, self.G, self.V
			gcd, x, y = self.euclid(self.G[(i,i)], self.G[(i+1,i+1)])

			if (gcd == self.G[(i+1,i+1)]):
				self.cSwap(i, i+1)
				self.rSwap(i, i+1)
			elif (gcd != self.G[(i,i)]):
				self.rLC(i,i, i+1, 1, 1) 
				self.cLC(i, i, i+1, x, y, gcd)
				self.cLC(i, i+1, i, 1, -self.G[(i,i+1)] / self.G[(i,i)])
				self.rLC(i, i+1, i, 1, -self.G[(i+1,i)] / self.G[(i,i)])
		return self.U, self.G, self.V
"""			
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
"""
class MatrixTests(unittest.TestCase):
	def setUp(self):
		self.c0 = Matrix([[0]])
		self.c1 = Matrix([[1]])
		self.c2 = Matrix([[10]])
		self.c3 = Matrix([[-1]])

		self.v0 = Matrix([[0, 0]])
		self.v1 = Matrix([[1, 2]])
		self.v2 = Matrix([[1, 2, 3]])
		self.v3 = Matrix([[1, -2, 3, -4]])
		
		self.v4 = Matrix([[0], [0]])
		self.v5 = Matrix([[1], [2]])
		self.v6 = Matrix([[1], [2], [3]])
		self.v7 = Matrix([[1], [-2], [3], [-4]])
		
		self.m0 = Matrix([[0, 0], [0, 0]])
		self.m1 = Matrix([[1, 0], [0, 1]])
		self.m2 = Matrix([[1, 2], [2, 3]])
		self.m3 = Matrix([[4, -2], [-2, 3]])
		
		self.m4 = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
		self.m5 = Matrix([[6, 1, 1], [4,-2, 5], [2, 8, 7]])
		self.m6 = Matrix([[1, 3, 3], [1, 4, 3], [1, 3, 4]])
		self.m7 = Matrix([[1, 2, 3], [0, 1, 4], [5, 6, 0]])
		
		self.m8 = Matrix([[4, 1, -7, 2], [-1, 9, 6, 3]])
		self.m9 = Matrix([[8, -3, 1], [4, -6, 2], [7, 3, 5], [-2, -5, 1]])
		
	def testEq(self):
		self.assertEqual(self.c0, 0)
		self.assertEqual(self.c1, 1)
		self.assertEqual(self.c2, 10)
		self.assertEqual(self.c3, -1)
		self.assertEqual(self.c0, self.c0)
		self.assertEqual(self.c3, -1)
		self.assertEqual(self.c3, -1)
		self.assertEqual(self.c3, -1)
			
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
			
	def testUnitCondition(self):
		self.assertTrue(self.numsys1.unit_condition())
		self.assertTrue(self.numsys2.unit_condition())
	
	def testPhi(self):
		self.assertEqual(self.numsys1.phi(123456789, 0), 123456789)
		self.assertEqual(self.numsys1.phi(123456789, 1), 12345678)
		self.assertEqual(self.numsys1.phi(123456789, 2), 1234567)
		self.assertEqual(self.numsys1.phi(123456789, 5), 1234)
		self.assertEqual(self.numsys1.phi(123456789, 8), 1)
		self.assertEqual(self.numsys1.phi(123456789, 9), 0)
	
	def testCheck(self):
		self.assertTrue(self.numsys1.necessary_condition())
		self.assertTrue(self.numsys2.necessary_condition())
		self.assertFalse(self.numsys3.necessary_condition())
		self.assertTrue(self.numsys4.necessary_condition())
		self.assertTrue(self.numsys5.necessary_condition())
		
	def testCalculatingBox(self):
		self.assertEqual(self.numsys2.calculate_fundamental_box_border(2), (Matrix([[0],[0]]),Matrix([[0],[0]])))
		self.assertEqual(self.numsys1.calculate_fundamental_box_border(2), (Matrix([[0]]),Matrix([[0]])))
		self.assertEqual(self.numsys7.calculate_fundamental_box_border(2), (Matrix([[-1]]),Matrix([[1]])))
		self.assertEqual(self.numsys6.calculate_fundamental_box_border(2), (Matrix([[-2],[-1]]),Matrix([[0],[0]])))
		self.assertEqual(self.numsys8.calculate_fundamental_box_border(2), (Matrix([[0],[-2]]),Matrix([[2],[0]])))
		
class RadixSystemTests(unittest.TestCase):
	def setUp(self):
		self.numsys1 = RadixSystem(10, {0,1,2,3,4,5,6,7,8,9})
		self.numsys2 = RadixSystem(Matrix([[-1,-1],[1,-1]]), {(0,0),(1,0)})
		self.numsys3 = RadixSystem(10, {0,1,2,3,4,5,6,7,8,18})
		self.numsys4 = RadixSystem(10, {0,1,2,3,4,5,6,7,8,39})
		self.numsys5 = RadixSystem(Matrix([[-1,-1],[1,-1]]), {(1,1),(1,0)})
		self.numsys6 = RadixSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(0,1),(0,-1),(-2,-3)})
		self.numsys7 = RadixSystem(3, {-2,0,2})
		self.numsys8 = RadixSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(2,0),(3,0),(4,0)})

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

class AlgorithmTests(unittest.TestCase):
	def setUp(self):
		self.v1 = Algorithm(Matrix([[1, 2, 3]]))
		self.v2 = Algorithm(Matrix([[4, 5, 6]]))
		self.v3 = Algorithm(Matrix([[1], [2], [3]]))
		
		self.c1 = Algorithm(Matrix([[10]]))
		self.c2 = Algorithm(Matrix([[0]]))
		
		self.m1 = Algorithm(Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
		self.m2 = Algorithm(Matrix([[4, 1, -7, 2], [-1, 9, 6, 3]]))
		self.m3 = Algorithm(Matrix([[8, -3, 1], [4, -6, 2], [7, 3, 5], [-2, -5, 1]]))
		self.m4 = Algorithm(Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
		self.m5 = Algorithm(Matrix([[6, 1, 1], [4,-2, 5], [2, 8, 7]]))
		self.m6 = Algorithm(Matrix([[1, 3, 3], [1, 4, 3], [1, 3, 4]]))
		self.m7 = Algorithm(Matrix([[1, 2, 3], [0, 1, 4], [5, 6, 0]]))
		
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

test = False

if __name__ == "__main__":
	if test:
		unittest.main()
	else:
		#print "Hooray"
		numsys0 = NumberSystem(Matrix([[10]]), {0,1,2,3,4,5,6,7,8,9}) 						#(0)-(0)
		numsys1 = NumberSystem(Matrix([[-1,-1],[1,-1]]), {(0,0),(1,0)}) 					#(0,0)-(0,0)
		numsys2 = NumberSystem(Matrix([[3]]), {-2,0,2})										#(-1)-(1)
		numsys3 = NumberSystem(Matrix([[0,-3],[1,0]]), {(0,0),(1,0),(-1,1)})				#
		#numsys4a = NumberSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(2,0),(3,0),(4,0)})	#(0,-2)-(2,0)
		numsys4a = NumberSystem(Matrix([[2,-1],[1,2]]), jcan=1)	#(0,-2)-(2,0)
		numsys4b1 = NumberSystem(Matrix([[3,0],[0,3]]), {(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)}) #(-0.5,-0.5)-(0.5,0.5)
		numsys4b2 = NumberSystem(Matrix([[3,0],[0,3]]), {(0,0),(1,0),(2,0),(0,1),(0,2),(1,2),(2,1),(-1,2),(-2,1)}) #(-1,0)-(1,1)
		numsys5 = NumberSystem(Matrix([[1,-2],[1,1]]), {(0,0),(1,0),(-1,0)})				#
		numsys7 = NumberSystem(Matrix([[-3,-1],[1,-3]]), {(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(6,0),(7,0),(8,0),(9,0)})				#?
		numsys8a = NumberSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(0,1),(0,-1),(-6,-5)})	#
		numsys8b = NumberSystem(Matrix([[2,-1],[1,2]]), {(0,0),(1,0),(0,1),(0,-1),(-2,-3)})
		
		#print Encoder("Phyllorhiza_punctata.txt", True)
		
		print NumberSystem(Matrix.companion([4, -4, 1]), jcan=1).necessary_condition()

