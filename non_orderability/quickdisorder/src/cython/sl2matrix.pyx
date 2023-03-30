"""
Fast versions of certain heavily called functions.
"""
cdef extern from "matrix.h":
    ctypedef struct GL2CMatrix:
        double real[2][2]
        double imag[2][2]

    void copy_GL2C(GL2CMatrix* A, GL2CMatrix* B)
    void zero_out_GL2C(GL2CMatrix* A)
    void identity_GL2C(GL2CMatrix* A)
    void multiply_GL2C(GL2CMatrix* A, GL2CMatrix* B, GL2CMatrix* C)
    void trace_GL2C(GL2CMatrix* A, double* real, double* imag)
    void hash_GL2C(GL2CMatrix* A, long* hash, int bits)
    void inverse_SL2C(GL2CMatrix* A, GL2CMatrix* B)
    double norm_GL2(GL2CMatrix* A)
    int is_one(GL2CMatrix* A, int bits)


cdef copy_to_GL2CMatrix(M, GL2CMatrix* N):
    cdef int i, j
    for i in range(2):
        for j in range(2):
            entry = complex(M[i][j])
            N.real[i][j] = entry.real
            N.imag[i][j] = entry.imag

cdef GL2CMatrix_to_array(GL2CMatrix* M):
    cdef int i, j
    ans = []
    for i in range(2):
        row = []
        for j in range(2):
            row.append(complex(M.real[i][j], M.imag[i][j]))
        ans.append(row)
    return ans

cdef class DoubleGroupElement(object):
    cdef GL2CMatrix matrix
    cdef int min_bits_accuracy
    cpdef readonly word
    cdef long[8] _hashed_matrix
    cdef _hash
    
    def __cinit__(self, matrix, min_bits_accuracy, word=''):
        self.min_bits_accuracy, self.word = min_bits_accuracy, word
        if matrix is not None: 
            copy_to_GL2CMatrix(matrix, &self.matrix)
            self.set_hash()
        else:
            zero_out_GL2C(&self.matrix)
            self._hash = None

    cdef set_matrix(self, GL2CMatrix* M):
        copy_GL2C(M, &self.matrix)

    def norm(self):
        return norm_GL2(&self.matrix)

    def __mul__(DoubleGroupElement self, DoubleGroupElement other):
        cdef DoubleGroupElement ans
        ans = DoubleGroupElement(None, self.min_bits_accuracy,
                                self.word + other.word)
        multiply_GL2C(&self.matrix, &other.matrix, &ans.matrix)
        ans.set_hash()
        return ans


    def __richcmp__(DoubleGroupElement self, DoubleGroupElement other, int op):
        cdef int i
        cdef long s, o
        for i in range(8):
            s = self._hashed_matrix[i]
            o = other._hashed_matrix[i]
            if op == 0:
                c = (s < o)
            elif op == 1:
                c = (s <= o)
            elif op == 2:
                c = (s == o)
            elif op == 3:
                c = (s != o)
            elif op == 4:
                c = (s > o)
            else:
                c = (s >= o)
            if not c:
                return False
        return True

    def trace(self):
        cdef double real, imag
        trace_GL2C(&self.matrix, &real, &imag)
        return complex(real, imag)

    cdef set_hash(self):
        cdef long *ans = self._hashed_matrix
        hash_GL2C(&self.matrix, ans, self.min_bits_accuracy)
        hash_tuple = (ans[0], ans[1], ans[2], ans[3], ans[4], ans[5], ans[6], ans[7])
        self._hash = hash(hash_tuple)

    def is_one(self):
        return bool(is_one(&self.matrix, self.min_bits_accuracy))
        
    def __hash__(self):
        return self._hash

    def inverse(self):
        cdef DoubleGroupElement ans
        ans = DoubleGroupElement(None, self.min_bits_accuracy,
                                 self.word.swapcase()[::-1])
        inverse_SL2C(&self.matrix, &ans.matrix)
        ans.set_hash()
        return ans

    def _set_word(self, word):
        self.word = word
    
    def __repr__(self):
        A = GL2CMatrix_to_array(&self.matrix)
        entries = A[0] + A[1]
        args = tuple([self.word] + entries)
        return "<NGE: %s; %r %r %r %r>" % args

cdef class ProjectiveDoubleGroupElement(DoubleGroupElement):
    def __cmp__(DoubleGroupElement self, DoubleGroupElement other):
        cdef int i
        cdef int result = 0
        cdef int equal_mod_center = 1
        for i in range(8):
            s = self._hashed_matrix[i]
            o = other._hashed_matrix[i]
            if s < o:
                result = -1
                break
            if s > o:
                result = 1
                break
        if result == 0:
            return 0
        for i in range(8):
            if self._hashed_matrix[i] != -other._hashed_matrix[i]:
                equal_mod_center = 0
                break
        if equal_mod_center:
            return 0
        else:
            return result
