#written by Chainplain at 2023-9-18 14:50:18
from cvxopt import matrix, solvers
import numpy as np
import sympy as sp

class quad_program_for_perching:
    """
    # quad_program_for_perching 
    ## Dependencies
    - numpy
    - cvxopt \n

    ## Descirptions
    - ...
    ## Notations
    In this class we use the cvxopt, the basic quadratic problem can be solved as
    \n \n
    min  x 1/2 * x^T Q x + p^T x \n
    subject to   G x <= h \n
    A x == b \n

    Then one can use:
    ```Python
    import numpy
    from cvxopt import matrix
    Q = 2*matrix([ [2, .5], [.5, 1] ])
    p = matrix([1.0, 1.0])
    G = matrix([[-1.0,0.0],[0.0,-1.0]])
    h = matrix([0.0,0.0])
    A = matrix([1.0, 1.0], (1,2))
    b = matrix(1.0)
    sol=solvers.qp(Q, p, G, h, A, b)
    ## if there are A, b 
    sol = solvers.qp(P,q,G,h,A,b)
    sol=solvers.qp(Q, p, None, None, A, b)
    ```

    """
    mu_vel = 1
    mu_acc = 1
    mu_jerk = 1
    mu_snap = 1


    def __init__(self, polynomial_order,  T_end):
        """
        We need to obtain the coefficients due to different order derivatives. \n
        Here we define a unit polynomial as 
        /
        X**2 + X + 1,
        /
        that is all the coefficient as one, \n
        this is quite representive, because we barely can multiply the corresponding  \n
        coefficients and get any polynomial we want. \n
        """
        self. t = sp.Symbol('t')
        # It is a time symbol, since the polynomial is a function of time, this is just 
        # all the symbol we need.

        marker_for_unit_polynomial = np.eye(polynomial_order + 1)
        
        list_of_unit_polynomials = []
        # Here we initialize them.
        for i in range(polynomial_order + 1):
            single_polynomial = sp.Poly.from_list( marker_for_unit_polynomial[i,:], self.t)
            list_of_unit_polynomials.append(single_polynomial)

        # horizontal_polynomial_mat shown as  [[Poly(1.0*t**2, t, domain='RR') Poly(1.0*t, t, domain='RR')
        # Poly(1.0, t, domain='RR')]] 
        self. pos_horizontal_polynomial_mat = np.mat( list_of_unit_polynomials )
        
        self. vel_horizontal_polynomial_mat = self.diff_mat(self. pos_horizontal_polynomial_mat)
        # print('self. vel_horizontal_polynomial_mat:', self. vel_horizontal_polynomial_mat)
        self. acc_horizontal_polynomial_mat = self.diff_mat(self. vel_horizontal_polynomial_mat)
        self. jerk_horizontal_polynomial_mat = self.diff_mat(self. acc_horizontal_polynomial_mat)
        self. snap_horizontal_polynomial_mat = self.diff_mat(self. jerk_horizontal_polynomial_mat)

        print('self. pos_horizontal_polynomial_mat:', self. pos_horizontal_polynomial_mat)
        print('self. snap_horizontal_polynomial_mat:', self. snap_horizontal_polynomial_mat)
        vel_square_mat = self. vel_horizontal_polynomial_mat.T * self. vel_horizontal_polynomial_mat
        acc_square_mat = self. acc_horizontal_polynomial_mat.T * self. acc_horizontal_polynomial_mat
        jerk_square_mat = self. jerk_horizontal_polynomial_mat.T * self. jerk_horizontal_polynomial_mat
        snap_square_mat = self. snap_horizontal_polynomial_mat.T * self. snap_horizontal_polynomial_mat

        total_square_mat =  self.mu_vel * vel_square_mat +\
                            self.mu_acc * acc_square_mat +\
                            self.mu_jerk * jerk_square_mat +\
                            self.mu_snap * snap_square_mat
        
        self. int_total_square_mat = self.Indef_int_mat(total_square_mat)
        
        self.Q = matrix(self. subs_mat(self. int_total_square_mat, T_end) - self. subs_mat(self. int_total_square_mat, 0), tc='d')
        self.p = matrix([0.0] * (polynomial_order + 1), tc='d')

        self.A = []
        self.b = []
    def set_Tend(self, T_end):
        self.Q = matrix(self. subs_mat(self. int_total_square_mat, T_end) - self. subs_mat(self. int_total_square_mat, 0), tc='d')

    def clear_constraints(self):
        self.A = []
        self.b = []
        
    def add_pos_constraint_1_by_1(self, tic, pos):
        unit_pos_at_tic = self. subs_mat(self. pos_horizontal_polynomial_mat, tic)
        self.A. append( unit_pos_at_tic[0,:].tolist() )
        self.b. append( pos )

    def add_vel_constraint_1_by_1(self, tic, vel):
        unit_vel_at_tic = self. subs_mat(self. vel_horizontal_polynomial_mat, tic)
        self.A. append( unit_vel_at_tic[0,:].tolist() )
        self.b. append( vel )

    def add_acc_constraint_1_by_1(self, tic, acc):
        unit_acc_at_tic = self. subs_mat(self. acc_horizontal_polynomial_mat, tic)
        self.A. append( unit_acc_at_tic[0,:].tolist() )
        self.b. append( acc )

    def opt(self):
        A_mat = matrix(np.mat(self.A))
        b_mat = matrix(np.mat(self.b).T)
        sol=solvers.qp(self.Q, self.p, None, None, A_mat, b_mat)
        return sol['x']

    def diff_mat(self, symbol_mat):
        [row,col]  = symbol_mat.shape
        t = sp.Symbol('t')
        output_mat = np.tile( t,(row,col))

        for i in range(row):
            for j in range(col):
                output_mat[i,j] = sp.diff( symbol_mat[i,j] , self. t)
        return output_mat
    
    def Indef_int_mat(self, symbol_mat):
        [row,col]  = symbol_mat.shape
        t = sp.Symbol('t')
        output_mat = np.tile( t,(row,col))
        # t = sp.Symbol('t')
        [row,col]  = output_mat.shape
        for i in range(row):
            for j in range(col):
                output_mat[i,j] = sp.Poly.integrate( symbol_mat[i,j] , self. t)
        return output_mat
    
    def subs_mat(self, the_mat, value):
    # output_mat = symbol_mat
        [row,col]  = the_mat.shape
        # t = sp.Symbol('t')
        output_mat = np.zeros((row,col))
        for i in range(row):
            for j in range(col):
                output_mat[i,j] =  the_mat[i,j].subs(self. t, value)
        return output_mat