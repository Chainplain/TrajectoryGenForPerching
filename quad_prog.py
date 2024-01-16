#written by Chainplain at 2023-9-18 14:50:18
from cvxopt import matrix, solvers
import numpy as np
import sympy as sp
from sympy import symbols, Expr
# from qpsolvers import solve_qp

class quad_program_for_return:
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
    sol=solvers.qp(Q, p, None, None, A, b)
    ```
    """
    mu_vel  = 1
    mu_acc  = 1
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
        self. Tend = T_end
        self. poly_order = polynomial_order
        # It is a time symbol, since the polynomial is a function of time, this is just 
        # all the symbol we need.

        marker_for_unit_polynomial = np.eye(polynomial_order + 1)
        
        list_of_unit_polynomials = []
        # Here we initialize them.
        for i in range(polynomial_order + 1):
            single_polynomial = self.t  ** (polynomial_order  - i)
            list_of_unit_polynomials.append(single_polynomial)

        list_of_unit_polynomials_minus = []
        for i in range(polynomial_order + 1):
            single_polynomial = self.t  ** (polynomial_order  - i) #(0.8 *(self.t - 0.5)) ** (polynomial_order  - i)
            list_of_unit_polynomials_minus.append(single_polynomial)
        print('list_of_unit_polynomials_minus:', list_of_unit_polynomials_minus)
        print('list_of_unit_polynomials_minus test:', sp.diff(list_of_unit_polynomials_minus[0],self.t))
        # horizontal_polynomial_mat shown as  [[Poly(1.0*t**2, t, domain='RR') Poly(1.0*t, t, domain='RR')
        # Poly(1.0, t, domain='RR')]] that is, it is a nx1 vector indeed
        self. pos_horizontal_polynomial_mat = np.mat( list_of_unit_polynomials )
        self. vel_horizontal_polynomial_mat = self.diff_mat(self. pos_horizontal_polynomial_mat)
        self. acc_horizontal_polynomial_mat = self.diff_mat(self. vel_horizontal_polynomial_mat)
        self. jerk_horizontal_polynomial_mat = self.diff_mat(self. acc_horizontal_polynomial_mat)
        self. snap_horizontal_polynomial_mat = self.diff_mat(self. jerk_horizontal_polynomial_mat)

        self. pos_horizontal_polynomial_mat_minus = np.mat( list_of_unit_polynomials_minus )
        self. vel_horizontal_polynomial_mat_minus = self.diff_mat(self. pos_horizontal_polynomial_mat_minus)
        # print('vel_horizontal_polynomial_mat_minus:', self. vel_horizontal_polynomial_mat_minus)
        self. acc_horizontal_polynomial_mat_minus = self.diff_mat(self. vel_horizontal_polynomial_mat_minus)
        self. jerk_horizontal_polynomial_mat_minus = self.diff_mat(self. acc_horizontal_polynomial_mat_minus)
        self. snap_horizontal_polynomial_mat_minus = self.diff_mat(self. jerk_horizontal_polynomial_mat_minus)

        # print('self. pos_horizontal_polynomial_mat:', self. pos_horizontal_polynomial_mat)
        # print('self. snap_horizontal_polynomial_mat:', self. snap_horizontal_polynomial_mat)
        self. int_vel_square_mat =   self.Calc_SQ_double(self. vel_horizontal_polynomial_mat,
                                                     self. vel_horizontal_polynomial_mat)
        self. int_acc_square_mat =   self.Calc_SQ_double(self. acc_horizontal_polynomial_mat,
                                                     self. acc_horizontal_polynomial_mat)
        self. int_jerk_square_mat =  self.Calc_SQ_double(self. jerk_horizontal_polynomial_mat,
                                                     self. jerk_horizontal_polynomial_mat)
        self. int_snap_square_mat =  self.Calc_SQ_double(self. snap_horizontal_polynomial_mat,
                                                     self. snap_horizontal_polynomial_mat)

        # print('self. vel_square_mat:',self. vel_square_mat)
        self. int_total_square_mat= self.mu_vel  * self. int_vel_square_mat +\
                                    self.mu_acc  * self. int_acc_square_mat +\
                                    self.mu_jerk * self. int_jerk_square_mat +\
                                    self.mu_snap * self. int_snap_square_mat
        
        # renew the definition of Q P A b G h
        self. I_4_avoid_singularity = np.mat(np.eye( (self.poly_order + 1) * 2 ) )* 1e3
        self.Q = matrix( self. subs_mat(self. int_total_square_mat, self.Tend) -\
                         self. subs_mat(self. int_total_square_mat, 0 ) ) +\
                              self. I_4_avoid_singularity
        
        # self.Q = matrix(self.Q, tc='d')
        print('Q_mention:', np.shape(self.Q ))

        self.p = np.mat([0.0] * (2 * (self.poly_order + 1))).T

        self.p =  matrix(self.p, tc='d')

        self.A = []
        self.b = []
        self.G = []
        self.h = []

    def set_Tend(self, T_end):
        self.Tend = T_end
        self.Q = self. subs_mat(self. int_total_square_mat, self.Tend) -\
                 self. subs_mat(self. int_total_square_mat, 0)+\
                              self. I_4_avoid_singularity
        # self.Q = matrix(self.Q, tc='d')
        self.p = np.mat([0.0] * (2 * (self.poly_order + 1))).T
        self.p =  matrix(self.p, tc='d')

    def add_ob_at_t(self, ob_px, ob_py, mu, t):
        pos_mat       = self. subs_mat(self. pos_horizontal_polynomial_mat, t)
        Q_non_trivial = - matrix( pos_mat. T * pos_mat ) #* (self.Tend - self.t) / self.Tend 
        Q_upper       = np.hstack((Q_non_trivial, np.zeros([self. poly_order  + 1,self. poly_order  + 1])))
        Q_lower       = np.hstack((np.zeros([self. poly_order  + 1,self. poly_order  + 1]), Q_non_trivial))
        Q_add         = matrix(np.vstack((Q_upper, Q_lower)), tc='d')

        p_add         = matrix( np.vstack(( ob_px * pos_mat.T, ob_py * pos_mat.T) ), tc='d')  #* (self.Tend - self.t) / self.Tend 
        self.Q        = self.Q + mu * Q_add
        self.p        = self.p + mu * p_add

    def add_att_at_t(self, ob_px, ob_py, mu, t):
        pos_mat       = self. subs_mat(self. pos_horizontal_polynomial_mat, t)
        Q_non_trivial = matrix( pos_mat. T * pos_mat ) #* (self.Tend - self.t) / self.Tend 
        Q_upper       = np.hstack((Q_non_trivial, np.zeros([self. poly_order  + 1,self. poly_order  + 1])))
        Q_lower       = np.hstack((np.zeros([self. poly_order  + 1,self. poly_order  + 1]), Q_non_trivial))
        Q_add         = matrix(np.vstack((Q_upper, Q_lower)), tc='d')

        p_add         = - matrix( np.vstack( (ob_px * pos_mat.T, ob_py * pos_mat.T)) , tc='d')  #* (self.Tend - self.t) / self.Tend 
        self.Q        = self.Q + mu * Q_add
        self.p        = self.p + mu * p_add

    def add_curv_Q(self, delta_T, mu):
        pos_now             =   self. pos_horizontal_polynomial_mat
        pos_at_delta_T      =   self. subs_mat(pos_now, self.t + delta_T)
        pos_at_delta_2T     =   self. subs_mat(pos_now, self.t + 2 * delta_T)

        first_mat_item      =   ((pos_now - pos_at_delta_T).T  *  (pos_now - pos_at_delta_T))
        second_mat_item     =    ((pos_at_delta_T - pos_at_delta_2T).T  *  (pos_at_delta_T - pos_at_delta_2T))
        third_mat_item      =    -((pos_now - pos_at_delta_2T).T  *  (pos_now - pos_at_delta_2T))

        total_mat           =   first_mat_item + second_mat_item + third_mat_item
        int_total_mat       =   self.Indef_int_mat( total_mat )
        int_total_mat_float       =   self. subs_mat(int_total_mat, self.Tend - 2 * delta_T) -\
                                self. subs_mat(int_total_mat, 0)

        Q_upper             =   np.hstack((int_total_mat_float, np.zeros([self. poly_order  + 1,self. poly_order  + 1])))
        Q_lower             =   np.hstack((np.zeros([self. poly_order  + 1,self. poly_order  + 1]), int_total_mat_float))
        Q_add               =   np.vstack((Q_upper, Q_lower))
       
        
        Q_add               =   matrix(Q_add, tc='d')
        print('Curve mu *Q_add:',mu *Q_add)
        self.Q              =   self.Q + mu * Q_add
        

    def set_coeff(self, m_vel, m_acc, m_jerk, m_snap):
        if m_vel is not None:
            self.mu_vel = m_vel
        if m_acc is not None:
            self.mu_acc = m_acc
        if m_jerk is not None:
            self.mu_jerk = m_jerk 
        if m_snap is not None:  
            self.mu_snap = m_snap

        self. int_total_square_mat= self.mu_vel * self.  int_vel_square_mat +\
                                    self.mu_acc * self.  int_acc_square_mat +\
                                    self.mu_jerk * self. int_jerk_square_mat +\
                                    self.mu_snap * self. int_snap_square_mat
        
        # renew the definition of Q P A b G h
        self.Q = self. subs_mat(self. int_total_square_mat, self.Tend) -\
                 self. subs_mat(self. int_total_square_mat, 0)+\
                 self. I_4_avoid_singularity
        # self.Q = self.Q + 1000 * np.eye( 2 * (self.poly_order + 1))
        
        self.Q = matrix(self.Q, tc='d')
        self.p = np.mat([0.0] * (2 * (self.poly_order + 1))).T
        self.p =  matrix(self.p, tc='d')

        
    
    def Calc_SQ_double(self,  horizontal_polynomial_mat,horizontal_polynomial_mat_minus):
        Transpose_mul_origin_Mat = horizontal_polynomial_mat.T @ horizontal_polynomial_mat
        Transpose_mul_origin_Mat_minus = horizontal_polynomial_mat_minus.T @ horizontal_polynomial_mat_minus
        Transpose_Mat_uppper = np.hstack( ( Transpose_mul_origin_Mat, 
                                           np.zeros((self.poly_order+1, self.poly_order+1))) )
        Transpose_Mat_lower = np.hstack( ( np.zeros((self.poly_order+1, self.poly_order+1)), 
                                           Transpose_mul_origin_Mat_minus) )
        Transpose_Mat = np.vstack( (Transpose_Mat_uppper, Transpose_Mat_lower) )

        Int_Mat = self.Indef_int_mat(Transpose_Mat)
        print('Int_Mat:', Int_Mat)
        return Int_Mat
    
    def Calc_SQ(self,ssymbol, opt_order):
        if not (opt_order == 1 or opt_order == 2 or \
                opt_order == 3 or opt_order == 4):
            raise ValueError("Opt_order must be 1, 2, 3, and 4.")
        t = ssymbol
        output_mat = np.tile( t,(self.poly_order+1, self.poly_order+1))
        for l in range(self.poly_order+1):
            for k in range(self.poly_order+1):
                if l >= opt_order and k >= opt_order:
                    gain = 1.0
                    for m in range(opt_order):
                        gain = gain * float(l - m) * float(k - m)
                    output_mat[self.poly_order-l,self.poly_order-k]=gain /float(l + k - 2 * opt_order + 1) \
                        * t ** float(l + k - 2 * opt_order + 1)
                else:
                    output_mat[self.poly_order-l,self.poly_order-k]=0. * t
        return output_mat

    def clear_constraints(self):
        self.A = []
        self.b = []
        self.G = []
        self.h = []
        
    def add_pos_constraint_1_by_1(self, tic, posx, posy):
        unit_pos_at_tic = self. subs_mat(self. pos_horizontal_polynomial_mat, tic)
        unit_posx_at_tic_with_zeros = np. hstack( (unit_pos_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_posy_at_tic = np. hstack( (np.zeros([1,self.poly_order +1]), unit_pos_at_tic))
        if  posx != []:
            self.A. append( unit_posx_at_tic_with_zeros[0,:].tolist() )
            self.b. append( posx )
        if  posy != []:
            self.A. append( zeros_with_unit_posy_at_tic[0,:].tolist() )
            self.b. append( posy )
        
    def add_upper_pos_constraint_1_by_1(self, tic, posx, posy):
        unit_pos_at_tic = self. subs_mat(self. pos_horizontal_polynomial_mat, tic)
        unit_posx_at_tic_with_zeros = np. hstack( (unit_pos_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_posy_at_tic = np. hstack( (np.zeros([1,self.poly_order +1]), unit_pos_at_tic))
        if posx != []:
            self.G. append( unit_posx_at_tic_with_zeros[0,:].tolist() )
            self.h. append( posx )
        if posy != []:
            self.G. append( zeros_with_unit_posy_at_tic[0,:].tolist() )
            self.h. append( posy )

    def add_lower_pos_constraint_1_by_1(self, tic, posx, posy):
        unit_pos_at_tic = self. subs_mat(self. pos_horizontal_polynomial_mat, tic)
        unit_posx_at_tic_with_zeros = np. hstack( (unit_pos_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_posy_at_tic = np. hstack( (np.zeros([1,self.poly_order +1]), unit_pos_at_tic))
        if posx != []:
            self.G. append( (-unit_posx_at_tic_with_zeros[0,:]).tolist() )
            self.h. append( -posx )
        if posy != []:
            self.G. append( (-zeros_with_unit_posy_at_tic[0,:]).tolist() )
            self.h. append( -posy )

    def add_symmetry_pos_constraint_1_by_1(self, tic, posx, posy):
        self.add_upper_pos_constraint_1_by_1(tic, posx, posy)
        self.add_lower_pos_constraint_1_by_1(tic, -posx, -posy)

    def add_vel_constraint_1_by_1(self, tic, velx, vely):
        unit_vel_at_tic = self. subs_mat(self. vel_horizontal_polynomial_mat, tic)
        unit_velx_at_tic_with_zeros = np. hstack( (unit_vel_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_vely_at_tic = np. hstack( (np.zeros([1,self.poly_order +1]), unit_vel_at_tic))
        if velx != []:
            self.A. append( unit_velx_at_tic_with_zeros[0,:].tolist() )
            self.b. append( velx )
        if vely != []:
            self.A. append( zeros_with_unit_vely_at_tic[0,:].tolist() )
            self.b. append( vely )


    def add_upper_vel_constraint_1_by_1(self, tic, velx, vely):
        unit_vel_at_tic = self. subs_mat(self. vel_horizontal_polynomial_mat, tic)
        unit_velx_at_tic_with_zeros = np. hstack( (unit_vel_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_vely_at_tic = np. hstack( (np.zeros([1,self.poly_order +1]), unit_vel_at_tic))

        if velx != []:
            self.G. append( unit_velx_at_tic_with_zeros[0,:].tolist() )
            self.h. append( velx )
        if vely != []:
            self.G. append( zeros_with_unit_vely_at_tic[0,:].tolist() )
            self.h. append( vely )

    def add_lower_vel_constraint_1_by_1(self, tic, velx, vely):
        unit_vel_at_tic = self. subs_mat(self. vel_horizontal_polynomial_mat, tic)
        unit_velx_at_tic_with_zeros = np. hstack( (unit_vel_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_vely_at_tic = np. hstack( (np.zeros([1,self.poly_order +1]), unit_vel_at_tic))
        
        if velx != []:
            self.G. append( (-unit_velx_at_tic_with_zeros[0,:]).tolist() )
            self.h. append( -velx )
        if vely != []:
            self.G. append( (-zeros_with_unit_vely_at_tic[0,:]).tolist() )
            self.h. append( -vely )

    def add_symmetry_vel_constraint_1_by_1(self, tic, velx, vely):
        self.add_upper_vel_constraint_1_by_1(tic, velx, vely)
        self.add_lower_vel_constraint_1_by_1(tic, -velx, -vely)

    def add_acc_constraint_1_by_1(self, tic, accx, accy):
        unit_acc_at_tic = self. subs_mat(self. acc_horizontal_polynomial_mat, tic)
        unit_accx_at_tic_with_zeros = np. hstack( (unit_acc_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_accy_at_tic = np. hstack( (np.zeros([1,self.poly_order +1]), unit_acc_at_tic))

        if accx != []:
            self.A. append( unit_accx_at_tic_with_zeros[0,:].tolist() )
            self.b. append( accx )
        if accy != []:
            self.A. append( zeros_with_unit_accy_at_tic[0,:].tolist() )
            self.b. append( accy )

    def add_upper_acc_constraint_1_by_1(self, tic, accx, accy):
        unit_acc_at_tic = self. subs_mat(self. acc_horizontal_polynomial_mat, tic)
        unit_accx_at_tic_with_zeros = np. hstack( (unit_acc_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_accy_at_tic = np. hstack( ( np.zeros([1,self.poly_order +1]), unit_acc_at_tic))

        if accx != []:
            self.G. append( unit_accx_at_tic_with_zeros[0,:].tolist() )
            self.h. append( accx )
        if accy != []:
            self.G. append( zeros_with_unit_accy_at_tic[0,:].tolist() )
            self.h. append( accy )

    def add_lower_acc_constraint_1_by_1(self, tic, accx, accy):
        unit_acc_at_tic = self. subs_mat(self. acc_horizontal_polynomial_mat, tic)
        unit_accx_at_tic_with_zeros = np. hstack( (unit_acc_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_accy_at_tic = np. hstack( (np.zeros([1,self.poly_order +1]), unit_acc_at_tic))
        if accx != []:
            self.G. append((-unit_accx_at_tic_with_zeros[0,:]).tolist() )
            self.h. append( -accx )
        if accy != []:
            self.G. append( -zeros_with_unit_accy_at_tic[0,:].tolist() )
            self.h. append( -accy )

    def add_symmetry_acc_constraint_1_by_1(self, tic, accx, accy):
        self.add_upper_acc_constraint_1_by_1(tic, accx, accy)
        self.add_lower_acc_constraint_1_by_1(tic, -accx, -accy)

    def add_jerk_constraint_1_by_1(self, tic, jerkx, jerky):
        unit_jerk_at_tic = self. subs_mat(self. jerk_horizontal_polynomial_mat, tic)
        unit_jerkx_at_tic_with_zeros = np. hstack( (unit_jerk_at_tic, np.zeros([1,self.poly_order +1])))
        zeros_with_unit_jerky_at_tic = np. hstack( (np.zeros([1,self.poly_order +1]), unit_jerk_at_tic))

        if jerkx != []:
            self.A. append( unit_jerkx_at_tic_with_zeros[0,:].tolist() )
            self.b. append( jerkx )
        if jerky != []:
            self.A. append( zeros_with_unit_jerky_at_tic[0,:].tolist() )
            self.b. append( jerky )

    def opt(self):
        # solvers.options['show_progress'] = False
        A_mat = matrix(np.mat(self.A))
        b_mat = matrix(np.mat(self.b).T)
     
        
        print('Q :', self. Q)
        print('p :',  self. p)
        print('A_mat :',  A_mat)
        print('b_mat :',  b_mat)
        if  not len(self.G):
            sol=solvers.qp(self.Q, self.p, None, None, A=A_mat, b=b_mat)
        else:
            G_mat = matrix(np.mat(self.G))
            h_mat = matrix(np.mat(self.h).T)
            print('G_mat :',  G_mat)
            print('h_mat :',  h_mat)

            sol=solvers.qp(self.Q, self.p, G_mat, h_mat, A_mat, b_mat)



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
                output_mat[i,j] = sp.integrate( symbol_mat[i,j] , self. t)
        return output_mat
    
    def subs_mat(self, the_mat, value):
    # output_mat = symbol_mat
        [row,col]  = the_mat.shape
        # t = sp.Symbol('t')
        if isinstance(value, Expr):
            t = sp.Symbol('t')
            output_mat = np.tile( t,(row,col))
        else:
            output_mat = np.zeros((row,col))
        for i in range(row):
            for j in range(col):
                output_mat[i,j] =  the_mat[i,j].subs(self. t, value)
        return output_mat