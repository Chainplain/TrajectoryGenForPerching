import numpy as np
import sympy as sp

def diff_mat(symbol_mat):
    t = sp.Symbol('t')
    [row,col]  = output_mat.shape
    output_mat = np.tile( t,(row,col))
    for i in range(row):
        for j in range(col):
            output_mat[i,j] = sp.diff( output_mat[i,j] , t)
    return output_mat

def Indef_int_mat(symbol_mat):
    output_mat = symbol_mat
    [row,col]  = output_mat.shape
    for i in range(row):
        for j in range(col):
            output_mat[i,j] = sp.Poly.integrate( output_mat[i,j] , t)
    return output_mat

def subs_mat(symbol_mat, value):
    # output_mat = symbol_mat
    [row,col]  = symbol_mat.shape
    output_mat = np.zeros((row,col))
    for i in range(row):
        for j in range(col):
            output_mat[i,j] =  symbol_mat[i,j].subs(t,value)
    return output_mat



polynomial_order = 6

marker_for_unit_polynomial = np.eye(polynomial_order + 1)
        
list_of_unit_polynomials = []
        # Here we initialize them.
for i in range(polynomial_order + 1):
    single_polynomial = sp.Poly.from_list( marker_for_unit_polynomial[i,:], t)
    list_of_unit_polynomials.append(single_polynomial)

horizontal_polynomial_mat = np.mat( list_of_unit_polynomials )
print('horizontal_polynomial_mat: ', horizontal_polynomial_mat)

vel_horizontal_polynomial_mat = diff_mat(horizontal_polynomial_mat)
# for i in range(polynomial_order + 1):
#     vel_single_polynomial = sp.diff( horizontal_polynomial_mat[0,i], t)
#     vel_horizontal_polynomial_mat.append(vel_single_polynomial)

print('vel_horizontal_polynomial_mat: ', vel_horizontal_polynomial_mat)

acc_horizontal_polynomial_mat = diff_mat(vel_horizontal_polynomial_mat)
print('acc_horizontal_polynomial_mat: ', acc_horizontal_polynomial_mat)

vel_horizontal_polynomial_mat_from_int = Indef_int_mat(acc_horizontal_polynomial_mat)
print('vel_horizontal_polynomial_mat_from_int: ', vel_horizontal_polynomial_mat_from_int)

# print('acc_matrix_square: ',acc_horizontal_polynomial_mat.T * acc_horizontal_polynomial_mat)

# acc_horizontal_polynomial_mat_value = subs_mat(acc_horizontal_polynomial_mat, 3)

# print('acc_horizontal_polynomial_mat_value: ',acc_horizontal_polynomial_mat_value)


