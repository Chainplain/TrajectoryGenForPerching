# too slow

import numpy as np
from quad_program import quad_program_for_perching as QPP
import matplotlib.pyplot as plt
import scipy.io
import functools
from scipy.optimize import minimize

# used for obstacle aboidance constraint -0.3， 0.3， -0.6， -0.5， -1， 0.1
obstacle_box = [-0.3,    0.3,    -0.6,   -0.5,   -1,    0.1] 
#               x_min   x_max   y_min   y_max   z_min   z_max

S_Terminal = 1

floor_height = -1

# used for terminal position constraint
staging_area_frontal = [ -1.5,    -0.5,     0.5,    -0.5,   0.5]
#                       x_max    y_min    y_max    z_min   z_max
# final_vel = np.mat([[0.2, 0, 0]]).T

kappa_max_here = 2



# used for define the best targeted flight begining position.
sweet_spot = [ 0,    0,     -1.5]

initial_pos = [ 0,    0.5,     0]

poly_order = 7
para_num_per_dim = poly_order + 1

qdx = QPP (poly_order, S_Terminal)
qdy = QPP (poly_order, S_Terminal)
qdz = QPP (poly_order, S_Terminal)

Coef = [1, 1, 1, 1]

qdx.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])
qdy.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])
qdz.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])


constraint_dict = []

def get_pos( para_vec, time):
    basis = [time**7, time**6, time**5, time**4, time**3, time**2, time**1, 1]
    basis_mat = np.mat(basis)
    Pos = para_vec.T * basis_mat.T
    return Pos

def get_pos_XYZ(para, t):
    para_x = []
    for i in range(para_num_per_dim):
        para_x.append( para[i])
    para_x_vec = np.mat( para_x ).T
    Pos_x_mat = get_pos( para_x_vec, t)
    Pos_x = Pos_x_mat[0,0]

    para_y = []
    for i in range(para_num_per_dim):
        para_y.append( para[i + para_num_per_dim] )
    para_y_vec = np.mat( para_y ).T
    Pos_y_mat = get_pos( para_y_vec, t)
    # print('Pos_y_mat:',Pos_y_mat)
    Pos_y = Pos_y_mat[0,0]

    para_z = []
    for i in range(para_num_per_dim):
        para_z.append( para[i + para_num_per_dim*2] )
    para_z_vec = np.mat( para_z ).T
    Pos_z_mat = get_pos( para_z_vec, t)
    Pos_z = Pos_z_mat[0,0]

    return [Pos_x, Pos_y, Pos_z]

def final_position_constraint_ineq(para):
    [Pos_x, Pos_y, Pos_z] = get_pos_XYZ(para, S_Terminal)

    if (Pos_x <= staging_area_frontal[0] and\
        Pos_y >= staging_area_frontal[1] and\
        Pos_y <= staging_area_frontal[2] and\
        Pos_z >= staging_area_frontal[3] and\
        Pos_z <= staging_area_frontal[4] ):
        print('IN Stage')
        return 1
    else:
        return -1
    
# import functools

# def add(x, y):
#     return x + y

# add_5 = functools.partial(add, 5)

# result = add_5(3)
    
def curvature_constraint_ineq_at_t(para, kappa_max, t_0, t_1, t_2):
    [Pos_x_0, Pos_y_0, Pos_z_0] = get_pos_XYZ(para, t_0)
    [Pos_x_1, Pos_y_1, Pos_z_1] = get_pos_XYZ(para, t_1)
    [Pos_x_2, Pos_y_2, Pos_z_2] = get_pos_XYZ(para, t_2)
    a = np.sqrt( (Pos_x_0 - Pos_x_1)**2 + (Pos_y_0 - Pos_y_1)**2)
    b = np.sqrt( (Pos_x_2 - Pos_x_1)**2 + (Pos_y_2 - Pos_y_1)**2)
    c = np.sqrt( (Pos_x_2 - Pos_x_0)**2 + (Pos_y_2 - Pos_y_0)**2)

    ab = (a + b) / 2
    return ( c - 2 * ab / ( kappa_max * np.sqrt( ab**2 + 1 / ( kappa_max**2 ) ) ) )


def obstacle_constraint_ineq_at_t(para, t):
    [Pos_x, Pos_y, Pos_z] = get_pos_XYZ(para, t)
    if (Pos_x <= obstacle_box[0] or\
        Pos_x >= obstacle_box[1] or\
        Pos_y >= obstacle_box[2] or\
        Pos_y <= obstacle_box[3] or\
        Pos_z >= obstacle_box[4] or\
        Pos_z <= obstacle_box[5] )\
    and Pos_z > floor_height :
        return  1
    else:
        return -1

def initial_pos_x_constraint_eq(para):
    [Pos_x, Pos_y, Pos_z] = get_pos_XYZ(para, 0)
    cons = (Pos_x - initial_pos[0])**2
    return cons

def initial_pos_y_constraint_eq(para):
    [Pos_x, Pos_y, Pos_z] = get_pos_XYZ(para, 0)
    cons = (Pos_y - initial_pos[1])**2
    return cons

def initial_pos_z_constraint_eq(para):
    [Pos_x, Pos_y, Pos_z] = get_pos_XYZ(para, 0)
    cons = (Pos_z - initial_pos[2])**2
    return cons
    
def objective(para):
    cost = 0

    para_x = []
    for i in range(para_num_per_dim):
        para_x.append( para[i])
    para_x_vec = np.mat( para_x ).T
    Qx = qdx. subs_mat(qdx. int_total_square_mat, qdx.Tend) - qdx. subs_mat(qdx. int_total_square_mat, 0)
    cost_x_mat = para_x_vec.T * Qx * para_x_vec
    cost = cost + cost_x_mat[0,0]

    para_y = []
    for i in range(para_num_per_dim):
        para_y.append( para[i + para_num_per_dim])
    para_y_vec = np.mat( para_y ).T
    Qy = qdy. subs_mat(qdy. int_total_square_mat, qdy.Tend) - qdy. subs_mat(qdy. int_total_square_mat, 0)
    cost_y_mat = para_y_vec.T * Qy * para_y_vec
    cost = cost + cost_y_mat[0,0]

    para_z = []
    for i in range(para_num_per_dim):
        para_z.append( para[i + para_num_per_dim*2] )
    para_z_vec = np.mat( para_z ).T
    Qz = qdz. subs_mat(qdz. int_total_square_mat, qdz.Tend) - qdz. subs_mat(qdz. int_total_square_mat, 0)
    cost_z_mat = para_z_vec.T * Qz * para_z_vec
    cost = cost + cost_z_mat[0,0]

    [Pos_x, Pos_y, Pos_z] = get_pos_XYZ(para, S_Terminal)
    # pos_show = [Pos_x, Pos_y, Pos_z]
    # print('pos_show:',pos_show)
    mu_ss = 10

    cost_sweet_spot = mu_ss * ( (Pos_x - sweet_spot[0])**2 +\
                                (Pos_y - sweet_spot[1])**2 +\
                                (Pos_z - sweet_spot[2])**2 )
    
    cost = cost + cost_sweet_spot
    
    return cost


sample_num = 5
# constraint_dict.append( {'type': 'ineq', 'fun': final_position_constraint_ineq} )

sample_interval = 1.0 / sample_num

# for i in range( 1, sample_num + 1 ):
#     obstacle_constraint_ineq = functools.partial(obstacle_constraint_ineq_at_t, t = sample_interval*i)
#     constraint_dict.append( {'type': 'ineq', 'fun': obstacle_constraint_ineq} )

# for i in range( 1, sample_num - 1):
#     curvature_constraint_ineq = functools.partial(curvature_constraint_ineq_at_t, 
#                                                 kappa_max = kappa_max_here,
#                                                 t_0 = sample_interval*i,
#                                                 t_1 = sample_interval* (i+1),
#                                                 t_2 = sample_interval* (i+1) )
#     constraint_dict.append( {'type': 'ineq', 'fun': curvature_constraint_ineq} )

constraint_dict.append( {'type': 'eq', 'fun': initial_pos_x_constraint_eq} )
constraint_dict.append( {'type': 'eq', 'fun': initial_pos_y_constraint_eq} )
constraint_dict.append( {'type': 'eq', 'fun': initial_pos_z_constraint_eq} )


N = para_num_per_dim*3  # Replace 10 with the number of random integers you want to generate
low = -1  # Replace 0 with the lowest possible integer value
high = 1  # Replace 100 with the highest possible integer value

initial_para = np.random.randint(low, high, size = N)

result = minimize(objective, initial_para,   constraints=constraint_dict)

x_opt = result.x
obj_value = result.fun

print("Optimized design variables:", x_opt)
print("Optimized objective value:", obj_value)    









    



        

    


