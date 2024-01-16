import numpy as np
from quad_program import quad_program_for_perching as QPP
from quad_prog import quad_program_for_return as QPR
import matplotlib.pyplot as plt
import scipy.io

obstacle_x_y = [0 ,  0] 

S_Terminal = 1

floor_height = -1

# used for terminal position constraint
staging_area_frontal = [ -1.5,    -0.5,     0.5,    -0.5,   0.5]
#                       x_max    y_min    y_max    z_min   z_max
# final_vel = np.mat([[0.2, 0, 0]]).T

kappa_max_here = 2
# used for define the best targeted flight begining position.
sweet_spot = [ -1.5,     0,     0]

initial_pos = [    0,   0.5,     0]
initial_vel = [  0.2,     0,     0]

poly_order = 7

qdxy = QPR (poly_order, 1)
qdz =  QPP (poly_order, 1)      

Coef = [0, 1, 1, 1]
vel_con = 0.15
#QPP:
def exert_con_on_vel_symmetry(qpp_thing:QPP, sample_num,  constraint_value):
    for i in range(sample_num):
        qpp_thing.add_symmetry_vel_constraint_1_by_1(qpp_thing.Tend *i /sample_num, constraint_value)

def exert_con_on_vel_upp(qpp_thing:QPP, sample_num,  constraint_value):
    for i in range(sample_num):
        qpp_thing.add_upper_vel_constraint_1_by_1(qpp_thing.Tend *i /sample_num, constraint_value)

def exert_con_on_vel_low(qpp_thing:QPP, sample_num,  constraint_value):
    for i in range(sample_num):
        qpp_thing.add_lower_vel_constraint_1_by_1(qpp_thing.Tend *i /sample_num, constraint_value)

def exert_con_on_ob(qpp_thing:QPP, sample_num, ob_p, mu):
    for i in range(sample_num):
        qpp_thing.add_ob_at_t(ob_p, mu, qpp_thing.Tend * i /sample_num)

#QPR:
def exert_con_on_vel_symmetry(qpp_thing:QPR, sample_num,  constraint_value_x, constraint_value_y):
    for i in range(sample_num):
        qpp_thing.add_symmetry_vel_constraint_1_by_1(qpp_thing.Tend *i /sample_num,
                                                     constraint_value_x, constraint_value_y)

def exert_con_on_vel_upp(qpp_thing:QPR, sample_num,  constraint_value_x, constraint_value_y):
    for i in range(sample_num):
        qpp_thing.add_upper_vel_constraint_1_by_1(qpp_thing.Tend *i /sample_num, 
                                                  constraint_value_x, constraint_value_y)

def exert_con_on_vel_low(qpp_thing:QPR, sample_num,  constraint_value_x, constraint_value_y):
    for i in range(sample_num):
        qpp_thing.add_lower_vel_constraint_1_by_1(qpp_thing.Tend *i /sample_num, 
                                                  constraint_value_x, constraint_value_y)

def exert_con_on_ob(qpp_thing:QPR, sample_num, ob_px, ob_py, mu):
    for i in range(sample_num):
        qpp_thing.add_ob_at_t(ob_px, ob_py, mu, qpp_thing.Tend * i /sample_num)

qdxy.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])
qdz. set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])

qdxy.add_pos_constraint_1_by_1(0., initial_pos[0], initial_pos[1])
# qdx.add_upper_pos_constraint_1_by_1(qdx.Tend, initial_pos[0])
qdxy.add_vel_constraint_1_by_1(0., initial_vel[0], initial_vel[1])
qdxy.add_vel_constraint_1_by_1(qdxy.Tend, initial_vel[0], initial_vel[1])

# qdxy.add_pos_constraint_1_by_1(qdxy.Tend, sweet_spot[0], sweet_spot[1])

# exert_con_on_ob(qdxy, 20, obstacle_x_y[0], obstacle_x_y[1], 1)

qdxy.add_att_at_t(sweet_spot[0],sweet_spot[1], 100, qdxy.Tend)
qdxy.add_curv_Q(0.01, 1e6)

qdz.add_pos_constraint_1_by_1(0., initial_pos[2])
# qdz.add_vel_constraint_1_by_1(0.,0.)

# qdz.add_att_at_t(sweet_spot[2], 20, qdz.Tend)

##
qdxy.add_upper_pos_constraint_1_by_1(qdxy.Tend, staging_area_frontal[0], staging_area_frontal[2])

qdxy.add_lower_pos_constraint_1_by_1(qdxy.Tend, [], staging_area_frontal[1])

qdz.add_lower_pos_constraint_1_by_1(qdz.Tend, staging_area_frontal[3])
qdz.add_upper_pos_constraint_1_by_1(qdz.Tend, staging_area_frontal[4])


xy_traj = qdxy.opt()

z_traj = qdz.opt()

print('xy_traj:',xy_traj)
print('z_traj:',z_traj)

data = {'xy_poly_res': xy_traj, 'z_poly_res': z_traj,'T':1}
scipy.io.savemat('NewReturnTest.mat', data)