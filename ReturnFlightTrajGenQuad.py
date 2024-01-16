import numpy as np
from quad_program import quad_program_for_perching as QPP
import matplotlib.pyplot as plt
import scipy.io
from testDubin import Waypoint, calcDubinsPath, dubins_traj
import math

obstacle_x_y = [0 ,  0] 

S_Terminal = 1

floor_height = -1

# used for terminal position constraint
staging_area_frontal = [ -1.5,    -0.5,     0.5,    -0.5,   0.5]
#                       x_max    y_min    y_max    z_min   z_max
# final_vel = np.mat([[0.2, 0, 0]]).T

kappa_max_here = 2
# used for define the best targeted flight begining position.
sweet_pos = [ -1.5,     0,     0]
sweet_vel  = [  0.2,     0,     0]  

initial_pos = [  -0.6,   0,     0.5]
initial_vel = [ 0.2,     0,     0]


radius = 0.3 # turn radius unit:m

initial_theta = math.atan2( initial_vel[1], initial_vel[0]) / math.pi * 180
initial_wpt = Waypoint(initial_pos[0],  initial_pos[1], initial_theta)

sweet_theta = math.atan2( sweet_vel[1], sweet_vel[0]) / math.pi * 180
sweet_wpt = Waypoint(sweet_pos[0],  sweet_pos[1], sweet_theta)

param = calcDubinsPath(initial_wpt, sweet_wpt, radius)
Dubin_length = (param.seg_final[0]+param.seg_final[1]+param.seg_final[2])*param.turn_radius

Dubin_path_step = 1e-2
path = dubins_traj(param, Dubin_path_step)
# print('path:',path)

poly_order = 11

T_end = Dubin_length / 0.2;
qdx = QPP (poly_order, T_end)
qdy = QPP (poly_order, T_end)
qdz = QPP (poly_order, T_end)

Coef = [1, 1, 1, 1]
vel_con = 0.15

def exert_con_on_vel_symmetry(qpp_thing, sample_num,  constraint_value):
    for i in range(sample_num):
        qpp_thing.add_symmetry_vel_constraint_1_by_1(qpp_thing.Tend *i /sample_num, constraint_value)

def exert_con_on_vel_upp(qpp_thing, sample_num,  constraint_value):
    for i in range(sample_num):
        qpp_thing.add_upper_vel_constraint_1_by_1(qpp_thing.Tend *i /sample_num, constraint_value)

def exert_con_on_vel_low(qpp_thing, sample_num,  constraint_value):
    for i in range(sample_num):
        qpp_thing.add_lower_vel_constraint_1_by_1(qpp_thing.Tend *i /sample_num, constraint_value)

def exert_con_on_ob(qpp_thing, sample_num, constraint_value, mu):
    i_plus = 2
    for i in range(sample_num):
        mu_use = 1/ sample_num * mu * ((sample_num+i_plus) / (i+i_plus))**4
        if mu_use > mu:
            mu_use = mu
        qpp_thing.add_ob_at_t(constraint_value, mu_use , qpp_thing.Tend /2 * i /sample_num)

qdx.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])
qdy.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])
qdz.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])

qdx.add_pos_constraint_1_by_1(0., initial_pos[0])
# qdx.add_upper_pos_constraint_1_by_1(qdx.Tend, initial_pos[0])
qdx.add_vel_constraint_1_by_1(0., initial_vel[0])
qdx.add_vel_constraint_1_by_1(qdx.Tend, sweet_vel[0])

ob_mu = 5e1
exert_con_on_ob(qdx, 50, obstacle_x_y[0], ob_mu)
# qdx.add_att_at_t(sweet_pos[0], 20, qdx.Tend)

qdy.add_pos_constraint_1_by_1(0., initial_pos[1])
qdy.add_vel_constraint_1_by_1(0.,initial_vel[1])
qdy.add_vel_constraint_1_by_1(qdy.Tend, sweet_vel[1])
# print('qdy.end:',qdy.end)
exert_con_on_ob(qdy, 50, obstacle_x_y[1], ob_mu)
# qdy.add_att_at_t(sweet_pos[1], 20, qdy.Tend)

qdz.add_pos_constraint_1_by_1(0., initial_pos[2])
qdz.add_vel_constraint_1_by_1(0.,0.)
qdz.add_att_at_t(sweet_pos[2], 2, qdz.Tend)



path_length = len(path);
attract_mu = 1e2
for i in range(path_length):
    qdx.add_att_at_t(path[i,0], attract_mu, Dubin_path_step * i / 0.2)
    qdy.add_att_at_t(path[i,1], attract_mu, Dubin_path_step * i / 0.2)

##
qdx.add_upper_pos_constraint_1_by_1(qdx.Tend, staging_area_frontal[0])

qdy.add_lower_pos_constraint_1_by_1(qdx.Tend, staging_area_frontal[1])
qdy.add_upper_pos_constraint_1_by_1(qdx.Tend, staging_area_frontal[2])

qdz.add_lower_pos_constraint_1_by_1(qdx.Tend, staging_area_frontal[3])
qdz.add_upper_pos_constraint_1_by_1(qdx.Tend, staging_area_frontal[4])

exert_con_on_vel_symmetry(qdz, 50, 0.2)
x_traj = qdx.opt()
y_traj = qdy.opt()
z_traj = qdz.opt()

print('x_traj:',x_traj)
print('y_traj:',y_traj)
print('z_traj:',z_traj)

data = {'Dubin_path_step':Dubin_path_step, 
        'Dubin_path': path, 
        'x_poly_res': x_traj, 
        'y_poly_res': y_traj, 
        'z_poly_res': z_traj,
        'T':T_end}

print('Path length:', len(path))
scipy.io.savemat('ReturnTest.mat', data) #_without_ob