# from cvxopt import matrix, solvers
# import numpy as np
# Q = 2*matrix([ [2, .5], [.5, 1] ])
# p = matrix([1.0, 1.0])
# G = matrix([[-1.0,0.0],[0.0,-1.0]])
# h = matrix([0.0,0.0])
# A = matrix(np.mat([1.0, 1.0]))
# b = matrix(1.0)
# sol=solvers.qp(Q, p, None, None, A, b)

# print(sol['x'])

# import numpy
# from cvxopt import matrix
# P = matrix(numpy.diag([1,0]), tc='d')
# q = matrix(numpy.array([3,4]), tc='d')
# G = matrix([])
# h = matrix([])
# sol = solvers.qp(P,q,G,h)

# print(sol['x'])
import numpy as np
from quad_program import quad_program_for_perching as QPP
import matplotlib.pyplot as plt
import scipy.io

X_start = -1.5
Y_start = [-0.5, 0., 0.5]
Z_start = [-0.5, 0., 0.5]

# terminal position is always set to [0,0,0]


Coef = [1, 0, 0, 0]
Exp_name = 'Coef_'+ str(Coef[0])+\
               '_'+ str(Coef[1])+\
               '_'+ str(Coef[2])+\
               '_'+ str(Coef[3])+'_OptRes'

Traj_name= [['Traj1', 'Traj2', 'Traj3'],
            ['Traj4', 'Traj5', 'Traj6'],
            ['Traj7', 'Traj8', 'Traj9']]

T_terminal = 5

for i in range(3):
    for j in range(3):
        print('Optimizing '+Traj_name[i][j]+'!!! \n')
        qdx = QPP (7, 1)
        qdy = QPP (7, 1)
        qdz = QPP (7, 1)

        qdx.add_pos_constraint_1_by_1(0.,X_start)
        qdx.add_vel_constraint_1_by_1(0.,0.2)
        qdx.add_pos_constraint_1_by_1(T_terminal,0.)
        qdx.add_vel_constraint_1_by_1(T_terminal,0.2)
        qdx.add_acc_constraint_1_by_1(T_terminal,0)

        qdy.add_pos_constraint_1_by_1(0.,Y_start[i])
        qdy.add_vel_constraint_1_by_1(0.,0.)
        qdy.add_pos_constraint_1_by_1(T_terminal,0.)
        qdy.add_vel_constraint_1_by_1(T_terminal,0.)
        qdy.add_acc_constraint_1_by_1(T_terminal,0)

        qdz.add_pos_constraint_1_by_1(0.,Z_start[j])
        qdz.add_vel_constraint_1_by_1(0.,0.)
        qdz.add_pos_constraint_1_by_1(T_terminal,0.)
        qdz.add_vel_constraint_1_by_1(T_terminal,0.)
        qdz.add_acc_constraint_1_by_1(T_terminal,0)
        
        qdx.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])
        qdy.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])
        qdz.set_coeff(Coef[0],Coef[1],Coef[2],Coef[3])

        x_traj = qdx.opt()
        y_traj = qdy.opt()
        z_traj = qdz.opt()

        data = {'x_poly_res': x_traj, 'y_poly_res': y_traj, 'z_poly_res': z_traj,'T':T_terminal}
        scipy.io.savemat(Exp_name + Traj_name[i][j] + '.mat', data)
        print('Saved ' + Exp_name + Traj_name[i][j] + '.mat')

# qd = QPP (7, 1)




# qd.add_pos_constraint_1_by_1(0.,0.)
# # qd.add_pos_constraint_1_by_1(0.2,-1.)
# # qd.add_pos_constraint_1_by_1(0.5,4.)
# qd.add_pos_constraint_1_by_1(1.,2.)
# # qd.add_symmetry_vel_constraint_1_by_1(0.1,0.5)
# # qd.add_symmetry_vel_constraint_1_by_1(0.3,0.5)
# # qd.add_symmetry_vel_constraint_1_by_1(0.5,0.5)
# qd.add_symmetry_vel_constraint_1_by_1(0.7,0.5)
# qd.add_symmetry_vel_constraint_1_by_1(0.9,0.5)
# # qd.add_pos_constraint_1_by_1(0.9,1.)
# # qd.add_pos_constraint_1_by_1(0.5,1.)
# # qd.add_pos_constraint_1_by_1(0.3,1.)
# qd.add_vel_constraint_1_by_1(0.,0.)
# qd.add_vel_constraint_1_by_1(1.,0.)
# # qd.add_vel_constraint_1_by_1(0.5,0.)
# qd.set_coeff(1,1,1,1)

# x = qd.opt()
# print(list(x))

# p = np.poly1d(list(x))
# x = np.linspace(0, 1, 100)
# y = p(x)

# plt.plot(x, y)
# plt.xlabel('x')
# plt.ylabel('p(x)')
# plt.title('Plot of p(x)')
# plt.show()