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

qd = QPP (7, 1)
qd.add_pos_constraint_1_by_1(0.,0.)
# qd.add_pos_constraint_1_by_1(0.2,-1.)
# qd.add_pos_constraint_1_by_1(0.5,4.)
qd.add_pos_constraint_1_by_1(1.,2.)
qd.add_vel_constraint_1_by_1(0.,0.)
qd.add_vel_constraint_1_by_1(1.,0.)
# qd.add_vel_constraint_1_by_1(0.5,0.)
qd.set_coeff(0,0,0,1)

print('qd.Q: ',qd.Q)
print('qd.A: ',qd.A)
print('qd.b: ',qd.b)

x = qd.opt()
print(list(x))

p = np.poly1d(list(x))
x = np.linspace(0, 1, 100)
y = p(x)

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('p(x)')
plt.title('Plot of p(x)')
plt.show()