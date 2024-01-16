from scipy.optimize import minimize

# Define the objective function
def objective(x):
    return (x[0] - 1)**2 + (x[1] - 2.5)**2

# Define the constraint function
def constraint(x):
    # if ( (x[0] + x[1] - 3) > 0):
    #     return -1
    # else:
    #     return 1;
    return x[0] + x[1] - 3

# constraint belike: constraint(x) >= 0
# Define the initial guess for the design variables
x0 = [0, 0]

# Define the bounds for the design variables
bounds = [(0, 10), (0, 10)]

# Define the constraint dictionary
constraint_dict = {'type': 'ineq', 'fun': constraint}

# Solve the optimization problem with constraints
result = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=constraint_dict)

# Retrieve the optimized design variables and objective value
x_opt = result.x
obj_value = result.fun

print("Optimized design variables:", x_opt)
print("Optimized objective value:", obj_value)