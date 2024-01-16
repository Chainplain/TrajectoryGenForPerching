from scipy.optimize import minimize

def objective(x):
    return x[0]**2 + x[1]**2

def constraint1(x):
    return x[0] + x[1] - 1



x0 = [0, 0]
cons = [{'type': 'eq', 'fun': constraint1}]
result = minimize(objective, x0, constraints=cons)
print(result.x)