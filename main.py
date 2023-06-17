from femsolver import FEMSolver

h_values = [0.2, 0.1, 0.05, 0.02, 0.01]
for h in h_values:
    print(f'>>> Current h = {h}, solver = gauss')
    solver = FEMSolver(h, solver='gauss')
    solver.solve()
    solver.plot()
    
    print(f'>>> Current h = {h}, solver = jacobi')
    solver = FEMSolver(h, solver='jacobi')
    solver.solve()
    solver.plot()