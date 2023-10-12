from pyomo.environ import*

Model = ConcreteModel()

#Parameters

Tp = 3

d = [0 , 40000 , 48000 , 44000]
# 50 packs in carton and each pack 50 gram

Capr = 6*24000
# 6 is the number of Machines and 24000 is the number of Minutes available in regular Time

Capo = 6*9000
# 6 is the number of Machines and 9000 is the number of Minutes available in regular Time
# 9000 = 60*6*25

I0 = 1000

Co = [0,50000,52500,55000]

h = [0,10000,10500,11000] 

b = 4

#Index

Model.p = RangeSet(Tp)

#Variables

Model.X = Var(Model.p , within = NonNegativeIntegers)

Model.O = Var(Model.p , within = NonNegativeIntegers)

Model.I = Var(Model.p , within = NonNegativeIntegers)

#Objective Function

Model.objecive = Objective(expr = sum(Co[p]*Model.O[p] + h[p]*Model.I[p]
                                      for p in Model.p) , sense = minimize)

#Constraints

Model.St = ConstraintList()

Model.St.add(I0 + Model.X[1] - Model.I[1] == d[1])

for p in range(2,4):
    Model.St.add(Model.I[p-1] + Model.X[p] - Model.I[p] == d[p])

for p in Model.p:
    Model.St.add(b*Model.X[p] - Model.O[p] <= Capr)
    Model.St.add(Model.O[p] <= Capo)

Solver = SolverFactory("cplex")

Solver.solve(Model)

display(Model)