from pyomo.environ import*
import matplotlib.pyplot as plt

Model = ConcreteModel()

#Parameters

Tp = 3

d = [40000 , 48000 , 44000]
# 50 packs in carton and each pack 50 gram

Capr = 6*24000
# 6 is the number of Machines and 24000 is the number of Minutes available in regular Time

Capo = 6*9000
# 6 is the number of Machines and 9000 is the number of Minutes available in regular Time
# 9000 = 60*6*25

I0 = 1000

Co = [50000,52500,55000]

h = [10000,10500,11000] 

b = 4

#Index

Model.p = RangeSet(Tp)

Model.s = RangeSet(0,Tp-1)

Ds = {}
Dh = {}
for i in range(len(d)):
    Ds[i] = d[i]
    Dh[i] = h[i]

# Parameters

Model.D = Param(Model.s , initialize = Ds , mutable = True)

Model.H = Param(Model.s , initialize = Dh , mutable = True)
#Variables

Model.X = Var(Model.p , within = NonNegativeIntegers)

Model.O = Var(Model.p , within = NonNegativeIntegers)

Model.I = Var(Model.p , within = NonNegativeIntegers)

#Objective Function

Model.objective = Objective(expr = sum(Co[p-1]*Model.O[p] + Model.H[p-1]*Model.I[p]
                                      for p in Model.p) , sense = minimize)

#Constraints

Model.St = ConstraintList()

Model.St.add(I0 + Model.X[1] - Model.I[1] == Model.D[0])

for p in range(2,4):
    Model.St.add(Model.I[p-1] + Model.X[p] - Model.I[p] == Model.D[p-1])

for p in Model.p:
    Model.St.add(b*Model.X[p] - Model.O[p] <= Capr)
    Model.St.add(Model.O[p] <= Capo)

X = []

Y = []

for i in range(10):
    Model.D[2] = Model.D[2] + 500
    Solver = SolverFactory("cplex")
    Solver.solve(Model)
    X.append(value(Model.D[2]))
    Y.append(value(Model.objective)/pow(10,9))
    print(value(Model.D[2]), ':',value(Model.objective))
    print("---------------------")
    
plt.plot(X , Y)
plt.title("Sensitivity Analysis on demand of Period 3")
plt.xlabel("Demand")
plt.ylabel("Cost in Billions")
#Model.pprint()