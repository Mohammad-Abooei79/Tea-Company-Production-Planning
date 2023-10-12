from pyomo.environ import*

Model = ConcreteModel()

#Parameters

T = 12

Tp = 3

k = 5

d = [[0 , 1000 , 1000 , 1000 , 1000 , 1200 , 1200 , 1200 , 1200 , 1100 , 1100 , 1100 , 1100],
     [0 , 1000 , 1000 , 1000 , 1000 , 1200 , 1200 , 1200 , 1200 , 1100 , 1100 , 1100 , 1100],
     [0 , 2000 , 2000 , 2000 , 2000 , 2400 , 2400 , 2400 , 2400 , 2200 , 2200 , 2200 , 2200],
     [0 , 4500 , 4500 , 4500 , 4500 , 5400 , 5400 , 5400 , 5400 , 4950 , 4950 , 4950 , 4950],
     [0 , 1500 , 1500 , 1500 , 1500 , 1800 , 1800 , 1800 , 1800 , 1650 , 1650 , 1650 , 1650]]


Xstar = [0, 49500 , 37500 , 44000]

Ostar = [0 , 54000 , 6000 , 32000]

Istar = [0 , 10500 , 0 , 0]

Capr = 6*24000

Cap_dummy = [(Capr + i)/4 for i in Ostar]
Cap = []
for i in Cap_dummy: Cap.extend([i,i,i,i])
Cap[0] = 0
for i in range(3):
    Cap.remove(36000)

Capo = 6*2250

I0 = [0 , 150 , 200 , 200 ,  250 , 200]

Cr = [0,160000,160000,160000,160000,168000,168000,168000,168000,177600,177600,177600,177600]

h = [0,2500,2500,2500,2500,2625,2625,2625,2625,2750,2750,2750,2750]

b = 4

miu = 0.2

#Index

Model.k = RangeSet(k)

Model.t = RangeSet(T)

#Variables

Model.Y = Var(Model.k , Model.t , within = NonNegativeIntegers)

Model.I = Var(Model.k , Model.t , within = NonNegativeIntegers)

#Objective Function

Model.objective = Objective(expr = sum(Cr[t]*Model.Y[k,t] + h[t]*Model.I[k,t] 
                                      for k in Model.k for t in Model.t) , sense = minimize)

#Constraints

Model.St = ConstraintList()

for k in Model.k:
    Model.St.add(I0[k] + Model.Y[k,1] - Model.I[k,1] == d[k-1][1])

for k in Model.k:
    for t in range(2,T+1):
        Model.St.add(Model.I[k,t-1] + Model.Y[k,t] - Model.I[k,t] == d[k-1][t])

for t in Model.t:
    Model.St.add(sum(b*Model.Y[k,t] for k in Model.k) <= Cap[t])

for t in Model.t:
    Model.St.add(sum(b*Model.Y[k,t] for k in Model.k) >= (1 - miu)*Cap[t])

for T in range(1,Tp+1):
    Model.St.add(sum(Model.Y[k,t] for k in Model.k for t in range(4*(T-1)+1 , 4*T+1)) == Xstar[T])

for T in range(1,Tp+1):
    Model.St.add(sum(Model.I[k,4*T] for k in Model.k) == Istar[T])

Solver = SolverFactory("cplex")

Solver.solve(Model)

display(Model)