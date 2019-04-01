# Optimization algorithms for two-dimensional and n-dimensional convex / nonconvex functions

This a data repository accompanying the paper "A gradient type algorithm with backward inertial steps for a nonconvex minimization
" by Cristian Daniel Alecsa, Szilárd Csaba László, Adrian Viorel.

Some optimization algorithms are implemented in Matlab, along with graphical representations of the decay in the objective function.

Objective mappings : 
- Beale
- Rastrigin
- Rosenbrock
- Logistic regression with L2 regularization

Optimization algorithms that are included in this repository :
- Polyak with non constant coefficients ( with the same / or modified coefficient appearing in Nesterov's algorithm)
- Nesterov
- Laszlo ( the algorithm form "Convergence rates for an inertial algorithm of gradient type associated to a smooth nonconvex minimization" by Szilárd Csaba László)
- ALV (the algorithm proposed in the present paper)

Parameters :
- tolerance (stopping criteria for the difference in the objective value of the iterates)
- initial conditions
- the objective / loss function (in the case of the logistic regression there are additional parameters)


One can modify the inertial coefficients of these algorithms inside their numerical implementation.
