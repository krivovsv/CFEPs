# The cut-based free energy profiles - CFEPs.

In this series of notebooks we will describe the cut-based free energy profiles (CFEP) and illustrate their usage in analysis of free energy landscapes (FEL), diffusive dynamics, and optimality of reaction coordinates (RC). In particular we will:
 - Introduce the family of ZC,α(x,Δt) profiles and describe how they can be computed from reaction coordinate timeseries r(iΔt).
 - Show how these profiles can be used to compute the conventional, histogram-based, free energy profile $F_H(x)$ and position dependent diffusion coefficient D(x).
 - Show how the profiles can be used to compute various properties of kinetics, e.g. the equilibrium flux, the mean first passage times, and the mean transition path times.
 - Show how these profiles, namely ZC,1(x,Δt) can be used to test the optimality of a putative reaction coordinate, i.e., how close this coordinate is to the committor, and which regions of the coordinate are most suboptimal.
 - Show how a combination of these profiles, namely the θ(x,Δt) function, can be used to check how close is the putative reaction coordinate to an eigenvector of the transfer operator, and what regions of the coordinate are most suboptimal.

The **RC_FEP** notebook provides initial introduction in to using reaction corodinates or indeces for the analysis/description of complex stochastic dynamcis.
