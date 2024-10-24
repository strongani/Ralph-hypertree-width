import pulp
import time

from hypergraph import *
from fractionalCover import *

import settings
import statistics



# Make the Balanced Separator LP (Incidence Graph Style)!! 
#
# H is the hypergraph.
# gamma is a dict : H.E --> [0,1]. Represents a fractional cover for the vertex set Z to be separated. 
#
def makeBalancedSeparatorHypergraphNormalLP(H, gamma):
    sum_gamma = sum(gamma[e] for e in H.E)

    # Important edges: F = all edges e in eset with gamma[e] > 0
    F = {e for e in H.E if gamma[e] > 0} 
    
    prob = pulp.LpProblem("Balanced_Clique_Covered_Separator", pulp.LpMinimize)

    # Variables: y_e (one per edge e in E), x_v (one per vertex v in V), d_e,v (one per e in F), d_e,e' (one per e in F and e' in E)
    Y = pulp.LpVariable.dicts("y", list(H.E), 0, 1)
    X = pulp.LpVariable.dicts("x", list(H.V), 0, None)
    D_ev = pulp.LpVariable.dicts("DFV", [(e,v) for e in F for v in H.V], 0, 1)
    D_ee = pulp.LpVariable.dicts("DFE", [(e1,e2) for e1 in F for e2 in H.E],0, 1)
    
    # Minimize sum of y_e's
    prob += (pulp.lpSum([Y[e] for e in H.E]), "Total_Covering_Cost")
    
    # x_v == sum of y_e over all e such that v is in e. 
    for v in H.V: 
        prob += pulp.lpSum([Y[e] for e in H.edgesOf[v]]) == X[v]
        
    #must contain at least one vertex not in the set to be split.
    #only true if G-Z is connected and every vertex in Z has a neighbor outside. 
    Z = coveredVertices(H, gamma)
    if len(H.inducedSubhypergraph(H.V - Z).connectedComponents()) == 1 and all(len(H.G[v] - Z) > 0 for v in Z):
        prob += pulp.lpSum([X[v] for v in H.V - Z]) >= 1
    
    # for every edge e in F, e' in E, vertex v in e': d_e,v <= [d_e,e' if e != e'] + x_v
    # for every edge e in F, e' in E, vertex v in e': d_e,e' <= d_e,v
    # So we pay x_v on the way from the edge to the vertex.
    for f in F:
        for e in H.E:
            for v in H.verticesOf[e]:
                if f != e:
                    prob += D_ev[(f,v)] <= D_ee[(f,e)] + X[v]
                else: 
                    prob += D_ev[(f,v)] <= X[v]
                prob += D_ee[(f,e)] <= D_ev[(f,v)]
    
    #from edge to edge, symmetric
    for f1 in F:
        for f2 in F:
            prob += D_ee[(f1, f2)] == D_ee[(f2, f1)]
    
    # for every edge e in F: sum over e' in F d_e,e' * gamma_e' >= Gamma / 2
    for f in F: 
        prob += pulp.lpSum([gamma[e] * D_ee[(f,e)] for e in F]) >= (sum_gamma / 2) - settings.epsilon
    
    return prob, X, Y


# Make the Balanced Separator LP (THE NORMAL GAIFMAN VERSION)!! 
#
# Also allows for X to be set to integral.
#
# H is the relevant hypergraph.
# gamma is a dict : H.E --> [0,1]. Represents a fractional cover for the vertex set Z to be separated. 
#
def makeBalancedSeparatorGaifmanNormalLP(H, gamma):
    if H.G is None: H.G = H.gaifmanGraph()

    sum_gamma = sum(gamma[e] for e in H.E)

    # Important edges: F = all edges e in eset with gamma[e] > 0
    F = {e for e in H.E if gamma[e] > 0} 
    
    prob = pulp.LpProblem("Balanced_Clique_Covered_Separator", pulp.LpMinimize)

    # Variables: y_e (one per edge e in E), x_v (one per vertex v in V), d_e,v (one per e in F), d_e,e' (one per e in F and e' in E)
    Y = pulp.LpVariable.dicts("y", list(H.E), 0, 1)
    
    #vertex variables are binary if solving ILP. For small instances only!
    if settings.useSeparatorLP == "GaifmanNormalILP":
        X = pulp.LpVariable.dicts("x", list(H.V), 0, None, cat="Binary")
    else:
        X = pulp.LpVariable.dicts("x", list(H.V), 0, None)

    D_ev = pulp.LpVariable.dicts("DFV", [(e,v) for e in F for v in H.V], 0, 1)
    D_ee = pulp.LpVariable.dicts("DFF", [(e1,e2) for e1 in F for e2 in F],0, 1)
    
    # Minimize sum of y_e's
    prob += (pulp.lpSum([Y[e] for e in H.E]), "Total_Covering_Cost")
    
    # x_v == sum of y_e over all e such that v is in e. 
    for v in H.V: 
        if settings.useSeparatorLP == "GaifmanNormalILP":
            prob += pulp.lpSum([Y[e] for e in H.edgesOf[v]]) >= X[v]
        else:
            prob += pulp.lpSum([Y[e] for e in H.edgesOf[v]]) == X[v]

    #must contain at least one vertex not in the set to be split.
    #only true if G-Z is connected and every vertex in Z has a neighbor outside.
    Z = coveredVertices(H, gamma)
    if len(H.inducedSubhypergraph(H.V - Z).connectedComponents()) == 1 and all(len(H.G[v] - Z) > 0 for v in Z):
        prob += pulp.lpSum([X[v] for v in H.V - Z]) >= 1

    
    #from edge to vertex in edge.
    for f in F:
        for v in H.verticesOf[f]:
            prob += D_ev[(f, v)] <= X[v]
    
    #from edge to vertex, triangle inequality:
    for f in F:
        for v in H.V:
            for u in H.G[v]:
                prob += D_ev[(f,u)] <= D_ev[(f, v)] + X[u]
    
    #from edge to edge, min vertex in edge
    for f1 in F:
        for f2 in F:
            for v in H.verticesOf[f2]:
                prob += D_ee[(f1, f2)] <= D_ev[(f1, v)]
    
    #from edge to edge, symmetric
    for f1 in F:
        for f2 in F:
            prob += D_ee[(f1, f2)] == D_ee[(f2, f1)]

    # for every edge e in F: sum over e' in F d_e,e' * gamma_e' >= Gamma / 2
    for f in F: 
        prob += pulp.lpSum([gamma[e] * D_ee[(f,e)] for e in F]) >= (sum_gamma / 2) - settings.epsilon
    
    return prob, X, Y


# Make the Balanced Separator LP (THE BOOSTED GAIFMAN VERSION)!! 
#
# the boosting only helps when X is fractional, so no point in having an ILP version of this. 
#
# H is the relevant hypergraph.
# gamma is a dict : H.E --> [0,1]. Represents a fractional cover for the vertex set Z to be separated. 
#
def makeBalancedSeparatorGaifmanBoostedLP(H, gamma):
    if H.G is None: H.G = H.gaifmanGraph()

    sum_gamma = sum(gamma[e] for e in H.E)

    # Important edges: F = all edges e in eset with gamma[e] > 0
    F = {e for e in H.E if gamma[e] > 0} 
    
    prob = pulp.LpProblem("Balanced_Clique_Covered_Separator", pulp.LpMinimize)

    # Variables: y_e (one per edge e in E), x_v (one per vertex v in V), d_e,v (one per e in F), d_e,e' (one per e in F and e' in E)
    Y = pulp.LpVariable.dicts("y", list(H.E), 0, 1)
    X = pulp.LpVariable.dicts("x", list((u, v) for u in H.V for v in H.G[u] | {u}), 0, None)
    D_ev = pulp.LpVariable.dicts("DFV", [(e,v) for e in F for v in H.V], 0, 1)
    D_ee = pulp.LpVariable.dicts("DFF", [(e1,e2) for e1 in F for e2 in F],0, 1)
    
    # Minimize sum of y_e's
    prob += (pulp.lpSum([Y[e] for e in H.E]), "Total_Covering_Cost")
    
    #x_(v,v) == sum of y_e over all e such that v is in e. 
    for v in H.V: 
        prob += pulp.lpSum([Y[e] for e in H.edgesOf[v]]) == X[(v, v)]
        
    #must contain at least one vertex not in the set to be split.
    #only true if G-Z is connected and every vertex in Z has a neighbor outside.
    Z = coveredVertices(H, gamma)
    if len(H.inducedSubhypergraph(H.V - Z).connectedComponents()) == 1 and all(len(H.G[v] - Z) > 0 for v in Z):
        prob += pulp.lpSum([X[(v, v)] for v in H.V - Z]) >= 1
        
    #x_(u,v) == sum of y_e over all e such that u is not in e and v is in e
    for v in H.V: 
        for u in H.G[v]:
            prob += pulp.lpSum( [Y[e] for e in H.edgesOf[v] if e not in H.edgesOf[u] ]) == X[(u, v)]
    
    #from edge to vertex in edge.
    for f in F:
        for v in H.verticesOf[f]:
            prob += D_ev[(f, v)] <= X[(v, v)]
    
    #from edge to vertex, triangle inequality:
    for f in F:
        for v in H.V:
            for u in H.G[v]:
                prob += D_ev[(f,u)] <= D_ev[(f, v)] + X[(v, u)]
    
    #from edge to edge, min vertex in edge
    for f1 in F:
        for f2 in F:
            for v in H.verticesOf[f2]:
                prob += D_ee[(f1, f2)] <= D_ev[(f1, v)]
    
    #check that this is ok
    #from edge to edge, symmetric
    for f1 in F:
        for f2 in F:
            prob += D_ee[(f1, f2)] == D_ee[(f2, f1)]

    # for every edge e in F: sum over e' in F d_e,e' * gamma_e' >= Gamma / 2
    for f in F: 
        #prob += pulp.lpSum([gamma[e] * D_ee[(f,e)] for e in F]) >= (sum_gamma / 2) - settings.epsilon
        #is the epsilon creating problems??
        prob += pulp.lpSum([gamma[e] * D_ee[(f,e)] for e in F]) >= (sum_gamma / 2)
    
    return prob, X, Y



    
def fractionalBalancedSeparator(H, gamma):  
    
    statistics.totalTimeSolvingBalSepLP -= time.time()

    if settings.useSeparatorLP == "gammaPlusFractionalCoverHeuristic":
        # 1/2gamma is a feasible solution for the LP. Add uniform-ish spread over all the vertices using fractional cover of all vertices. 
        fCovV, coverCoeff = fractionalCover(H, H.V)
        sumGamma = sum(gamma[e] for e in H.E)
        
        #scale cover coefficients so they add to sumGamma * settings.heuristicFCovWeight
        coverCoeff = {e : coverCoeff[e] * settings.heuristicFCovWeight * sumGamma / fCovV for e in H.E}
        
        # 0.5*gamma[e] is a feasible solution. Add 0.5 * coverCoeff to get uniform(ish) spread on all other vertices. 
        eVal = {e : 0.5 * (gamma[e] + coverCoeff[e]) for e in H.E}
        X_proxy = {v : sum(eVal[e] for e in H.edgesOf[v]) for v in H.V}
        vVal = {(u, v) : X_proxy[v] for v in H.V for u in (H.G[v] | {v})}

        # this LP feasilble solution can not be used as a lower bound on the fhtw, hence 0 for lpVal
        lpVal = 0
    else:
        if settings.useSeparatorLP in {"GaifmanNormal", "GaifmanNormalILP"}:
            prob, X, Y = makeBalancedSeparatorGaifmanNormalLP(H, gamma)
            prob.solve(settings.getSolver())    
            vVal = {(u, v) : X[v].varValue for v in H.V for u in (H.G[v] | {v})}    
        elif settings.useSeparatorLP == "HypergraphNormal":
            prob, X, Y = makeBalancedSeparatorHypergraphNormalLP(H, gamma)
            prob.solve(settings.getSolver())
            vVal = {(u, v) : X[v].varValue for v in H.V for u in (H.G[v] | {v})}
        elif settings.useSeparatorLP == "GaifmanBoosted":
            prob, XX, Y = makeBalancedSeparatorGaifmanBoostedLP(H, gamma)
            prob.solve(settings.getSolver())
            vVal = {(u, v) : XX[(u, v)].varValue for u in H.V for v in (H.G[u] | {u})}
            
        eVal = {e : Y[e].varValue for e in H.E}
        lpVal = sum(eVal[e] for e in H.E)

    statistics.totalTimeSolvingBalSepLP += time.time()
    
    return vVal, eVal, lpVal
    
