import sys
import random
import time

import settings
import statistics
import turboShrinkComponent

from fractionalBalancedSeparator import *
from hypergraph import *
from fractionalCover import *


# Round LP
#
# gamma is the dict vset --> [0,1] to separate.
# X is dict V x V --> [0,1] = cost to pass through an edge. 
# Y is dict E --> [0,1] = coverage coefficient of edge. 
#     
# while gamma-heaviest component C of G-S has > 1/2 of gamma: 
#       grow S to shrink C.

def roundSeparatorLP(H, gamma, X, Y, *args):
    statistics.totalTimeRoundingBalSepLP -= time.time()

    sumGamma = sum(gamma[e] for e in H.E)
    
    if len(args) > 0:
        careAboutRootCover = args[0]
    else:
        careAboutRootCover = False

    Z = coveredVertices(H, gamma)
    assert len(Z - H.V) == 0, "Z contains vertices outside H.V!"
    S = set()
    C = set(H.V)
    largeComponentZCover = sum(gamma[e] for e in H.E)
    largeComponentNewBagCover = float('inf')
   
    roundingLoop = 0
    while largeComponentZCover * 2 > sumGamma and (not careAboutRootCover or largeComponentNewBagCover >= sumGamma):
        S, C, largeComponentZCover = shrinkLargeComponent(H.inducedSubhypergraph(S | C | Z), gamma, X, Y, S, C, Z, careAboutRootCover)

        if settings.useZCoverEstimator: largeComponentZCover, _ = fractionalCover(H, Z & C)
        largeComponentNewBagCover, _ = fractionalCover(H, (Z & C) | H.openNeighborhood(C))
        
        assert len(S - H.V) == 0, "S contains vertices outside H.V!"
        assert len(C - H.V) == 0, "C contains vertices outside H.V!"
        assert len(Z - H.V) == 0, "Z contains vertices outside H.V!"

        roundingLoop += 1
    
    statistics.maxRoundingLoops = max(statistics.maxRoundingLoops, roundingLoop)
    statistics.totalTimeRoundingBalSepLP += time.time()
    return S
    
    
#computes a list of candidate separator sets to try to remove from current C_hat.
#one candidate is the y-heaviest edge. The others are "X-BFS"-layers from the gamma-heaviest edge. 
def computeCandidateSets(H, gamma, X, Y):
    heaviestYEdge = largestElementOf(H.E, Y)
    ans = [H.verticesOf[heaviestYEdge]]
    
    assert len(H.verticesOf[heaviestYEdge]) > 0, "Heaviest Y-edge is empty!"

    # do the add remove sequence from the heaviest gamma edge. Can try other options here!!
    if settings.candidateSetRoot == "heaviestGammaEdge":
        heaviestGammaEdge = largestElementOf(H.E, gamma)
        ans.extend(computeCandidateSetsFromEdge(H, X, heaviestGammaEdge))
    elif settings.candidateSetRoot == "allNonzeroGammaEdge":
        for e in H.E:
            if gamma[e] < settings.epsilon: continue
            ans.extend(computeCandidateSetsFromEdge(H, X, e))
    elif settings.candidateSetRoot == "randomNonzeroGammaEdge":
        randomGammaEdge = selectRandomElementOutside(H.E, set(), gamma)
        ans.extend(computeCandidateSetsFromEdge(H, X, randomGammaEdge))
    
    return ans
    

def computeCandidateSetsFromEdge(H, X, sourceEdge):
    ans = []
    vertexOrdering = addRemoveSequence(H, X, sourceEdge)
    
    #print("Add remove sequence: ", end="")
    #self.printVertexSet(vertexOrdering)
    
    previousD = 0
    candS = set()
    for (d, v) in vertexOrdering:
        if len(candS) > 0 and (not settings.pruneCandidateSets or d > previousD + settings.epsilon):
            ans.append(set(candS))   
            previousD = d
        if v in candS: candS.remove(v); 
        else: candS.add(v)
    
    if len(candS) > 0 and (not settings.pruneCandidateSets or d > previousD + settings.epsilon): ans.append(set(candS))   
    return ans

    
        
#V is vertex set, X is X-values and D are D-values per vertex. Sort the endpoints of all intervals [D[v]-X[v], D[v]]
def addRemoveSequence(H, X, startEdge):
    endD = H.dijkstra(X, H.verticesOf[startEdge])
    startD = {v : 0 if v in H.verticesOf[startEdge] else min(endD[u] for u in H.G[v]) - (settings.epsilon*settings.epsilon) for v in H.V}
    pairList = [(endD[v], v) for v in H.V if X[(v,v)] > settings.epsilon] + [(startD[v], v) for v in H.V if X[(v,v)] > settings.epsilon]
    pairList.sort()
    return pairList


#returns the first element of container whose weight is largest. 
def largestElementOf(container, weight):
    ans = None
    for elt in container: 
        if ans is None or weight[elt] > weight[ans]: ans = elt
    return ans



def newLossGainPairIsBetter(newLoss, newGain, oldLoss, oldGain):
    if oldLoss == None or newLoss == None:
        return True
    
    if newLoss < 0:
        return newLoss < oldLoss
    elif oldLoss < 0:
        return False
    
    return newLoss * oldGain < oldLoss * newGain


# for each candidate set S: 
#    compute gamma-heaviest component C of G-(S union S_hat). C is empty if none > gamma/2. 
#    update S to N(C) if C is non-empty
#    loss = fcov(Z cup S) - fcov(Z cup S_hat)
#    gain = decrease in y_e's in C_hat to C. 
# Return S with minimum loss/gain (if loss negative then stop and keep) (if gain is 0 then skip)  
def shrinkLargeComponent(H, gamma, X, Y, S, C, Z, careAboutRootCover):

    #it is important that the sum here is taken over gamma and not H.E since H changes but gamma does not,
    sumGamma = sum(gamma[e] for e in gamma)
    coverVal = getCoverVal(H, gamma)

    assert len(C - H.V) == 0, "C contains vertices outside H!"
    assert len(C) > 0, "C is empty"
    
    if settings.turboRounding:
        heaviestYEdge = largestElementOf(H.incidentHyperedges(C), Y)
        candidateSets = [H.verticesOf[heaviestYEdge] & C]

        #add support for settings selecting type of edge to round from in turbo mode.
        assert settings.candidateSetRoot == "randomNonzeroGammaEdge"
        randomGammaEdge = selectRandomElementOutside(H.incidentHyperedges(C), set(), gamma)
        bestCandidateSet, claimedGain, claimedLoss = turboShrinkComponent.bestCandidateSetFromEdge(H, gamma, X, Y, S, C, Z, randomGammaEdge, careAboutRootCover)
        candidateSets.append( bestCandidateSet )
    else:
        candidateSets = computeCandidateSets(H.inducedSubhypergraph(C), gamma, X, Y)
    
    statistics.maxCandidateSets = max(statistics.maxCandidateSets, len(candidateSets))
    
    assert len(candidateSets) > 0, "No candidate sets!"
    assert all(len(candS) > 0 for candS in candidateSets), "Empty candidate set!"
    
    #Will Z be added to the separator or not? This affects whether we include fcov of Z in the loss. 
    if careAboutRootCover:
        baseLoss, _ = fractionalCover(H, (Z | S))
        baseLossEstimator = sum(Y[e] * max((1-coverVal[v])/X[(v,v)] for v in H.verticesOf[e] if v in S) for e in H.incidentHyperedges(S) if Y[e] > 100 * settings.epsilon)
    else:
        baseLoss, _ = fractionalCover(H, S)
        baseLossEstimator = sum(Y[e] * max(1/X[(v,v)] for v in H.verticesOf[e] if v in S) for e in H.incidentHyperedges(S) if Y[e] > 100 * settings.epsilon)       

    sumY = sum(Y[e] for e in H.incidentHyperedges(C))

    bestGain, bestLoss, bestSeparator, bestComponent, bestComponentZCover = None, None, None, None, None
    
    for i in range(len(candidateSets)):
        candS = candidateSets[i]
        assert len(candS - C) == 0, "candS contains vertices outside of C!"
        
        possibleHeavyComponent, possibleHeavyComponentZCover = findHeaviestComponent(H.inducedSubhypergraph(C - candS), gamma)

        candExtendedS, newHeavyComponent, newComponentZCover = inferredSeparatorAndHeavyComponent(H, sumGamma, S, C, Z, candS, possibleHeavyComponent, possibleHeavyComponentZCover)
        
        #calculate loss. Loss can be negative (in that case great). 
        loss = calculateLoss(H, X, Y, Z, coverVal, baseLoss, baseLossEstimator, candExtendedS, careAboutRootCover)
        
        # Calculate gain. Gain is non-negative, can be 0.
        gain = calculateGain(H,Y,sumY,newHeavyComponent)
        
        #update gain/loss thing to fix this.
        #this seems to be set off by floating point errors relatively often, so decreasing sensitivity here by 100. 
        #if settings.turboRounding and i == 1:
        #        assert abs(loss-claimedLoss) < settings.epsilon*100, "Claimed loss not equal to computed loss: (%.3f, %.3f)" % (loss, claimedLoss)
        #        assert abs(gain-claimedGain) < settings.epsilon*100, "Claimed gain not equal to computed gain: (%.3f, %.3f)" % (gain, claimedGain)

        #adjust gain to include decrease in number of vertices. 
        vertexGain = (len(C) - len(newHeavyComponent)) / len(C)
        gain = (1 - settings.vertexGainWeight) * gain + (settings.vertexGainWeight * vertexGain)
        
        if newLossGainPairIsBetter(loss, gain, bestLoss, bestGain): 
            bestGain = gain
            bestLoss = loss
            bestSeparator = set(candExtendedS)
            bestComponent = set(newHeavyComponent)
            bestComponentZCover = newComponentZCover

    return bestSeparator, bestComponent, bestComponentZCover
    

def inferredSeparatorAndHeavyComponent(H, sumGamma, S, C, Z, candS, possibleHeavyComponent, possibleHeavyComponentZCover):

    if possibleHeavyComponentZCover * 2 > sumGamma + settings.epsilon and not settings.useZCoverEstimator:
        possibleHeavyComponentZCover, _ = fractionalCover(H, Z & possibleHeavyComponent)
        
    #determine whether the heaviest component is light enough to be done (according to gamma)
    if possibleHeavyComponentZCover * 2 <= sumGamma + settings.epsilon:
        candExtendedS = S | candS
        newHeavyComponent = set() 
        newComponentZCover = 0
    else: 
        # update additive set to N(new heavy component). Can only do if there is a heavy component. 
        candExtendedS = H.openNeighborhood(possibleHeavyComponent)
        newHeavyComponent = possibleHeavyComponent
        newComponentZCover = possibleHeavyComponentZCover
    
    assert len(candExtendedS - H.V) == 0, "candExtendedS contains vertices outside H.V!"
    
    return candExtendedS, newHeavyComponent, newComponentZCover

   
#calculate loss. Loss can be negative (in that case great). 
def calculateLoss(H, X, Y, Z, coverVal, baseLoss, baseLossEstimator, candExtendedS, careAboutRootCover):

    if settings.useLossEstimator:
        if careAboutRootCover:
            loss = sum(Y[e] * max((1-coverVal[v])/X[(v,v)] for v in H.verticesOf[e] if v in candExtendedS) for e in H.incidentHyperedges(candExtendedS) if Y[e] > 100 * settings.epsilon) - baseLossEstimator
        else:
            loss = sum(Y[e] * max(1/X[(v,v)] for v in H.verticesOf[e] if v in candExtendedS) for e in H.incidentHyperedges(candExtendedS) if Y[e] > 100 * settings.epsilon) - baseLossEstimator
    else:
        if careAboutRootCover:
            loss = fractionalCover(H, Z | candExtendedS)[0] - baseLoss
        else:
            loss = fractionalCover(H, candExtendedS)[0] - baseLoss
            
    return loss

def calculateGain(H,Y,sumY,newHeavyComponent):
    gain = (sumY - sum(Y[e] for e in H.incidentHyperedges(newHeavyComponent))) / sumY
    assert gain > 100 * -settings.epsilon, "gain was " + str(gain)
    assert gain <= 1, "gain was more than 1!"
    return gain



def findHeaviestComponent(H, gamma):
    if len(H.V) == 0: return set(), 0
    candComponents = H.connectedComponents()
    candCompGamma = [sum(gamma[e] for e in H.incidentHyperedges(candComponents[i])) for i in range(len(candComponents))]
    maxComponentIndex = largestElementOf(range(len(candCompGamma)), candCompGamma)
    return candComponents[maxComponentIndex], candCompGamma[maxComponentIndex]
    

#use depth=1 to not reduce the balanced separator with respect to the root bag before returning.
def computeBalancedSeparator(H, root, depth, *args):
    _, gamma = fractionalCover(H, root)
    
    if len(args) > 0:
        acceptIfBelow = args[0]
    else:
        acceptIfBelow = 0
        
    if len(args) > 1:
        careAboutRootCover = args[1]
    else:
        careAboutRootCover = False
        
    if settings.printProgress: print("| "*(depth-1) + "+-", end="--", file=sys.stderr)
    
    X, Y, optLP = fractionalBalancedSeparator(H, gamma)
    
    statistics.maxNonZeroY = max(statistics.maxNonZeroY, len([e for e in H.E if Y[e] > 100 * settings.epsilon]))
    
    if settings.printProgress: print("LP:", optLP, file=sys.stderr)

    if settings.useSeparatorLP == "GaifmanNormalILP":
        assert all(X[v,v] < 100 * settings.epsilon or X[(v,v)] > 1 - 100 * settings.epsilon for v in H.V), "X was not integral! " + str(X)
        balSep = {v for v in H.V if X[(v,v)] > 1 - 100 * settings.epsilon}

        if depth > 1: 
            #if not making progress, add one vertex outside of S with probability proportional to its X(v,v) value.
            # if balSep.issubset(root): balSep.add(selectRandomElementOutside(H.V, root, {v : X[(v,v)] for v in H.V}))
            assert not balSep.issubset(root)
            balSep = reduceBalSepWrtRootBag(H, balSep, root, gamma)

    else:
    
        balSep = None
        balSepVal = float('inf')
        for attemptNr in range(settings.roundingAttempts):
            if balSepVal <= acceptIfBelow: continue
            
            newBalSep = roundSeparatorLP(H, gamma, X, Y, careAboutRootCover) 

            #if not making progress, add one vertex outside of S with probability proportional to its X(v,v) value.
            if depth > 1: 
                if newBalSep.issubset(root): newBalSep.add(selectRandomElementOutside(H.V, root, {v : X[(v,v)] for v in H.V}))
                newBalSep = reduceBalSepWrtRootBag(H, newBalSep, root, gamma)

            #newBalSepVal, _ = fractionalCover(H, newBalSep)
            newBalSepVal, _ = fractionalCover(H, root | newBalSep)

            if newBalSepVal < balSepVal:
                balSep = newBalSep
                balSepVal = newBalSepVal
        
    assert len(balSep - H.V) == 0, "Balanced separator contains vertices outside H.V!"

    return balSep, optLP


def selectRandomElementOutside(container, setToAvoid, weight):
    sumOutsideX = sum(weight[v] if v not in setToAvoid else 0 for v in container)
    assert sumOutsideX > 0, "weight has no mass outside of set to avoid!"
    rValue = random.uniform(0, sumOutsideX)
    for v in container:
        if v in setToAvoid: continue
        rValue -= weight[v]
        if rValue <= 0:
            return v
    assert False, "Failed to select vertex"
    
def reduceBalSepWrtRootBag(H, S, root, gamma):
    CC = H.inducedSubhypergraph(H.V - (S|root)).connectedComponents()
    gammaOfComponent = {i : sum(gamma[e] for e in H.incidentHyperedges(CC[i])) for i in range(len(CC))}
    componentIndexSorted = sorted(range(len(CC)), key=lambda i : gammaOfComponent[i], reverse=True)
    

    sumGamma = sum(gammaOfComponent[i] for i in gammaOfComponent)
    totGamma = 0
    newS = set()
    
    for i in componentIndexSorted:
        newS.update(H.openNeighborhood(CC[i]) & S)
        totGamma += gammaOfComponent[i]
        if totGamma >= sumGamma / 2: break

    assert newS.issubset(S), "newS is not subset of S during trimming!"
    
    #if broke progress on graph, restore some.
    if len(newS - root) == 0 and len(S - root) > 0: newS.add(random.choice(list(S-root)))

    #if fractionalCover(H, root | newS)[0] < fractionalCover(H, root | S)[0]: print("REDUCED!")
    return newS    