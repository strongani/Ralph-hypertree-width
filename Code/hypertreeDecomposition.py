import random
import sys
import time

import settings
import statistics
from hypergraph import *
from fractionalCover import *
from balancedSeparator import *

#tree decomposition of hypergraph
class treeDecomposition:
    
    #bag is dict node ID --> set of nodes
    #children is map from ID to set of ID's (the child nodes). Empty childen mean leaves.
    #makes empty tree decomposition
    def __init__(self, rootBag):
        self.bag = {0 : rootBag}
        self.children = {0 : set()}
        self.firstFreeID = 1
        self.fhtwUpper = -1
        self.fhtwLower = -1
        self.ghwUpper = -1
        self.fhtwLowerWitness = None

    #number of nodes in the tree decomposition
    def numNodes(self):
        return len(self.children)

    #adds a new child with a new ID and given bag. Returns ID of new child.
    def add_leaf(self, parentID, newBag):
        childID = self.firstFreeID
        self.firstFreeID += 1
        self.children[parentID].add(childID)
        self.children[childID] = set()
        self.bag[childID] = newBag
        return childID
    
    #if parent is a subset of child or vice versa, contract to parent.
    def trim(self, parent, child): 
        if self.bag[parent].issubset(self.bag[child]) or self.bag[child].issubset(self.bag[parent]):
            self.bag[parent] = self.bag[parent] | self.bag[child]
            self.children[parent] |= self.children[child]
            self.children[parent].remove(child)
            del self.children[child]
            del self.bag[child]
            return True

    #pretty slow but whatever, will be way faster than the other stuff.
    def trimTD(self):
        success = True
        while(success): 
            success = False
            for v in self.children:
                for child in self.children[v]:
                    if self.trim(v, child):
                        success = True
                        break
                if success: break    
                
    def __str__(self, ):
        ans = [str(len(self.children)) + ", " + str(self.fhtwLower) + ", " + str(self.fhtwUpper) + ", " + str(self.ghwUpper)]
        self.DFS_string(0, ans, 0)        
        return '\n'.join(ans)
                
    def DFS_string(self, pos, ans, depth):
        ans.append(("| " * depth) + "+-- , " + str(pos) + " , " + self.children[pos].__str__() + " , " + self.bag[pos].__str__())
        for c in self.children[pos]: self.DFS_string(c, ans, depth + 1)
        return
    
    def verifyDecomposition(self, H):
        for v in H.V:
            bagsOfv = {b for b in self.bag if H.vName[v] in self.bag[b]}
            assert len(bagsOfv) > 0, H.vName[v] + " not found in any bag!"
            
            bagsOfvMinimum = min(bagsOfv)
            reachableFromRoot = set()
            self.visitFromRoot(bagsOfvMinimum, bagsOfv, reachableFromRoot)
            
            assert bagsOfv.issubset(reachableFromRoot), "Set of bags containing " + H.vName[v] + " is not connected! " + str(bagsOfv) + " " + str(reachableFromRoot)
            
        for e in H.E:
            bagsOfe = {b for b in self.bag if all(H.vName[v] in self.bag[b] for v in H.verticesOf[e])}
            assert len(bagsOfe) > 0, H.eName[e] + " not found in any bag!"
    
    
    def visitFromRoot(self, root, allowedBags, ans):
        ans.add(root)
        for v in self.children[root] & allowedBags: self.visitFromRoot(v, allowedBags, ans)
    

def decomposeWithWarmStart(H):
    outputDec = treeDecomposition(set())

    for component in H.connectedComponents():
        H_c = H.inducedSubhypergraph(component)

        numVerticesToAdd = 5
        separate = set()
        maxCC = len(H_c.V)

        numTriesLeft = 20
        
        while maxCC > len(H_c.V) * 0.7 and numTriesLeft > 0:
            separate = separate | {random.choice(list(H_c.V)) for _ in range(numVerticesToAdd)}
            numVerticesToAdd = min(len(H_c.V), int(numVerticesToAdd * 2))
            root, optLP = computeBalancedSeparator(H_c, separate, 1)
            if optLP > outputDec.fhtwLower:
                outputDec.fhtwLower = optLP
                outputDec.fhtwLowerWitness = ({H.vName[v] for v in H_c.V}, {H.vName[v] for v in separate})
            CC = H_c.inducedSubhypergraph(H_c.V - root).connectedComponents()
            maxCC = 0 if len(CC) == 0 else max(len(C) for C in CC)
            numTriesLeft -= 1
            
        decompose(H_c, root, outputDec, 0, 2)   
    outputDec.trimTD()
    
    return outputDec

# somewhat slower, not clear if better?
def decomposeWithWarmStartFcov(H):
    outputDec = treeDecomposition(set())

    for component in H.connectedComponents():
        H_c = H.inducedSubhypergraph(component)
    
        fcovVal, fcovCoeff = fractionalCover(H_c, H_c.V)
        
        numEdgesToAdd = 5
        separate = set()
        maxCC = len(H_c.V)

        numTriesLeft = 15
        
        while maxCC > fcovVal * 0.7 and numTriesLeft > 0:
            for _ in range(numEdgesToAdd): separate = separate | H_c.verticesOf[selectRandomElementOutside(H_c.E, set(), fcovCoeff)]
            numEdgesToAdd = min(len(H_c.E), int(numEdgesToAdd * 2))
            root, optLP = computeBalancedSeparator(H_c, separate, 1)
            if optLP > outputDec.fhtwLower:
                outputDec.fhtwLower = optLP
                outputDec.fhtwLowerWitness = ({H.vName[v] for v in H_c.V}, {H.vName[v] for v in separate})
            CC = H_c.inducedSubhypergraph(H_c.V - root).connectedComponents()
            maxCC = 0 if len(CC) == 0 else max(fractionalCover(H_c, C)[0] for C in CC)
            numTriesLeft -= 1
            
        decompose(H_c, root, outputDec, 0, 2)   
    outputDec.trimTD()
    
    return outputDec




def decompose(H, root, outputDecomposition, decompositionPos, depth):
    

    assert len(root - H.V) == 0, "root contains vertices outside H.V" 
    
    fcovVal, fcovCoeff = fractionalCover(H, root)
    
    #greedy set cover rounding    
    intCovVal, _ = greedyIntegralCover(H.inducedSubhypergraph(root), fcovCoeff)

    outputDecomposition.fhtwUpper = max(outputDecomposition.fhtwUpper, fcovVal)
    outputDecomposition.ghwUpper = max(outputDecomposition.ghwUpper, intCovVal)

    
    if settings.printProgress:
        print("| "*(depth-1) + "+-", end="--", file=sys.stderr)
        print(fcovVal, intCovVal, H.printVertexSet(root), file=sys.stderr)

    if outputDecomposition.fhtwUpper >= settings.killSwitch: return
     
    rootDecompositionPos = outputDecomposition.add_leaf(decompositionPos, {H.vName[v] for v in root})

    if H.V == root: return
            
    CC = (H.inducedSubhypergraph(H.V - root)).connectedComponents()
    if(len(CC) > 1):
        for C in CC:
            newRoot = H.openNeighborhood(C)
            decompose(H.inducedSubhypergraph(C | newRoot), newRoot, outputDecomposition, rootDecompositionPos, depth+1)
        return
        
    C = CC[0]
    newRoot = H.openNeighborhood(C)
    if not newRoot == root: 
        decompose(H.inducedSubhypergraph(C | newRoot), newRoot, outputDecomposition, rootDecompositionPos, depth+1)
        return

    newRoot = increaseToMaxWithSameCover(H, root)
    if not newRoot == root:
        decompose(H, newRoot, outputDecomposition, rootDecompositionPos, depth+1)
        return

    #print("| "*depth+"BALSEP!")
    S, balSepLPValue = computeBalancedSeparator(H, root, depth, outputDecomposition.fhtwUpper - 2, True)
    if balSepLPValue > outputDecomposition.fhtwLower:
        outputDecomposition.fhtwLower = balSepLPValue
        outputDecomposition.fhtwLowerWitness = ({H.vName[v] for v in H.V}, {H.vName[v] for v in root})

    if len(S - root) > 0:

        # not clear this is best done here, or after each individual bal sep is picked. 
        #newRoot = reduceBalSepWrtRootBag(H, S, root, fcovCoeff)

        decompose(H, root | S, outputDecomposition, rootDecompositionPos, depth + 1)
        return
    
    if settings.printProgress: 
        print("| "*(depth-1) + "+-", end="--", file=sys.stderr)
        print("FAIL", file = sys.stderr)

    assert False , "Bal sep failed to make progress!"
    #decompose(H, root | {random.choice(list(H.V - root))}, outputDecomposition, rootDecompositionPos, depth + 1)
    

def main():
    #for debugging, 
    if settings.debug:
        sys.stdin = open('Hypergraphs/s820.hg')   
        #sys.stdin = open('input\s820.hg')    
        #sys.stdin = open('testInputs/c12.hg')    
        #sys.stdin = open('testInputs/danielGraph.hg')    

    #doing this even in not debug mode, rember to remove:)
    #sys.stdin = open('input/s820.hg')    
    # sys.stdin = open('Hypergraphs/s820.hg')   
            
    sys.setrecursionlimit(100000) 

    H = hypergraph()

    #H.printMe()
    #H.printVertexSet(H.V)

    start_time = time.time()

    #outputDec = treeDecomposition(set())
    #decompose(H, {random.choice(list(H.V))}, outputDec, 0, 1)

    outputDec = decomposeWithWarmStart(H)

    print("--- %.3f seconds ---" % (time.time() - start_time), file=sys.stderr)
    if outputDec.fhtwUpper < settings.killSwitch:
        print("%.3f <= fhtw <= %.3f --- ghw <= %d" % (outputDec.fhtwLower , outputDec.fhtwUpper , outputDec.ghwUpper), file=sys.stderr)
    else:
        print("%.3f <= fhtw, killed by kill switch." % (outputDec.fhtwLower), file=sys.stderr)

    statistics.printStatistics()

    if len(sys.argv) >= 2 and outputDec.fhtwUpper < settings.killSwitch: 
        writeToFile = open(sys.argv[1], 'w')
        print(outputDec, file = writeToFile)
        writeToFile.close()

    if outputDec.fhtwUpper < settings.killSwitch:   
        outputDec.verifyDecomposition(H)

if __name__ == '__main__':
    main()


def runDecomposition(H):
    statistics.resetStatistics()

    start_time = time.time()
    outputDec = decomposeWithWarmStart(H)
    
    total_time = time.time() - start_time
    stats = statistics.getStatistics()
    stats['total_time'] = total_time

    return outputDec, stats