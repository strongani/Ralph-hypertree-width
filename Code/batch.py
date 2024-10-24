import multiprocessing as mp
import time
import subprocess
import sys
from sample_batch import sample_batch
from batch_less_than_six_fhtw import batch_less_than_six_fhtw
from batch_less_than_1000_vertices import batch_less_than_1000_vertices
from batch_more_than_100_vertices import batch_more_than_100_vertices
from batch_more_than_1000_vertices import batch_more_than_1000_vertices
from new_instance_batch import new_instance_batch

import os
import argparse

def run_fhtw(name):
    print(f'[{time.time()}]: Running {name} (index {names.index(name)}) on node {nodeID}.', flush=True)
    subprocess.run(['python', 'fhtw.py', '-n 60', '-t 3600', name], stdout=f, stderr=sys.stderr)
    print(f'[{time.time()}]: Finished {name} (index {names.index(name)}) on node {nodeID}.', flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='batch',
                    description='Do a batch run of fhtw')
    
    parser.add_argument('numNodes', type=int, help='The number of nodes on the compute server')
    parser.add_argument('numCoresPerNode', type=int, help='The number of cores per node on the compute server')
    parser.add_argument('nodeID', type=int, help='The id of this node')

    args = parser.parse_args()

    # make directory
    if not os.path.exists(f'Output/'):
        try:
            os.makedirs(f'Output/')
        except Exception as e:
            print(f'[{time.time()}]: Caught exception {e} for batch run', file=sys.stderr)
            print(f'OS Path Exists: { os.path.exists("Output/") }')
            traceback.print_exc(file=sys.stderr)

    f = open(f'Output/complete_summary.txt', 'a')

    # set up batch
    numCoresPerNode = args.numCoresPerNode
    numNodes = args.numNodes
    nodeID = args.nodeID

    batch = batch_less_than_1000_vertices

    if nodeID == 0:
        start_index = 2053
        end_index = 3026
    elif nodeID == 1:
        start_index = 400
        end_index = 1500
    elif nodeID == 2:
        batch = batch_more_than_1000_vertices
    elif nodeID == 3:
        batch = new_instance_batch
    else:
        assert False, f'nodeID is {nodeID}?'

    if nodeID == 1:
        graph_stats = [batch[graph] for graph in batch]
        graph_stats.sort(key = lambda graph : (graph[1], graph[2], graph[0]))
        names = [graph[0] for graph in graph_stats]

        still_need_to_run = [58, 133, 134, 172, 275, 276, 373, 391, 401, 414, 415, 448, 450, 452, 473, 474, 475, 476, 477, 490, 494, 512, 516, 517, 856, 1337, 1338, 1386, 1493, 1495, 1497, 1499, 1895, 1928, 1956, 2227, 2766, 2947, 2959, 2967, 2969, 2971, 2973, 2975, 2977, 2979, 2981, 2983, 2985, 2987, 2989, 2991, 2993, 2995, 2997, 2999, 3001, 3003, 3005, 3007, 3009, 3011, 3013, 3015, 3017, 3019, 3021, 3023, 3025]
    
        namesThisNodeRuns = [names[i] for i in still_need_to_run]
    elif nodeID == 2:
        graph_stats = [batch[graph] for graph in batch]
        graph_stats.sort(key = lambda graph : (graph[1], graph[2], graph[0]))
        names = [graph[0] for graph in graph_stats]

        still_need_to_run = [16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28]
        
        namesThisNodeRuns = [names[i] for i in still_need_to_run]
    elif nodeID != 3:
        graph_stats = [batch[graph] for graph in batch]
        graph_stats.sort(key = lambda graph : (graph[1], graph[2], graph[0]))
        names = [graph[0] for graph in graph_stats]

        namesThisNodeRuns = [names[i] for i in range(start_index + nodeID, min(end_index, len(names)), numNodes)]
    else:
        names = list(new_instance_batch)
        names = ['Kakuro-hard-015-ext.xml.hg', 'Kakuro-medium-025-ext.xml.hg', 'Kakuro-medium-018-ext.xml.hg', 'Kakuro-easy-015-ext.xml.hg', 'Kakuro-hard-035-ext.xml.hg', 'rand_q0240.hg', 'Nonogram-128-table.xml.hg', 'rand_q0284.hg', 'Kakuro-medium-015-ext.xml.hg', 'Kakuro-easy-018-ext.xml.hg', 'Kakuro-hard-036-ext.xml.hg', 'Kakuro-easy-025-ext.xml.hg', 'Nonogram-105-table.xml.hg', 'dubois27.hg', 'Kakuro-hard-043-ext.xml.hg', 'Kakuro-hard-044-ext.xml.hg', 'rand_q0478.hg']
        
        namesThisNodeRuns = names

    # sleep
    time.sleep(2)

    print(f'Starting node {nodeID}\n', file=f, flush=True)

    start_time = time.time()

    with mp.Pool(numCoresPerNode) as p:
        p.map(run_fhtw, namesThisNodeRuns, 1)

    end_time = time.time()

    print(f'Finished node {nodeID}, total time taken: {end_time - start_time}\n', file=f, flush=True)