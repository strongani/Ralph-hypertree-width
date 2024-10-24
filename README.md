The code, hypergraph instances, and our compiled outputs for Ralph (Randomized Approximation using Linear Programming for Hypertree-Decompositions)

For running the code:
- the requirements for packages are listed in _requirements.txt_
- besides that, you either need _Gurobi_ installed/with a license, or use change _settings.py_ to use a different solver (many available in _PuLP_)
- The pilot file to run is _fhtw.py_, the hypergraph ran needs to be in _Hypergraphs/_. For batch runs, look at _batch.py_ (very specificized to our HPC)
