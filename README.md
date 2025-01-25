The code, hypergraph instances, and our compiled outputs for Ralph (Randomized Approximation using Linear Programming for Hypertree-Decompositions)

For running the code:
- the requirements for packages are listed in _requirements.txt_
- besides that, you either need _Gurobi_ installed/with a license, or change _settings.py_ to use a different solver (many available in _PuLP_)
- many other hyperparameters can be set in _settings.py_
- The pilot file to run is _fhtw.py_, the hypergraph ran needs to be in _Hypergraphs/_. For batch runs, look at _batch.py_ (very specificized to our HPC)


For the compiled results:
- We have made our raw output available, as well as compiled data from our raw output into multiple CSV files
- For each entry in any CSV file with some algorithm's upper/lower bound on an instance, a positive real value corresponds to the upper/lower bound found within the time limit, a -1 or 'inf' corresponds to a timeout, and no entry corresponds to no data/not attempted.
