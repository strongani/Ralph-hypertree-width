### This repository contains the code, instances, and compiled outputs for https://dl.acm.org/doi/10.1145/3725296

We present our algorithm _Ralph_ (_**R**andomized **A**pproximation using **L**inear **P**rogramming for **H**ypertree-Decompositions_)

**Minimum requirements for running the code**
- Python 3.12
- 2 GB RAM

**Quick Summary**:
- the requirements for packages are listed in [Code/requirements.txt](Code/requirements.txt)
- besides that, you either need _Gurobi_ installed/with a license, or change [Code/settings.py](Code/settings.py) to use a different solver (many available in _PuLP_)
- many other hyperparameters can be set in [Code/settings.py](Code/settings.py)
- The pilot file to run is [Code/fhtw.py](Code/fhtw.py), the hypergraph ran needs to be in [Hypergraphs/](Hypergraphs/). For batch runs, look at [Code/batch.py](Code/batch.py)(very specificized to our HPC)

**Detailed steps:**
- Download / clone the repository
- In the repository, install the python requirements in [Code/requirements.txt](Code/requirements.txt)
  - This can be done in many ways, an easy way is to use a virtual environment
  - Here is an example using virtualenv:
    - `cd {repository}/Code`
    - `pip install virtualenv`
    - `virtualenv -p 3.12 venv`
    - `source venv/bin/activate` (or `.\venv\Scripts\activate.bat` on Windows)
    - `pip install -r requirements.txt`
- Add a Gurobi license or edit the settings file
  - The default settings use the Gurobi solver, which requires a license
  - If you do not have a Gurobi license, you can change this in [Code/settings.py](Code/settings.py)
    - Specifically, comment line 27 (`solver = pulp.GUROBI(msg=False, timeLimit=3600)`) and uncomment line 29 (`solver = pulp.HiGHS(msg=False, timeLimit=3600)`)
- Run the [Code/fhtw.py](Code/fhtw.py) file as follows: `python fhtw.py {hypergraph} -n {numRuns} -t {timeLimit}` 
  - For a quick test, run `python fhtw.py 1.dtl` (hypergraphs can be found in the [Hypergraphs/](Hypergraphs/) folder)
    - This will create an `{repository}/Code/Output/1.dtl` folder, in which there will be 10 output files, and one summary file

**Format of the output files**
- The format of the output file can be found at [Results/output_sample.out](Results/output_sample.out)

**For running multiple hypergraphs**
- As running the entire dataset of hypergraph can take days, we used a HPC setup to run multiple hypergraphs in batch. Our code can be found at [Code/batch.py](Code/batch.py). Parts of it are very specific to our setup, but the main ideas should be transferrable.

**For the compiled results:**
- We have made our raw output available, as well as compiled data from our raw output into multiple CSV files
- For each entry in any CSV file with some algorithm's upper/lower bound on an instance, we use the following convention: a positive real value corresponds to the upper/lower bound found within the time limit, a -1 or 'inf' corresponds to a timeout, and no entry corresponds to no data/not attempted.  
