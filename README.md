# Screening Template  
Skeleton code for systematically initializing, equilibrating, and analyzing systems  
1. `initialize.py`
    * Use mBuild to construct systems based on composition and parameters specified
    * Briefly equilibrate via EM, NVT, NPT
2. `rwmd.py`
    * Perform RWMD 
3. `production.py`
    * NPT production
4. `compute*.py`
    * Analyze each simulation
5. `process*.py`
    * Gather statistics for each composition from its constituent simulations
6. `plot*.py`
    * Plot properties for each composition with appropriate errors
