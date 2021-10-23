This repository contains all code needed to run the simulation models in "Autopolyploid establishment depends on life history strategy and the mating outcomes of clonal architecture" by W.E. Van Drunen and J. Friedman. bioRxiv [link](https://doi.org/10.1101/2021.10.21.465190)  

Description:

- **run.R** contains instructions and code for running the simulations in parallel on a personal computer, by calling on **main.R**. Simulations can be run in batches over a set of parameter combinations, or one at a time per parameter set. Basic code for visualizing simulation outcomes is provided. Place all four r-script files supplied in the same directory before running.
- **main.R** is the hub of the simulations. It initializes the population, manages calls to **reproduction.R** and **recruitment.R** to perform these tasks each generation, and keeps track of all data outputs.
- **reproduction.R** performs the reproduction steps in the model, including pollen dispersal, ovule fertilization, and the production of sexual/clonal offspring.
- **recruitment.R** performs the survival and recruitment steps in the model, including death of individuals, offspring dispersal, and offspring recruitment.
