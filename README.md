# Unconstrained Binary Models of the Travelling Salesman Problem Variants for Quantum Optimization

Person responsible for data: Özlem Salehi (osalehi[at]iitis.pl).

The scripts necessary for generating the results provided in the "Unconstrained Binary Models of the Travelling Salesman Problem Variants for Quantum Optimization".

## Data used in the publication

Data used in the publication is located inside the data folder. The folder ```instances``` contains the instances. ```dwave_results_ilp``` and ```dwave_results_edge```contain the results of the experiments obtained from D-Wave. For an instance named ```complete_0_3_0```, there are 3 files:

```complete_0_3_0_{p1}_{p2}_qa_{chain}_{time}```

```complete_0_3_0_{p1}_{p2}_qa_{chain}_{time}_sample```

```complete_0_3_0_{p1}_{p2}_qa_{chain}_{time}_energy```

```p1``` and ```p2``` are the penalty parameters, ```chain``` is the chain strength and it is not indicated in the file name if it is 2 and ```time``` is the annealing time. The first file contains the energies of the sampleset obtained, route coresponding to each sample, cost of the route and information about whether the found route is an Hamiltonian cycle, whether it satisfies the time windows whether it is feasible, or optimal. 
The second file is the raw experiment output returned by D-Wave. Third file stores the energy of the bit assignment encoding the optimal route with no penalties, for the penalty parameters ```p1``` and ```p2```.

The folder named ```AFG``` contains the instances from the dataset presented in [1].

## Reproducing data

### Generating random instances  

To generate random instances, run the following command in the main directory:

```
python tsptw_generator.py instances 10 3 
```

 ```instances``` is the name of the folder to store instances. ```10``` is the number of instances to be generated. ```3``` is the number of cities including the depot. The command above uses the default values for the upper and lower bounds of the time windows and cost matrix. The default values are the values which are used to generate the instances used in the publication. 
 
 The generated output is stored in the npz format with the variables n (the number of cities), weight (cost matrix) and window (time-windows).

The details of the optional keywords are described below:

```-l```: Lower bound for the entries in the cost matrix. Default is 1.
                        
```-u```: Upper bound for the entries in the cost matrix. Default is 10.

```-et```: Upper bound for the earliest start time. Default is 20.

```-lt```: Upper bound for the latest start time. Default is 40.


If you want to generate symmetric instances please add the keyword ```-sym```. If you want to store the instances in the txt format as well, then add the keyword ```-txt```.

### Penalty parameter search using simulated annealing

To search for penalty parameters, run the following code in the main directory.

```
python parameter_search.py instances edge_params e 
```
The results of the annealing experiments are stored as Pandas dataframe. In addition, the probability of the samples encoding optimal routes is calculated for different penalty values and the penalty values which maximize the probability for eacn instance is stored as a dictionary. ```instances``` is the name of the folder in which the instances are located. ```edge_params``` is the name of the folder to store the outputs. ```e``` represents the edge-based formulation and type ```i``` for the ILP formulation.

The details of the optional keywords are described below:

```-lp1```: Lower bound for the parameter p1. Default is 1.

```-up1```: Upper bound for the parameter p1. Default is 5.

```-lp2```: Lower bound for the parameter p2. Default is 1.

```-lp2```: Upper bound for the parameter p2. Default is 100.

You may change the simulated annealing parameters using the following keywords:

```-reads```: Number of reads for simulated annealing. Default is 100.

```-steps```: Number of steps for simulated annealing. Default is 1000.

```-lb```: Lower bound for beta range. Default is 5.

```-ub```: Upper bound for beta range. Default is 100.

```-schedule```: Schedule type. Default is geometric


### Experiments

You can run simulated annealing and quantum annealing experiments using edge-based or ilp formulations. D-Wave hybrid solver can be also used to run the experiments. In order to use quantum annealing or the hybrid solver, you should have the necessary access to D-Wave.

To run experiments using a specific instance, please type the following command.

```
python experiment.py instances dwave_results e qa -ins complete_0_3_0 -p1 1 -p2 50
``` 

The above command uses the instance named ```complete_0_3_0.npz``` located inside the folder ```instances``` and uses quantum annealing and edge-based formulation (type ```i``` for ILP formulation) to perform the experiments. The penalty parameters p1 is set to 1 and p2 is set to 50. Both the summary of experiment results and the output obtained from D-Wave is stored inside the folder named ```dwave_results```.

To run the experiment for all the instances in a folder, instead of specifying instance name and the penalty parameters from the command line, you should provide the location and the name of the dictionary storing the instance names as the they keys and penalty parameters as the values as a npz file. 

```
python experiment.py inst dwave_results i sa -dictf params -dict pdict
``` 
```params``` is the folder storing the dictionary named ```pdict```. This time, ILP formulation is used with simulated annealing.

For the simulated annealing experiments, the parameters ```-reads```, ```-steps```, ```-lb```, ```-ub``` and ```-schedule``` can be specified. For the quantum annealing experiments the following parameters can be used:

```-time```: Annealing time. Default is 50.

```-reads```: Number of reads for simulated and quantum annealing. Default is 1000.

```-chain```: Chain strength. Default is 2.

### Plotting the results

To plot the heat map of probabilities resulting from the parameter search, please type the following:

```
python  parameter_plots.py edge_params e mp
```

```edge_params``` is the name of the folder storing the results of the parameter search, ```e``` indicates that edge formulation is used (type ```i``` for ILP formulation), and ```mp``` indicates the the probability will be set to 0 if the lowest energy sample does not encode an optimal route (type ```p``` otherwise). If you want to plot a specific instance, us the following keyword:

``` -ins ```:  Name of the instance.

To plot the chain lengths, type the following:

```
python chain_lengths_plot.py
```
The above command plots the chain lengths for the experiment results located under the directories ```data/dwave_results_edge``` and ```data/dwave_results_ilp```. If you want to provide a different directory name, then it is possible using the following commands:

```-resi```: Path of the folder containing the D-Wave experiment results (sample files) for the ILP formulation.

```-rese```: Path of the folder containing the D-Wave experiment results (sample files) for edge-based formulation.  

---

To plot the histogram of energies for the experiment results you should use the following command.

```
python histogram_plots.py e
``` 
```e``` indicates that the edge-based formulation is used. All results are plotted which are located under the default folders ``data/dwave_results_edge``` and ```data/dwave_results_ilp``` depending on the selected formulation. You can also specify a specific result or a different folder using the following keywords.

``` -res```: Path of the folder containing the results.

```-ins```: Name of the file for which the histogram will be plotted. (It should be the complete name of the file e.g. ```complete_0_3_0_1_2_qa_1.8_20```.

---

To plot the probabilities of observing samples encoding feasible and optimal routes for a specific number of cities, use the following command:

```
python probability_plots.py e 3 
```
```e``` indicates that the edge formulation is used and ```3``` is the number of cities. The results are taken from the default folders ```data/dwave_results_edge``` and ```data/dwave_results_ilp``` depending on the selected formulation. You can also specify a different folder using the following keyword:

``` -res```: Path of the folder containing the results.

### Resource analysis for real instances
The real instances used in the paper are located in the directory AFG. To obtain the resource analysis presented in the paper, type the following.  

```
python real_instances.py data AFG tw vars
```

```data``` is the name of the folder storing the instances, ```AFG``` is the name of the folder, ```tw``` is the extension of the instance files, ```vars``` is the name of the folder to store the results. The instances are converted to npz format and stored under the folder ```vars/AFG_npz/```. A file named ```analysis``` containing the summary of the resources is created under the directory ```vars/AFG_count/```. If you add the keyword ```-plot```, then a plot is generated under the directory named ```plots```, comparing the number of variables required for different formulations. 


[1] N. Ascheuer. Hamiltonian Path Problems in the On-line Optimization of Flexible Manufacturing Systems. PhD thesis, Technische Universität Berlin, Berlin, Germany, 1995.