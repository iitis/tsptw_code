# Unconstrained Binary Models of the Travelling Salesman Problem Variants for Quantum Optimization

Person responsible for data: Özlem Salehi (ozlemsalehi [at] gmail.com).

The scripts necessary for generating the results provided in the "Unconstrained Binary Models of the Travelling Salesman Problem Variants for Quantum Optimization".

## Reproducing data
The repository contains the data used in the publication. To generate new data, please follow the instruction below.

### Generating random instances  

To generate random instances, run the following command in the main directory:

```
python tsptw_generator.py -dir instances -no 10 -n 3 -l 0 -u 10 -et 20 -lt 40
```
The generated output is stored in the npz format with the variables n (the number of cities), weight (cost matrix) and window (time-windows). The details of the keywords are described below:

```-dir```: Name of the folder to store instances. Default is "instances".

```-no```: Number of instances to be generated. Default is 10.

```-n```: Number of cities including the depot. Default is 3.

```-l```: Lower bound for the entries in the cost matrix. Default is 0.
                        
```-u```: Upper bound for the entries in the cost matrix. Default is 10.

```et```: Upper bound for the earliest start time. Default is 20.

```lt```: Upper bound for the latest start time. Default is 40.

If you want to generate symmetric instances please add the keyword ```-sym```. If you want to store the instances in the txt format as well, then add the keyword ```-txt```.