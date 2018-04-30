# GeneticNeuralNetwork
A simple implementation of a genetic algorithm to optimize a neural network


## Usage
The program accepts a Genetic Algorithm configuration file and one file for the input and for the output respectively.  
The program runs a genetic algorithm to refine a neural network.  
The genetic algorithm parameters have to be specified into the file:  

- GA_config.txt file  

The input and output must be a (nxm) matrix of double numbers named:  

- GA_inputs.txt and GA_outputs.txt.  

If the files are given in a folder the program runs by giving the option:

- -d path/to/folder/

The program can run also by giving the files separately with the options:  

- -c path/to/config/file.txt  

- -i path/to/inputs/file.txt  

- -o path/to/outputs/file.txt  

The program runs a self generated test and saves self generated input, outputs and a config file if the program is run with option:  

- -t  

Any other option will call for this help (also available with the option -h).  

## Example
![Network Results picture](https://github.com/Pella86/GeneticNeuralNetwork/blob/master/read_me_images/NN_analysis.png)

    
