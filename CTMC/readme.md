This file contains codes to implement the continuous-time Markov chain model for different layout configurations under full and partial access controls. 
The detailed information about the model can be found in our preprint. 
To view and experiment with any specific configurations, go to the desired folder.

****
In each sub-directory, to run the CTMC model and compute the steady state distribution, in your computer or server, first make sure that you are under the correct directory, then use the command 'python markov.py'.
To run the CTMC-based simulation, use the command 'python simulation.py'.

The input parameters are initiated in the file 'params.py'. For a full list of input to the model, please refer to our paper. We will only include here a list of input that is easy to change and experiment for the general audience who hopes to see the numerical results under a spectrum of settings. Changes to the other items in the input list will requre thorough understanding towards the model, careful calibration of the input data, and very likely but not necessarily some debugging. For a hassle-free start, you can work with the following input parameters. Note that the first alternative for each item is the default as written in the 'params.py' file.

side = 'single' for a single-sided configuration of your choice, and you may change to a double-sided configuration by replacing the value with 'double'.

control = 'partial' for partial access control, and you may change to full access control by replacing the value with 'full'.

meanSERV = 60 (seconds) for mean service time.

simType = 'mc' for input distributions that are consistent with the Markovian assumptions.

SIM_HOUR = 20 for duration of simulation in each iteration being 20 hrs. You may change to any integer of choice. Please make sure that your input is in integral form.

SIM_ITER = 20 for the number of iterations in the simulation being 20. You may change to any integer of choice. Please make sure that your input is in integral form.

dirname = '' for any desired directory to save the output. You may change to any directory of choice.

suffix = '' for any desired suffix to differentiate the output files of multiple experimental runs with different input combinations. You may change to any suffix of choice.

****

The 'inputSetDef.py' specifies the input values for the sets specified in the paper.

The 'transitionDef.py' specifies the two types of transitions in the CTMC model.

The 'eventDef.py' file specifies the transition definitions as events in the simulation.

The 'valueIter.py' file specifies the functions to calculate the discounted long-term rewards of full and partial access controls.

The files 'utils.py' and 'SLinkedList.py' store the helper functions and objects.
