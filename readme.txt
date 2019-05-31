readme.txt

Simulation code accompanying the manuscript:
"Neuronal coupling by endogenous electric fields: 
Cable theory and applications to coincidence detector neurons in the auditory brainstem"
By JH Goldwyn and J Rinzel
Manuscript available on the arXiv

Matlab (R2012b) simulation code written by JH Goldwyn
Simulation code posted to ModelDB on 8/5/2015

This code makes use of SUNDIALS (Suite of Nonlinear and Differential
Algebraic Equation Solvers) and its interface to Matlab (sundialsTB).

These can be downloaded at the website:

http://computation.llnl.gov/casc/sundials/main.html

Documentation and installation instructions for SUNDIALS and
sundialsTB are also available at that address.

Contents:

MSO_dae.m: A function file that defines and solves the system of
                    equations that model the membrane potential of a
                    MSO neuron (Vm), the extracellular potential in a
                    one-dimensional volume conductor surrounding the
                    neuron (Ve), and the membrane potential of a "test
                    neuron" that does not contribute to Ve but can be
                    influenced by it through ephaptic coupling (Vm
                    TEST). See manuscript for details.  This function
                    is runEphapticMSO.m

runEphapticMSO.m: An m-file that reproduces Figures 8, 9, and 10
		    (excluding panel F) from the manuscript.  Also
		    includes an example showing that ephaptic coupling
		    can alter spike threshold

CableModelGui.m: An m-file that can be used to plot amplitude
	            profiles, phase profiles, and animations of
	            evolution of spatial profiles of passive cable
	            model.  Parameter changes and simulation control
	            is implemented in a graphical user interface (gui)
 
CableModelGui.fig:  A fig-file that is called by CableModelGui.m


Instructions:
* Download m files and install sundialsTB 
* For cable model: execute CableModelGui from matlab command line

* For MSO model: Execute runEphapticMSO from matlab command line or
                 edit m-file and run "sections" of code as desired
                 with cut and paste onto the command line.
