# SPEBNR: a Simple Python Environment for statistical estimation of Biochemical Network Robustness

SPEBNR is a simple Python tool that permits estimating the distance between two biochemical systems. 

This tool finds inspiration from SPEAR - a Simple Python Environment for statistical estimation of Adaptation and Reliability - by Castiglioni-Loreti-Tini.

Each system consists of a main component: 
  * a *process* describing the behaviour of the biochemical system. Essentially, a solution with *n* species is modelled by the parallel composition of *n* processes, where each process represents one species and its concentration level.

Two systems are compared via a (pseudo)metric, called the *evolution metric*. Thanks to the possibility of extrapolating process behaviour from that of the system typical of our model, this metric allows us to
  * verify how well a biochemical system is fulfilling its tasks by comparing it with other biochemical systems obtained by perturbing its initial state;
  * compare the behaviour of one biochemical system with respect to different initial process and changes in the initial concentrations.

Via the metric we introduce a dependability property of a biochemical system, called *robustness*.

Robustness is the ability of the program to tolerate internal and external perturbation while preserving the original behaviour. Intuitively, *the more robust the system is, the less the gap between the original and the perturbed behaviour is*.

[spebnr.py](./spebnr.py) is used in order to model the evolution sequence of the system and calculate the distance over systems and the distance between system distributions. Then, it is used to calculate the robustness of the biochemical system.

In [createSystem.py](./createSystem.py) uses SPEBNR to create and model the biochemical system.

We consider three different case studies: the *EnvZ/OmpR Osmoregulatory Signaling System*, the *Bacterial Chemotaxis System* and the *Enzyme Activity at Saturation System*.

In [caseStudies/envZOmpR/envZOmpR.py](./caseStudies/envZOmpR/envZOmpR.py) createSystem is used to model a simple biochemical system in which the *EnvZ/OmpR Osmoregulatory Signaling System* regulates two proteins in the *Escherichia Coli* bacterium. We use our algorithm to evaluate the differences between the original EnvZ/OmpR system and two other systems which are obtained by perturbing its initial state. Then, we use our algorithm to evaluate the differences between the original EnvZ/OmpR system *S* and *20* systems which have an input distance in (*eta_1*-0.1,*eta_1*], for five different values for *eta_1*: 0.3, 0.4, 0.5, 0.6, 0.7.

In [caseStudies/bactChem/bactChem.py](./caseStudies/bactChem/bactChem.py) createSystem is used to model a simple biochemical system in which the *Bacterial Chemotaxis System* in *Escherichia Coli* bacterium. *E. Coli* senses the concentration level of attractant *L*, in order to decide whether to keep moving in the current direction or to make a tumble, thus changing direction randomly. We use our algorithm to evaluate the differences between the original Bacterial Chemotaxis system and two other systems which are obtained by perturbing its initial state -- considering two different input: *L* and *CheR*. Then, for each input - *L* and *CheR* - we use our algorithm to evaluate the differences between the original Bacterial Chemotaxis system *S* and *20* systems which have an input distance in (*eta_1*-0.1,*eta_1*], for five different values for *eta_1*: 0.3, 0.4, 0.5, 0.6, 0.7.

In [caseStudies/enzyme/enzyme.py](./caseStudies/enzyme/enzyme.py) createSystem is used to model a simple biochemical system in which the *Enzyme Activity at Saturation System* where an enzyme *R* produces a molecule *X*. We use our algorithm to evaluate the differences between the original Enzyme system and two other systems which are obtained by perturbing its initial state. Then, we use our algorithm to evaluate the differences between the original Enzyme system *S* and *10* systems which have an input distance in (*eta_1*-0.1,*eta_1*], for five different values for *eta_1*: 0.3, 0.4, 0.5, 0.6, 0.7.

## Download 

To download SPEBNR you have just to clone GitHub project:

```
git clone https://github.com/dmanicardi/spebnr.git
```

Run this command in the folder where you want to download the tool.

## How to run experiments

To run experiments Python3 (>= 3.11) is needed. Moreover, the following Python packages must be available:
  * [numpy](https://numpy.org) >= 1.23.4
  * [scipy](https://scipy.org/) >= 1.10.1
  * [matplotlib](https://matplotlib.org) >= 3.6.2
  * [statistics](https://github.com/digitalemagine/py-statistics) >= 1.0.3.5
  
To install all the required packages you can just run:

```
pip install -r requirements.txt
```

If all the needed packages are installed in your system, you have to add [spebnr.py](./spebnr.py) and [createSystem.py](./createSystem.py) in the directory of the chosen case study. Then, you have to execute the chosen case study e.g.:

```
EnvZOmpR.py
```
