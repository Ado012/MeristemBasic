# MeristemBasic: 3D dynamics of Arabidopsis thaliana WUS/CLV3 network

MeristemBasic is a hybrid ODE/stochastic based signaling simulator that tracks the changes in concentration of different effectors of the WUS/CLV3 network over time under user defined conditions. Several utilties are provided to prepare and process data both before and after simulation. 

## Getting Started

### Dependencies and requirements

MeristemBasic has been tested to run natively on Ubuntu Linux and Windows through Windows Subsystem for Linux. 


### Installing

Linux

unzip repository to directory of your choice

Install Windows Subsystem for Linux
https://docs.microsoft.com/en-us/windows/wsl/install-win10 (You just need WSL 1 although WSL 2 is also good) 
On updated systems WSL can be installed with a single command through Powershell
wsl --install
Restart system to complete installation
Open Powershell or open and run Ubuntu through the Start Menu
Download the repository from Bitbucket and extract it.



### Usage

Linux

Open Terminal and browse to build-OrgSimBasic-Desktop_Qt_5_15_2_GCC_64bit-Debug

sample run command
./MeristemBasic parameterlist032621FcAcWUS1-1.model StandardSteadyStateWUSGradientEmpty.init rk7 > results032621FcAcWUS1-1

A results file will be generated in the directory (in this case results 021321DcTrWt)


Windows Subsystem for Linux (Windows 10)


Open Powershell and navigate to the extracted repository 
Example (change to wherever the directory is on your system) 
cd C:\Users\Ace\Downloads\aldo012-orgsimbasic-88e21b32508f\aldo012-orgsimbasic-88e21b32508f\build-OrgSimBasic-Desktop_Qt_5_15_2_GCC_64bit-Debug
run the Simulator as you would a local Linux install but with wsl at the beginning if you are using Powershell
example

wsl ./MeristemBasic parameterlist_021321DcTrWt.model StandardSteadyStateWUSGradientEmpty.init rk7 > results021321DcTrWt



## Acknowledgement
The base simulator code is forked from Organism Simulator by Henrik Jonsson used in [1]

[1] Yadav RK, Perales M, Gruel J, et al. Plant stem cell maintenance involves direct transcriptional repression of differentiation program. Mol Syst Biol. 2013;9:654. doi:10.1038/msb.2013.8
