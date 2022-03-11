# Beyond Research 2021-2

## To execute in MyBinder:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Miagarciaru/Beyond_Research_2021_2/main)

## Introduction:
This repository contains a set of tools for analysis of high energy physics and particle physics. This project was done during the semester 2021-2 at Universidad Nacional de Colombia in a program named "Beyond Research". This program allows students to work with international universities and research centers, having the support of professionals in the field of research over the internship. In this case, we mostly worked on tools on jupyter notebooks that implement analysis on particle physics using the [Open Data of ATLAS Experiment](http://opendata.atlas.cern/release/2020/documentation/) of Run 2 (13 TeV Data). 

Here, you can find two main folders in which you will be able to see some analysis written in python or c++ languages. All of these notebooks contains a brief description of the physics behind the processes in SM (Standard Model) or BYS (Beyond Standar Model), as well as a description of some computational tools implemented along the analysis.

### notebooks-cpp folder:

These notebooks make use of the [C++ language](http://www.cplusplus.com/doc/tutorial/) and is interfaced with [ROOT](https://root.cern.ch/), and is available under this [Github link](https://github.com/Miagarciaru/Beyond_Research_2021_2/tree/main/notebooks-cpp). After cloning/downloading the repository, the only things you need to setup are: you need to have the ROOT framework (see [here](https://root.cern.ch/building-root#quick-start) for a quick start on ROOT) and a [gcc compiler](https://gcc.gnu.org/). 

All of these notebooks have different levels of complexity, so people who are interested on learn about techniques and tools for physics analysis on HEP can do it progressively. All of these notebooks use ROOT tools to perform the analysis, and they show you how to read ROOT files, define histograms, functions, variables, fits, use some cuts for the criteria selection and applied these cuts to get the interested events of each process. At the end of these analysis, you will see some final graphs with the information of any physical variable, such as the missing tranverse energy, the invariant mass of Higgs, Z or W bosons, the number of jets and so on. 

Furthermore, some of these analysis compare data from MC (Monte Carlo) samples and you will see how to scale histograms or stack them in order to see an accurate comparison and how effective the selection criteria was in the analysis. 

### notebooks-python folder:

In this folder you will see at least three notebooks written in python language. However, all of them use alternative tools for ROOT framework

