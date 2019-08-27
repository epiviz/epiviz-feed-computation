===================
Epiviz Feed Compute
===================


This repository sets up the backend computational system for automatic statistical 
analysis in Epiviz Feed. It is implemented in Python. Each instance of the server 
is configured using a JSON files which lists datasets to be included in the computations. 
The statistical analysis methods are designed as modules in this python package that implement 
a specific API that allows them to be called by the computation server. New modules can be 
included by the analysts for their specific analyses. 


Installation
============

The current repos uses Python 2. 7. The `requirements.txt` contains all the dependencies and
 can be installed through pip by 

.. code-block:: python
    pip install -r requirements.txt

Configuration File (JSON format)
================================

The Epiviz Feed computational system uses a JSON configuration file to setup statistical tests 
and datasets for analysis. This JSON file contains 

* Measurements/datasets
* list of statistical tests
* Statistical threshold (e.g., p-value) to filter results (defaults to 0.1) and
* Annotations/Links to the data (this information is displayed in the "information button" from the interface)

An example json file is setup in the repository. This configuration contains dataset hosted at 
University of Maryland and example computations. 

Setup datasets
==============

Epiviz data provider is a python package that allows users to query genomic data stored in 
files or database. To setup an instance of the data provider, visit http://epiviz.org/data-provider.

Implement new statistical tests
===============================

The statistical analysis methods are modules and independent. 
The base of all currently implemented statistical methods lies in the `StatMethod` 
class and developers can extend this class to implement customized computational modules. 
The module also needs to implement a `compute` function that performs the statistical 
analysis in a specific genomic region.

Running the server
==================

Epiviz Feed uses websockets to communicate between the user interface and 
the computational server. To run the server, use

.. code-block:: python
    python run.py 
    
Note:
    make sure you update the configuration file name in this script

Epiviz Feed Application
=======================

This ([Epiviz Feed Polymer](https://github.com/epiviz/epiviz_feed_polymer)) 
repository contains the code and instructions for running the web interface. 
You have to run both to have a fully functional application.

Note
====

This project has been set up using PyScaffold 3.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.
