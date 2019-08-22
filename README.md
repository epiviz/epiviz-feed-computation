# Epiviz Feed Computation

This repository sets up the backend computational system for automatic statistical analysis in Epiviz Feed. It is implemented in Python. Each instance of the server is configured using a JSON files which lists datasets to be included in the computations. The statistical analysis methods are designed as modules in this python package that implement a specific API that allows them to be called by the computation server. New modules can be included by the analysts for their specific analyses. 

## Getting Started

### Install Python package Dependencies

The `requirements.txt` contains all the dependencies and can be installed through pip by 

`pip install -r requirements.txt`

### Configuration File (JSON)

The Epiviz Feed computational system uses a JSON configuration file to setup statistical tests and datasets to analyze. This JSON file contains 

- Measurements/datasets to use
- list of statistical tests to compute
- pvalue threshold to filter results (defaults to 0.1) and
- Any annotations/links to the data (this information is displayed in the UI)

An example json file is setup in the repository. This uses datasets hosted at UMD and computations available through this package. 

### Setup datasets

Epiviz data provider is a python package that allows users to query genomic data stored in files or database. To setup an instance of the data provider, visit http://epiviz.org/data-provider.

### Implement new statistical tests

The statistical analysis methods are designed as modules in this python package and implement a specific API that allows them to be called by the computation server.The base of all currently implemented statistical methods is the `StatMethod` class and developers can extend this class to implement custom computational modules. The module also needs to implement a `compute` function that performs the statistical analysis in a specific genomic region.

### Run the server

Epiviz Feed uses websockets to communicate between the UI and the computational server. To run the server, use

`python run.py` (make sure you update the configuration file name in this script)

### Epiviz Feed Application

This ([Epiviz Feed Polymer](https://github.com/epiviz/epiviz_feed_polymer) ) repository contains the code and instructions for running the web application.
