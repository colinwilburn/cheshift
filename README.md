# Introduction

CheShift is a software for prediction of  13Cα and 13Cβ chemical shifts and validation of protein structures.
\
This fork has updated it to be called by a python script on a machine with PyMOL installed. It has been tested on PyMOL 1.8.4.0 and Python 2.7.17.

# Running

To generate predictions of 13Cα and 13Cβ shifts for a protein structure, in Python 2 do:
\
'''
from cheshift_predictions import predict
predict(pdb_filename, output_filename)
'''
