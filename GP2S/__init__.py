"""GPTwoSample

A python package for performing global and time-local two sample tests using Gaussian processes.

Directory structure

demo/
 - a demo script showing how to run GPTwoSample form some demo csv data.
 All useful switches are explained.

doc/
 - Package documentation.
  Later

src/ 
 - actual implementation of GPTwoSample
 There are a number of outdated and up2date versions in here.
 Currently used are twosample.py and twosample_interval_smooth.py
 The first performs the "basic" two sample test; the latter is the Two Sample Test published at GCB 2009.

pygp/
 - Library for Gaussian process models using pythong.
 This directory contains the basic regression model; regression using EP (robust) and a whole family of covariance functions.

mlib/
 - General purpose helper library
 Plotting, statisics, priors, optimisation
 This will be refined in a later release.
"""
