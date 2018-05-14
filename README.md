# Neural Models

This repository contains code I wrote during graduate school to model the responses of V1 neurons. 

PD_Estimate_Variance.m: This script generates fake neuronal responses that have different color tuning in the long-wavelength sensitive and medium-wavelength sensitive (LM) plane of cone contrast space. By fitting the fake responses, we can estimate bias and confidence intervals in the estimated preferred directions. 

ComputeModel.m: This function computes the responses of fake neuron based on the specified model parameters. User must specify the type of model to be fit, the parameters of specified model, and the stimuli to which responses should be generated.

FitModel.m: This function computes an error metric between predicted neural responses and real neural responses. This code is used in an iterative numerical minimization procedure designed to determine the parameter values that minimize the error metric specified by the user (e.g., sum of squared error, negative Poisson log-likelihood, Bernoulli error, etc.)
