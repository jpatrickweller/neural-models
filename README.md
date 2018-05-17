# Neural Models

This repository contains some code for generating fake V1 neural responses from an LN model, then estimating the parameters of that model by fitting it to the dummy data.

PD_Estimate_Variance.m: This script generates fake neuronal responses that have different color tuning in the long-wavelength sensitive and medium-wavelength sensitive (LM) plane of cone contrast space. By fitting the fake responses, we can estimate bias and confidence intervals in the estimated preferred directions. 

The user can specify some feature of the simulated neurons: top firing rate (upper asymptote), baseline firing rate (baseline), slope of the contrast response function (exp), and (negative binomial) variance of the responses (kappa).

The user may also specify how many color directions to sample in the LM Plane, and how many fake neurons should be sampled in each of these directions.

After running, the generated data (including stimuli, responses, and model fits) can be browsed by selecting data points. The current dataset is denoted have a red asterick, and the stimuli, responses, and fitted model are displayed above. This plot is 3D and may be rotated.

The generated data, (including stimuli, responses, and model fits) are saved in the current directory. These data will also re-populate the figure when reloaded.


ComputeModel.m: This function computes the responses of fake neuron based on the specified model parameters. User must specify the type of model to be fit, the parameters of specified model, and the stimuli to which responses should be generated.


FitModel.m: This function computes an error metric between predicted neural responses and real neural responses. This code is used in an iterative numerical minimization procedure designed to determine the parameter values that minimize the error metric specified by the user (e.g., sum of squared error, negative Poisson log-likelihood, Bernoulli error, etc.)
