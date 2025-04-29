# nonlinear-kink-oscillations
This project includes analytic formula of nonlinear kink oscillations damped by KHI-induced turbulence (written in IDL/Python) and procedures to fit the data by this function.
The fitting ptocedure is demonstarted by mcmc_models_fitting.pro, which uses the SoBAT Toolkit (https://github.com/Sergey-Anfinogentov/SoBAT).

The original function consists of 9 parameters, see nonlinear_model.pro. Since there are dependencies between some parameters, we reduce the number of parameters from 9 to 7, please see nonlinear_model_v2.pro and nonlinear_model_v3.pro.
Details are presented in a manuscript which will be submitted soon.
