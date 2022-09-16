The MATLAB file is involving the data extraction from the live script  https://de.mathworks.com/matlabcentral/fileexchange/74797-reading-and-visualizing-covid-19-data

The section titled fitting the data to the model to find the nearest beta and gamma uses the data from t =67 to the end to be used for the parameter estimation. The MATLAB function fminsearch was used to estimate the infection parameters. We later use these parameters to make comparison between data and the model. 

In section improving fit we have made use of a method available in https://www.frontiersin.org/articles/10.3389/fams.2022.645614/full to make use of time-dependent (data-dependent) contagion rate. 

The function cost = obj_func gives us the error functions to be used in fminsearch. We are using fminsearch to find the best fitting parameters that minimises this error function. 

The function SIR_model_1 is the ODE solver which gives us the infected number of people from the model. 

The function SIR_time_model is the ODE solver which takes the data and uses this data to make use of a time-dependent beta.

The function obj_func_new is for finding the error function for the data from the ODE solver SIR_time_model

If you run the script all together it will return the figures for each four states. 

First figure will be the data

second to fifth will be normal SIR model fitting

sixth to ninth will be the improved SIR model fitting
