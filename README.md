Hello! 

This code is designed for the purposes of finding patterns (as specified by the user) in astrophysical signals.
An example of the codes use can be found in Gravitational-Wave-SETI/main.py (the user will need to supply their own data for foreground and/or background data).
The primary function in this package is the likelihood function which (with the use of the other functions) calculates a Bayesian likelihood ratio of your data fitting a certain pattern.
The strength of the methods used in this code is that it searches all the possibilities where the given pattern can appear in the data - providing the data satisfies the constraints input by the user.

Data format:
  The likelihood function take in the input of a pandas DataFrame (from the python package pandas) with specific labels on its columns.
  The label that the code looks for are gpstime: 'time0', right accession (degrees): 'phi0', declination (degrees): 'theta0'
  If your DataFrame doesn't have these columns the code will give a key error, so change your column names with the pandas rename function as necessary.
  If you don't have the sky locations the n you can still run the code by simply putting in 'None' in to the distance window parameter but this is only recommended for data segments and hence signal patterns that are short in time.
  The likelihood function also requires a sequence function that must be coded in a particular format, which is a s follows:
    sequence_function(a,b), where a and b are the start and end of the sequence list output respectively.

How the code finds patterns:
  In order to locate patterns the likelihood function takes two data points that are within a given angular distance on he sky (defined by the distance window parameter) and uses this pair as the starting point to draw the desired sequence around.
  The code is designed so that the sequence designated up to the limit specified (both sequence function and max sequence are parameters of the likelihood function) gets projected forward, backward, and in between the starting data pair.
  If a sequence length of 1 is chosen then the pattern will not place itself between the pair. 
  The sequence repeats itself up until the data ends (or begins when repeating itself backwards).
  This projection that is created is then used as the model via plotting a Gaussian (of specified width) on each sequence points. 
  This model (specifically the gaussians on at the sequence points) is then compared to the data - in time (further sky positions are ignored) by taking the value of each Gaussian at the location of the nearest data point.
  Using this combined statistic of all the Gaussian values (that are above a certain threshold) a likelihood is calculated.
  This likelihood is compared to a uniform distribution in time to give a the likelihood ratio that the code produces.

Methods strengths and weaknesses:
  This code is good for instances where you have a very periodic signal, as missing data points in your sequence will penalise the likelihood ratio.
  Furthermore if there is a low ratio of foreground data to background data the likelihood ratio will also be punished strictly.
  Given that the code looks at all pairs that are within the given parameters (distance window and minimum considered time separation (minDelta)) it is very though in it's search over a piece of data. 
  Hence, it will not miss a sequence if it is present, but it might assign it a low likelihood ratio if the conditions are not adequate (to many background points).
  This can lead to long computation runs, but as long as the number of data points is under about 1000 this can be expectably computationally cost wise.

likelihood function output: 
   The likelihood function produces the maximum likelihood as well as a pandas DataFrame that gives the following information:
   'likelihood_ratio' - bayesian likelihood ratio (not an odds ratio as prior odds are not considered in this code) calculated by L = 1/(number of combinations the sequence points can go into the data points) * (sum of Gaussian statistics created) * (observation length)**(No. Gaussian statistics created)
   'ith_datum' -  the index of the first data point in the trigger pair.
   'jth_datum' - the index of the second data point in the trigger pair.
   'seq_length' - the number of sequence points searched over the data.
   '#_additional_signal_candidates' - the number of sequence points (aside from the first two) that had a data point close enough in time to be counted.
   'true_trigger_flag' - if #_additional_signal_candidates > 0 then this is 1 other wise it is 0.
   'sequence_period' - the time between the ith and jth data points (the pair).

Extra functionality of likelihood function:
  Included in the likelihood function is the ability to plot the sequence functions that you are testing against the data. In oder to do this simply set the plot parameter to 'True'.
  This will have the effect of plotting the sequence defined by each generating data pair, the sequence function, and the maximum sequence defined by the maxSeq parameter for each sequence plotted (data is not plotted).
  As mentioned before the code can run without any positional information and just search for patterns in the time data that is given to it. For this set distanceWindow = None.
  In order to reduce computational time the an optional parameter has been added for the likelihood function - minDelta. 
  This has the effect of not allowing the function to create sequence pairs that are within the number of seconds specified by this parameter - which will have the effect of preventing fine and computationally expensive searches.

Further improvements to be made to the code:
  Creating a bank of sequences to search that would mimic the stretching in a periodic signal during the year (which will compensate for not having sky positions).
  Changing the sequence model in the code to only stretch over a set sequence amount (as this is more likely to find Astrophysical signals). 
  Designing the code so that it can be run on a GPU to reduce the large computation times (pandas can be run on an INVIDA GUP using the cuDF package)
  Create a new likelihood approximation that can deal with a large number of data points and background data points.
  Create a single likelihood statistic that can be used to judge (against a population of noise models) the significance of a piece of data.
  Include other parameters that can be used to constrain the search and decrease computational times (such as false alarm rate of a signal or chirp mass of a gravitational wave candidate).
  
