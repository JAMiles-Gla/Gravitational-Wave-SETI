import numpy as np
import math
from scipy import stats

# local imports
from coordinate_conversions import solarSystemBarycentre, angularMidpoint
from similarity_checks import similarityDistance
from time_functions import timeLocations
from trigger_search import search
from statistics import combinations
from plotting_functions import plotTimeLocs, histListLengths


def likelihood(dataframe, distanceWindow, timeWindow, sequence, maxSeq, params=None, paramMidpoints=None, paramWindows=None, plot=False, verbose=False):
  '''
  Function defines trigger pairs and loops through rest of dataframe
  to determine triggers that are in sequence.

  Params:     dataframe:       background or foreground data as pandas dataframe
              distanceWindow:  allowed distance uncertainty window
              timeWindow:      allowed time uncertainty
              sequence:        sequence function (e.g.: primes, Fibonacci, ...)
              maxSeq:          max number of steps before sequence restarts

  Optional:   params:          list of other parameters to check similarity on,
                               given as a list of strings
              paramMidpoints:  midpoint in parameter space with which to compare
                               the parameter space of new triggers
              paramWindows:    allowed uncertainty windows for each parameter,
                               stored as a dataframe with entries corresponding
                               to given params
              plot:            boolean, determines if lists necessary for plots
                               need to be filled, False by default
              verbose:         boolean, prints status updates to simplify debugging,
                               False by default

  Returns:    logLikelihoodValues:   list of log likelihoods for each sequence
              maxLogLikelihood:      maximum log likelihood value
  '''

  # initialise list of log likelihood values which store the cumulative
  # likelihoods that each signal sequence is extraterrestrial
  logLikelihoodValues = []

  if plot == True:

    # keep track of number of templates that are tried
    numberOfTemplates = 0

    # initialise list plotting variable, if requested
    plotList = []

  # rescale all times to Solar System Barycentre
  df = solarSystemBarycentre(dataframe)

  # sort data by time for clear forward and backward directions for location
  # search
  df = df.sort_values(by='baryTime', ignore_index=True)
  if verbose == True:
    print('Sorted dataframe:\n', df)

  # calculate and store longitude and latitude of each point by converting
  # phi to longitude and theta to latitude
  df['long0'] = np.where(df['phi0'] > 180, df['phi0'] - 360, df['phi0'])
  df['lat0'] = 90 - df['theta0']

  # minimum and maximum global time
  minTime = df['baryTime'].min()
  maxTime = df['baryTime'].max() + 1
  totalTime = maxTime - minTime
  print('Total time:', totalTime / (24 * 3600), 'days')

  # find minimum sequence step length for computational feasibility
  minDelta = totalTime / 250
  if verbose == True:
    print('Min delta:', minDelta)

  # initialised log likelihood value at each iteration (to avoid calculating it
  # each time); two true triggers are assigned likelihoods later, hence the -2
  logLikelihoodInit = - len(df.index - 2) * math.log(totalTime)

  # loop over rows (triggers)
  for i in range(len(df.index)-1):

    # loop over every other row (trigger)
    for j in range(i+1, len(df.index)):

      if verbose == True:
        print()
        print('i:', i, 'j:', j)

      # define trigger pair
      t1 = df.iloc[i]
      t2 = df.iloc[j]

      # check similarity based on given params
      if verbose == True:
        print('Initial similarity check...', end=' ')
      if similarityDistance(t1, t2, distanceWindow, verbose=verbose) == True:

        # find midpoint along Great Circle arc
        midDist = angularMidpoint(t1, t2)

        # time of two triggers
        t1Time = t1['baryTime']
        t2Time = t2['baryTime']

        # find time locations to search in
        timeLocs, timeWindows, sequenceLocs = timeLocations(t1Time, t2Time, minTime, maxTime, sequence, maxSeq,
                                                            minDelta, timeWindow, plot=plot, verbose=verbose)

        # find log likelihood for intial trigger pair
        logLikelihoodT1 = stats.norm.logpdf(t1Time, loc=t1Time, scale=timeWindow)
        logLikelihoodT2 = stats.norm.logpdf(t2Time, loc=t2Time, scale=timeWindow)

        # for injection data
        #if i == 18 and j == 28:
          #print(timeLocs)
          #print(midDist)

        # loop over sequence list
        for seqList, seqWindows in zip(timeLocs, timeWindows):

          # initialise log likelihood value with background
          logLikelihoodStart = logLikelihoodInit + logLikelihoodT1 + logLikelihoodT2

          # search for triggers in each sequence list
          logLikelihood, signalCandidates = search(logLikelihoodStart, totalTime, df,
                                                   seqList, seqWindows, midDist, distanceWindow,
                                                   params=params, paramMidpoints=paramMidpoints,
                                                   paramWindows=paramWindows, verbose=verbose)

          # subtract combination statistic
          combinStat = combinations(len(seqList)+2, len(df.index))
          logLikelihood -= math.log(combinStat)

          # add log likelihood to list, specifying by 0 if all triggers come
          # from background and by 1 if some are true triggers
          if logLikelihood != logLikelihoodStart - math.log(combinStat):
            logLikelihoodValues.append([logLikelihood, i, j, len(seqList) + 2, signalCandidates, 1])
          else:
            logLikelihoodValues.append([logLikelihood, i, j, len(seqList) + 2, signalCandidates, 0])

        # additional operations for plotting time locations
        if plot == True:

          # extend list of time locations
          plotList.extend(timeLocs)

          # plot time locations
          if len(timeLocs) != 0:
            plotTimeLocs(t1Time, t2Time, timeLocs, sequenceLocs, timeWindows)

            # increase number of plots
            numberOfTemplates += len(timeLocs)

  if plot == True:

    # plot list lengths if prompted
    histListLengths(plotList)

    # print number of templates
    print('Number of templates trialled in this run:', numberOfTemplates)

  # make sure to add failsafe in case list is completely empty
  if not logLikelihoodValues:
    return np.zeros(6), 0

  # find maximum likelihood value
  maxLogLikelihood = np.max(np.array(logLikelihoodValues)[:, 0])

  # return list of log likelihoods and maximum likelihood value
  return np.array(logLikelihoodValues), maxLogLikelihood