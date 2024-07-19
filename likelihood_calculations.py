import numpy as np
import math
from pandas import DataFrame, Series

# local imports
from GravitationalWaveSETI.coordinate_conversions import solarSystemBarycentre, angularMidpoint
from GravitationalWaveSETI.similarity_checks import similarityDistance
from GravitationalWaveSETI.time_functions import timeLocations
from GravitationalWaveSETI.trigger_search import search
from GravitationalWaveSETI.plotting_functions import plotTimeLocs, histListLengths


def likelihood(dataframe, distanceWindow, timeWindow, sequence, maxSeq, activeFraction = 1, minDelta = None, params=None, paramMidpoints=None, paramWindows=None, plot=False, verbose=False):
  '''
  Function defines trigger pairs and loops through rest of dataframe
  to determine triggers that are in sequence.

  Params:     dataframe:       background or foreground data as pandas dataframe
              distanceWindow:  allowed distance uncertainty window
              timeWindow:      allowed time uncertainty
              sequence:        sequence function (e.g.: primes, Fibonacci, ...)
              maxSeq:          max number of steps before sequence restarts

  Optional:   activeFraction:  The fraction of time you detector (or detector array)
                               was active for between your first datapoint and your last
              minDelta:        minimum time separation, in seconds, that you would like
                               your sequences points to be created
                               (lower minDelta can lead to long run times)
              params:          list of other parameters to check similarity on,
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

  if distanceWindow == None:  ## adhawk way of removing the similarity distance parameter
    distanceWindow = 360 #the whole sky so that all combinations of points are accepted

    df = dataframe

    df = df.rename(columns={'time0': 'baryTime'})  # thus no solar system baryTime will be needed
    df['phi0'] = Series(0, index=df.index)
    df['theta0'] = Series(0, index=df.index)

  else:  # code proceeds as originally designed
    # rescale all times to Solar System Barycentre
    df = solarSystemBarycentre(dataframe)

  df['long0'] = np.where(df['phi0'] > 180, df['phi0'] - 360, df['phi0'])
  df['lat0'] = 90 - df['theta0']

  # sort data by time for clear forward and backward directions for location
  # search
  df = df.sort_values(by='baryTime', ignore_index=True)
  if verbose == True:
    print('Sorted dataframe:\n', df)

  # calculate and store longitude and latitude of each point by converting
  # phi to longitude and theta to latitude

  # minimum and maximum global time
  minTime = df['baryTime'].min() - timeWindow
  maxTime = df['baryTime'].max() + timeWindow

  totalTime = (maxTime - minTime) * activeFraction # to compensate in the likelihood ratio calculation for the detector being inactive at times
  print('Total time:', totalTime / (24 * 3600), 'days')

  # find minimum sequence step length for computational feasibility
  if minDelta == None: # allowing for a change in minDelta outwith the function
    minDelta = totalTime / 250 # for plot etc, should be small enough that in a day we could find a signal the next day
  if verbose == True:
    print('Min delta:', minDelta, 'seconds')

  logLikelihoodInit = 0

  # loop over rows (triggers)
  for i in range(len(df.index)-1):

    # loop over every other row (trigger)
    for j in range(i+1, len(df.index)):

      if verbose:
        print()
        print('i:', i, 'j:', j)

      # define trigger pair
      t1 = df.iloc[i]
      t2 = df.iloc[j]

      # check similarity based on given params
      if verbose:
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

        # initialise log likelihood value with background
        logLikelihoodStart = logLikelihoodInit

        # loop over sequence list
        for seqList, seqWindows in zip(timeLocs, timeWindows):

          # search for triggers in each sequence list
          logLikelihood, signalCandidates = search(logLikelihoodStart, totalTime, df,
                                                   seqList, seqWindows, midDist,distanceWindow,
                                                   params=params, paramMidpoints=paramMidpoints,
                                                   paramWindows=paramWindows, verbose=verbose)

          # subtract combination statistic
          combinStat = math.perm(len(df.index)-2,len(seqList))

          if verbose:
            print( 'The combination statistic is:', combinStat)

          # removes excess likelihood from the many possible combinations
          # - more combintions implies more chance of randomnly finding a siganl
          # dividing (or subtracting in log space) by the number of posible combinations has the
          # effect of marginalising over the unknown time possitions of sequential points,
          # and meaning sequances with higer point density are less likely.
          if combinStat > 0: # stops code breaking from taking in log(0)
              logLikelihood -= math.log(combinStat)
          else:
              combinStat = 1e-20 # small number ~0

          # add log likelihood to list, specifying by 0 if all triggers come
          # from background and by 1 if some are true triggers
          if logLikelihood != logLikelihoodStart - math.log(combinStat):
            logLikelihoodValues.append([logLikelihood, i, j, len(seqList) + 2, signalCandidates, 1, t2Time - t1Time,])
          else:
            logLikelihoodValues.append([logLikelihood, i, j, len(seqList) + 2, signalCandidates, 0, t2Time - t1Time,])

        # additional operations for plotting time locations
        if plot:

          # extend list of time locations
          plotList.extend(timeLocs)

          # plot time locations
          if len(timeLocs) != 0:
            plotTimeLocs(t1Time, t2Time, timeLocs, sequenceLocs, timeWindows)

            # increase number of plots
            numberOfTemplates += len(timeLocs)

  if plot:

    # plot list lengths if prompted
    histListLengths(plotList)

    # print number of templates
    print('Number of templates trialled in this run:', numberOfTemplates)

  # make sure to add failsafe in case list is completely empty
  if not logLikelihoodValues:
    return np.zeros(6), 0

  # find maximum likelihood value
  maxLogLikelihood = np.max(np.array(logLikelihoodValues)[:, 0])

  # create a dataframe to return
  results = np.transpose(np.array(logLikelihoodValues))
  LogLikelihoodtable = DataFrame({'logLikelihood_ratio':  results[0],
                                  'ith_datum':  results[1],
                                  'jth_datum':  results[2],
                                  'seq_length':  results[3],
                                  '#_additional_signal_candidates': results[4],
                                  'true_trigger_flag': results[5],
                                  'sequence_period': results[6]})
  # return data frame and maximum likelihood value
  return LogLikelihoodtable , maxLogLikelihood

'''

Changes from previous version:

- minDelta changing option
- changed code output to pandas dataframe
- ± timeWindow on lower and upper limit on sequence
- verbose shows the number of combinations of 2 points
- len of log likeliehoodInit changed, as ordering of brackets was wrong
- added a functionality where None can be an input for the distance window which means the code dosn't consider sky possitions
- added a functionality where the fraction of the time the detectors were active can be input into modifying the total time
- code optimised by removing redundant functions:  stats.norm.logpdf(t1Time, loc=t1Time, scale=timeWindow) appears twice but producess the same value
- added a column into the results table that gives the difference in time between between the ith and jth points
- minimum and maximum global times extended by timeWindow
- logLikelihoodStart was set to 0 as first two data points should not be considered in likelihood calculation

'''
