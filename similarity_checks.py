import math
from scipy import stats

# local imports
from GravitationalWaveSETI.similarity_checks import similarityParams, similarityDistance


def search(logLikelihoodStart, totalTime, triggers, sequenceTimes, timeWindows, midDist, distanceWindow, params=None, paramMidpoints=None, paramWindows=None, verbose=False):
  '''
  Searches for triggers in each input time sequence.

  Params:   logLikelihoodStart: default likelihood value if triggers are all from
                                background except initial two
            totalTime:          full time range of data segment
            triggers:           all triggers in data
            sequenceTimes:      central times around which trigger should be found
            timeWindows:        allowed uncertainty windows around given times
            midDist:            coordinates of midpoint of Great Circle distance
                                between initial trigger pair
            distanceWindow:     allowed distance uncertainty window (if None then no similarity distance critera is considered
            params:             list of parameters to check similarity on,
                                given as a list of strings
            paramMidpoints:     midpoint in parameter space with which to compare
                                the parameter space of new triggers
            paramWindows:       allowed uncertainty windows around given parameters,
                                contains one value for each parameter

  Returns:  logLikelihood:      cumulative likelihood that a sequence of signals is
                                extraterrestrial
            signalCandidates:   number of signal candidates in the sequence
  '''

  # initialise log likelihood
  logLikelihood = logLikelihoodStart

  # initialise number of signal candidates in a sequence
  signalCandidates = 0

  # loop over time list:
  for time, window in zip(sequenceTimes, timeWindows):

    # find trigger where distance between trigger and true time location is
    # minimal
    triggers['timeDifference'] = abs(triggers['baryTime'] - time)
    closestTriggerIndex = triggers['timeDifference'].idxmin()
    closestTrigger = triggers.iloc[closestTriggerIndex]

    # check distance similarity with original triggers

    if verbose == True:
      print('Similarity check...', end=' ')

    # evaluate Gaussian at time location of closest trigger
    gaussianStatistic = stats.norm.logpdf(closestTrigger['baryTime'], loc=time, scale=window)

    # additionally check similarity in alternative parameter spaces,
    # if prompted
    if params != None:
      if similarityParams(closestTrigger, params, paramMidpoints, paramWindows) == True:

         # check if statistic is greater than background value
        if gaussianStatistic > - math.log(totalTime):

          # alert that all checks have been passed
          if verbose == True:
            print('Signal candidate!')

          # update number of signal candidates
          signalCandidates += 1

          # add log likelihood for closest trigger to default sum and make
          # sure to offset by the background value term again
          logLikelihood += gaussianStatistic + math.log(totalTime) # why are we adding the background?

    else:

      # check if statistic is greater than background value
      if gaussianStatistic > - math.log(totalTime):

        # alert that all checks have been passed
        if verbose == True:
          print('Signal candidate!')
          #print('baryTime of candidate:',time)
        # update number of signal candidates
        signalCandidates += 1

        # add log likelihood for closest trigger to default sum and make
        # sure to offset by the background value term again
        logLikelihood += (gaussianStatistic + math.log(totalTime)) #should this second sign not be a minus?

  # otherwise, return log likelihood value
  return logLikelihood, signalCandidates

'''
Change log

- optimisation: similarity has already been checked in likelihood function so no need to check it here

'''
