# local imports
from coordinate_conversions import greatCircleDistance
from uncertainty_windows import windows


def similarityDistance(newTrigger, midDist, distanceWindow, verbose=False):
  '''
  Checks if trigger is within distance uncertainty bounds.

  Params:   newTrigger:        new trigger to check if it lies within the
                               defined uncertainty window
            midDist:           median distance coordinate between initial
                               trigger pair, stored as a dataframe containing
                               léongitude and latitude
            distanceWindow:    allowed distance uncertainty window

  Returns:  boolean:           similarity check yields either True or False
  '''

  # check if Great Circle distance is greater than allowed window
  if greatCircleDistance(newTrigger, midDist, verbose=verbose) > distanceWindow:
    return False

  return True


def similarityParams(newTrigger, params, paramMidpoints, paramWindows):
  '''
  Performs similarity check on further parameters besides distance.

  Params:   newTrigger:     trigger to perform similarity check on
            params:         list of parameters to check similarity for, given
                            as a list of strings
            paramMidpoints: midpoint in parameter space with which to compare
                            the parameter space of new triggers
            paramWindows:   allowed uncertainty windows around given parameters,
                            contains one value for each parameter

  Returns:  boolean:        similarity check yields either True or False
  '''

  # loop over desired parameters
  for p in params:

    # check if new trigger lies in an acceptable range compared to midpoint
    # of original trigger pair
    if abs(newTrigger[p] - paramMidpoints[p][0]) > 2 * windows[p][0]:
      return False

  # if none of the checks evaluate to False
  return True