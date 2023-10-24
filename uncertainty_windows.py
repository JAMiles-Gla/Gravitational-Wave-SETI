import numpy as np


def windows(trigger1Time, trigger2Time, timeLocs, timeErr, delta):
  '''
  Returns list of 1-sigma error window on each trigger location.

  Params:   trigger1Time:   time of first trigger in pair
            trigger2Time:   time of second trigger in pair
            timeLocs:       two dimensional list containing time locations
                            to be checked for
            errTime:        inherent error on each potential trigger time value,
                            stems mostly from signal reconstruction errors
            delta:          sequence unit, (e.g., the number 2 in a sequence
                            would be marked as 2 * delta)

  Returns:  timeWindows:    list of 1-sigma uncertainty windows around given
                            time locations
  '''

  # initialise uncertainty windows list
  timeWindows = []

  # find number of time steps separating triggers
  triggerSteps = round(abs(trigger1Time - trigger2Time) / delta)

  # loop over list of possible time values
  for time in timeLocs:

    # check which trigger time is closer and calculate uncertainties
    # with respect to the closer trigger
    if abs(time - trigger1Time) > abs(time - trigger2Time):

      # find the number of time steps that separates the new time from the
      # nearest trigger
      timeSteps = round(abs(time - trigger2Time) / delta)

      # add inherent uncertainty around time location to cumulative error
      # propagated through each time step
      timeWindow = timeErr * np.sqrt( (2 * timeSteps**2) / triggerSteps**2 + 1 )

      # add uncertainty value to list
      timeWindows.append(timeWindow)

    else:

      # do the same as above but for the other trigger
      timeSteps = round(abs(time - trigger1Time) / delta)
      timeWindow = timeErr * np.sqrt( (2 * timeSteps**2) / triggerSteps**2 + 1 )
      timeWindows.append(timeWindow)

  # return list of uncertainty windows
  return timeWindows