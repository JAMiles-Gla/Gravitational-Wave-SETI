import numpy as np
import pandas as pd

# to calculate angular distances
from sklearn.metrics.pairwise import haversine_distances
from math import degrees, radians

# for solar system barycentre conversion
from astropy import time, coordinates as coord, units as u

def greatCircleDistance(trigger1, trigger2, verbose=False):
  '''
  Calculate the Great Circle distance between the sky locations of two triggers.

  Params:   trigger1, trigger2:  the trigger pair to find the distance between
  Returns:  greatCircleDist:     Great Circle (haversine) distance in degrees
  '''

  # find haversine distance (used to avoid floating point errors)
  haversine = haversine_distances([[radians(trigger1['lat0']), radians(trigger1['long0'])],
                                   [radians(trigger2['lat0']), radians(trigger2['long0'])]])

  # express distance as a single number in degrees
  greatCircleDist = degrees(haversine[0, 1])

  if verbose == True:
    print(greatCircleDist)

  # return distance
  return greatCircleDist


def angularMidpoint(trigger1, trigger2):
  '''
  Calculate the midpoint along the Great Circle distance between two sky
  locations.

  Params:   trigger1, trigger2:  trigger pair to find midpoint of
  Returns:  midDist:             pseudo-trigger at midpoint location, returned
                                 as a dataframe slice
  '''

  # define longitude and latitude of triggers (in radians)
  lat1, lon1 = radians(trigger1['lat0']), radians(trigger1['long0'])
  lat2, lon2 = radians(trigger2['lat0']), radians(trigger2['long0'])

  # switch to cartesian coordinates
  x1, y1, z1 = np.cos(lat1) * np.cos(lon1), np.cos(lat1) * np.sin(lon1), np.sin(lat1)
  x2, y2, z2 = np.cos(lat2) * np.cos(lon2), np.cos(lat2) * np.sin(lon2), np.sin(lat2)

  # find midpoint in cartesian coordinates
  xmid, ymid, zmid = (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2

  # retransform into latitude and longitude
  lon = np.arctan2(ymid, xmid)
  hyp = np.sqrt(xmid**2 + ymid**2)
  lat = np.arctan2(zmid, hyp)

  # change back to degrees again
  latmid, lonmid = degrees(lat), degrees(lon)

  # create dataframe slice
  midDist = {'lat0': [latmid], 'long0': [lonmid]}
  midDist = pd.DataFrame(data=midDist)

  # return dataframe slice
  return midDist.iloc[0]


def solarSystemBarycentre(dataframe):
  '''
  Converts and stores signal time of arrival at the Solar System Barycentre

  Params:   dataframe:   data for which to find Solar System Barycentre times
  Returns:  df:          updated dataframe containing Barycentre times
  '''

  df = dataframe

  # get lists of RA and DEC coords
  ra = df['phi2'].tolist()
  dec = df['theta2'].tolist()

  # get lists of signal arrival times
  times = df['time0'].tolist()

  # express signal locations in appropriate format
  triggerLoc = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')

  # define location at which signal is received, here, the LIGO site is given
  # although the actual location is the centre of the Earth, but the correction
  # is minimal so this should not have too much of an effect
  LIGO = coord.EarthLocation.of_site('LIGO Hanford Observatory')

  # express time of arrival of signals
  earthTimes = time.Time(times, format='gps', scale='utc', location=LIGO)

  # find time (in seconds) that needs to be added to signal arrival times at
  # Earth to find signal arrival times at Solar System Barycentre
  baryAdd = earthTimes.light_travel_time(triggerLoc).to('second')

  # store barycentre times
  df['baryTime'] = (np.array(times) + np.array(baryAdd)).tolist()

  # return updated dataframe
  return df