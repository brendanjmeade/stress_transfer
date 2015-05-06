from __future__ import division
import math
import code
import datetime
import urllib
import utm
import os.path
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.basemap.pyproj as pyproj
from scipy import io as sio
from okada_wrapper import dc3d0wrapper, dc3dwrapper
from matplotlib import cm

# For buffering around fault shape
from shapely.ops import cascaded_union
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


# Function to check to see if a string can be converted into a float.  Useful for error checking.
def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# Function to either load or request and load ISC event catalog between two different dates
def getIscEventCatalog(startDateTime, endDateTime, catalogType):
    # catalogType must be either 'COMPREHENSIVE' or 'REVIEWED'

    print 'ISC earthquake catalog data from ' + \
          'start day: ' + startDateTime.strftime('%d/%m/%Y %H:%M:%S') + ' to '\
          'end day', endDateTime.strftime('%d/%m/%Y %H:%M:%S')
    outputFileNameCsvDated = 'iscsearch_' + catalogType + '_' + \
                             startDateTime.strftime('%d_%m_%Y_%H_%M_%S') + \
                             endDateTime.strftime('_%d_%m_%Y_%H_%M_%S') + '.csv'

    # Does file exist for these search parameters?  If so use.  If not pull data from ISC server.
    if os.path.isfile(outputFileNameCsvDated):
        print '    File: ' + outputFileNameCsvDated + ' already exists locally. Using local file.'

    else:
        if catalogType == 'EHB':
            print '    File: ' + outputFileNameCsvDated + ' does not exist locally. Requesting data from ISC server.'
            composedUrl = 'http://colossus.iris.washington.edu/cgi-bin/web-db-v4?request=' + \
                          'COLLECTED&req_agcy=EHB'+ \
                          '&out_format=CATCSV&bot_lat=&top_lat=&left_lon=&right_lon=&ctr_lat=&ctr_lon=&radius=&max_dist_units=deg&searchshape=GLOBAL&srn=&grn=' + \
                          '&start_year=' + startDateTime.strftime('%Y') + \
                          '&start_month=' + startDateTime.strftime('%m') + \
                          '&start_day=' + startDateTime.strftime('%d') + \
                          '&start_time=' + startDateTime.strftime('%H') + '%3A' + startDateTime.strftime('%M') + '%3A' +  startDateTime.strftime('%S') + \
                          '&end_year=' + endDateTime.strftime('%Y') + \
                          '&end_month=' + endDateTime.strftime('%m') + \
                          '&end_day=' + endDateTime.strftime('%d') + \
                          '&end_time=' + endDateTime.strftime('%H') + '%3A' + endDateTime.strftime('%M') + '%3A' +  endDateTime.strftime('%S') + \
                          '&min_dep=&max_dep=&min_mag=&max_mag=&null_mag=on&req_mag_type=Any&req_mag_agcy=Any&include_links=off' + \
                          '&table_owner=ehb'
            urllib.urlretrieve(composedUrl, outputFileNameCsvDated)
            print '    File: ' + outputFileNameCsvDated + ' generated from ISC server request'

        else:
            print '    File: ' + outputFileNameCsvDated + ' does not exist locally. Requesting data from ISC server.'
            composedUrl = 'http://colossus.iris.washington.edu/cgi-bin/web-db-v4?request=' + \
                          catalogType + \
                          '&out_format=CATCSV&bot_lat=&top_lat=&left_lon=&right_lon=&ctr_lat=&ctr_lon=&radius=&max_dist_units=deg&searchshape=GLOBAL&srn=&grn=' + \
                          '&start_year=' + startDateTime.strftime('%Y') + \
                          '&start_month=' + startDateTime.strftime('%m') + \
                          '&start_day=' + startDateTime.strftime('%d') + \
                          '&start_time=' + startDateTime.strftime('%H') + '%3A' + startDateTime.strftime('%M') + '%3A' +  startDateTime.strftime('%S') + \
                          '&end_year=' + endDateTime.strftime('%Y') + \
                          '&end_month=' + endDateTime.strftime('%m') + \
                          '&end_day=' + endDateTime.strftime('%d') + \
                          '&end_time=' + endDateTime.strftime('%H') + '%3A' + endDateTime.strftime('%M') + '%3A' +  endDateTime.strftime('%S') + \
                          '&min_dep=&max_dep=&min_mag=&max_mag=&null_mag=on&req_mag_type=Any&req_mag_agcy=Any&include_links=off'
            urllib.urlretrieve(composedUrl, outputFileNameCsvDated)
            print '    File: ' + outputFileNameCsvDated + ' generated from ISC server request'

    # Read in downloaded file  
    lines = [line.strip() for line in open(outputFileNameCsvDated)]
    if catalogType == 'COMPREHENSIVE':
        startLineIdx = 29 # Start processing at this line to avoid header
    elif catalogType == 'REVIEWED':
        startLineIdx = 34 # Start processing at this line to avoid header
    elif catalogType == 'EHB':
        startLineIdx = 26 # Start processing at this line to avoid header
    endLineIdx = np.shape(lines)[0]-5 # Stop processing at this line to avoid footer

    # Dictionary for storing catalog of events
    catalog = dict()
    catalog['yr'] = []
    catalog['mon'] = []
    catalog['day'] = []
    catalog['hr'] = []
    catalog['min'] = []
    catalog['sec'] = []
    catalog['latitude'] = []
    catalog['longitude'] = []
    catalog['depth'] = []
    catalog['magnitude'] = []
    catalog['datetime'] = []

    for i in range(startLineIdx, endLineIdx):
        firstCommaIndex = lines[i].index(',')
        firstCommaOffset = firstCommaIndex - 7 # 7 is the default value that all other indices are relative to

        # Make sure there is both a depth and a magnitude given
        depthGiven = (isNumber(lines[i][61 + firstCommaOffset:65 + firstCommaOffset]) == True)
        magnitudeGiven = (isNumber(lines[i][91 + firstCommaOffset:94 + firstCommaOffset]) == True)
        processCurrentLine = depthGiven & magnitudeGiven

        if processCurrentLine: # If depth and magnitude are present process the event
            catalog['yr'].append(lines[i][18 + firstCommaOffset:22 + firstCommaOffset])
            catalog['mon'].append(lines[i][23 + firstCommaOffset:25 + firstCommaOffset])
            catalog['day'].append(lines[i][26 + firstCommaOffset:28 + firstCommaOffset])
            catalog['hr'].append(lines[i][29 + firstCommaOffset:31 + firstCommaOffset])
            catalog['min'].append(lines[i][32 + firstCommaOffset:34 + firstCommaOffset])
            catalog['sec'].append(lines[i][35 + firstCommaOffset:37 + firstCommaOffset])
            catalog['latitude'].append(float(lines[i][42 + firstCommaOffset:49 + firstCommaOffset]))
            catalog['longitude'].append(float(lines[i][50 + firstCommaOffset:59 + firstCommaOffset]))
            catalog['depth'].append(-float(lines[i][61 + firstCommaOffset:65 + firstCommaOffset])) # negative sign here to make depth ready for Okada conventions
            catalog['magnitude'].append(float(lines[i][91 + firstCommaOffset:94 + firstCommaOffset]))
            catalog['datetime'].append(datetime.datetime.strptime(catalog['yr'][-1] + ' ' + catalog['mon'][-1] + ' ' + catalog['day'][-1] + ' ' + 
                                                                  catalog['hr'][-1] + ' ' + catalog['min'][-1] + ' ' + catalog['sec'][-1], 
                                                                  '%Y %m %d %H %M %S'))
    return(catalog)


def readSrcmodFile(fileName):
    # Constants for changing dimensions
    km2m = 1e3 # Convert kilometers to meters
    cm2m = 1e-2 # Convert centimeters to meters

    # Load the .mat (HDF5-ish) version of the model geometry and slip distribution
    F = sio.loadmat(fileName + '.mat')
    F = F[fileName]
    F = F[0]

    # Declare fields to populate
    EventSrcmod = dict() # Empty dictionary
    EventSrcmod['angle'] = []
    EventSrcmod['dip'] = []
    EventSrcmod['length'] = []
    EventSrcmod['rake'] = []
    EventSrcmod['slip'] = []
    EventSrcmod['slipDip'] = []
    EventSrcmod['slipStrike'] = []
    EventSrcmod['strike'] = []
    EventSrcmod['width'] = []
    EventSrcmod['x1'] = []
    EventSrcmod['x2'] = []
    EventSrcmod['x3'] = []
    EventSrcmod['x4'] = []
    EventSrcmod['y1'] = []
    EventSrcmod['y2'] = []
    EventSrcmod['y3'] = []
    EventSrcmod['y4'] = []
    EventSrcmod['z1'] = []
    EventSrcmod['z2'] = []
    EventSrcmod['z3'] = []
    EventSrcmod['z4'] = []
    EventSrcmod['patchLongitude'] = []
    EventSrcmod['patchLatitude'] = []

    # Calculated in UTM coordinates after everything else has been calculated
    EventSrcmod['x1Utm'] = []
    EventSrcmod['x2Utm'] = []
    EventSrcmod['x3Utm'] = []
    EventSrcmod['x4Utm'] = []
    EventSrcmod['y1Utm'] = []
    EventSrcmod['y2Utm'] = []
    EventSrcmod['y3Utm'] = []
    EventSrcmod['y4Utm'] = []
    EventSrcmod['z1Utm'] = []
    EventSrcmod['z2Utm'] = []
    EventSrcmod['z3Utm'] = []
    EventSrcmod['z4Utm'] = []

    # Extract values that are universal for the entire rupture
    EventSrcmod['epicenterLatitude'] = F['evLAT'][0][0][0]
    EventSrcmod['epicenterLongitude'] = F['evLON'][0][0][0]
    EventSrcmod['depth'] = F['evDPT'][0][0][0]
    EventSrcmod['date'] = F['evDAT'][0][0]
    EventSrcmod['tag'] = F['evTAG'][0][0] # pretty much the same as the filename
    EventSrcmod['description'] = F['event'][0][0]
    EventSrcmod['magnitude'] = F['srcMwMoS'][0][0][0]
    EventSrcmod['moment'] = F['srcMwMoS'][0][0][1]
    EventSrcmod['flinnEngdahlRegionNo'] = F['FEregionNo'][0][0]
    EventSrcmod['flinnEngdahlRegion'] = F['FEregion'][0][0]

    # Convert epicenter to UTM coordintes
    # Note strange *lat, lon* ordering of arguments next line only
    _, _, EventSrcmod['zoneNumber'], EventSrcmod['zoneLetter'] = utm.from_latlon(EventSrcmod['epicenterLatitude'], EventSrcmod['epicenterLongitude'])
    EventSrcmod['projEpicenter'] = pyproj.Proj(proj='utm', zone=str(EventSrcmod['zoneNumber']) + EventSrcmod['zoneLetter'], ellps='WGS84')
    EventSrcmod['epicenterXUtm'], EventSrcmod['epicenterYUtm'] = EventSrcmod['projEpicenter'](EventSrcmod['epicenterLongitude'], EventSrcmod['epicenterLatitude'])

    # Convert date string to datetime object
    EventSrcmod['datetime'] = datetime.datetime.strptime(EventSrcmod['date'], '%m/%d/%Y')

    # Extract the fault geometry and slip into a single stucture of ungrouped patches
    nPanel = int(F['invSEGM'][0][0][0]) # Count number of panels
    for iPanel in range(1, nPanel + 1): # This index starts at 1 because of naming convention of the fields in the dictionary from the .mat file
        # Extract geometric parameters from this panel common to all patches
        strike = 0
        strike = F['seg' + str(iPanel) + 'AStke'][0][0][0]
        angle = -(strike-90)
        if angle < 0:
            angle = angle + 360
        
        # calculate the length and wide if individual patch elements in current panel
        L = F['seg' + str(iPanel) + 'DimWL'][0][0][1] / np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[1]
        W = F['seg' + str(iPanel) + 'DimWL'][0][0][1] / np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[0]
        rTemp = np.array([[math.cos(math.radians(angle)), -math.sin(math.radians(angle))], [math.sin(math.radians(angle)), math.cos(math.radians(angle))]])
        xTempOrig = np.array([[L/2.0], [0.0]])
        xTempRot = np.dot(rTemp, xTempOrig)
        xTopOffset = xTempRot[0];
        yTopOffset = xTempRot[1];
        zTopOffset = 0;
        xTopBottomOffset = F['seg' + str(iPanel) + 'geoX'][0][1][0] - F['seg' + str(iPanel) + 'geoX'][0][0][0]
        yTopBottomOffset = F['seg' + str(iPanel) + 'geoY'][0][1][0] - F['seg' + str(iPanel) + 'geoY'][0][0][0]
        zTopBottomOffset = F['seg' + str(iPanel) + 'geoZ'][0][1][0] - F['seg' + str(iPanel) + 'geoZ'][0][0][0]
        nDownDip = np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[0]
        nAlongStrike = np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[1]

        # Loops over the down-dip and along-strike patches of the current panel
        for iDownDip in range(0, nDownDip):
            for iAlongStrike in range(0, nAlongStrike):
            
                # Extract top center coordinates of current patch
                xTopCenter = F['seg' + str(iPanel) + 'geoX'][0][iDownDip][iAlongStrike]
                yTopCenter = F['seg' + str(iPanel) + 'geoY'][0][iDownDip][iAlongStrike]
                zTopCenter = F['seg' + str(iPanel) + 'geoZ'][0][iDownDip][iAlongStrike]
                EventSrcmod['patchLongitude'].append(F['seg' + str(iPanel) + 'geoLON'][0][iDownDip][iAlongStrike])
                EventSrcmod['patchLatitude'].append(F['seg' + str(iPanel) + 'geoLAT'][0][iDownDip][iAlongStrike])
            
                # Calculate location of top corners and convert from km to m
                EventSrcmod['x1'].append(km2m * (xTopCenter + xTopOffset))
                EventSrcmod['y1'].append(km2m * (yTopCenter + yTopOffset))
                EventSrcmod['z1'].append(km2m * (zTopCenter - zTopOffset)) # not sure this is right?!
                EventSrcmod['x2'].append(km2m * (xTopCenter - xTopOffset))
                EventSrcmod['y2'].append(km2m * (yTopCenter - yTopOffset))
                EventSrcmod['z2'].append(km2m * (zTopCenter - zTopOffset))

                # Calculate location of bottom corners and convert from km to m
                EventSrcmod['x3'].append(km2m * (xTopCenter + xTopBottomOffset + xTopOffset))
                EventSrcmod['y3'].append(km2m * (yTopCenter + yTopBottomOffset + yTopOffset))
                EventSrcmod['z3'].append(km2m * (zTopCenter + zTopBottomOffset - zTopOffset))
                EventSrcmod['x4'].append(km2m * (xTopCenter + xTopBottomOffset - xTopOffset))
                EventSrcmod['y4'].append(km2m * (yTopCenter + yTopBottomOffset - yTopOffset))
                EventSrcmod['z4'].append(km2m * (zTopCenter + zTopBottomOffset - zTopOffset))

                # Create UTM version of the same
                xTopCenterUtm, yTopCenterUtm = EventSrcmod['projEpicenter'](EventSrcmod['patchLongitude'][-1], EventSrcmod['patchLatitude'][-1], inverse=False)
                EventSrcmod['patchXUtm'] = xTopCenterUtm
                EventSrcmod['patchYUtm'] = yTopCenterUtm
                EventSrcmod['x1Utm'].append(xTopCenterUtm + km2m * xTopOffset)
                EventSrcmod['y1Utm'].append(yTopCenterUtm + km2m * yTopOffset)
                EventSrcmod['z1Utm'].append(km2m * (zTopCenter - zTopOffset)) # not sure this is right?!
                EventSrcmod['x2Utm'].append(xTopCenterUtm - km2m * xTopOffset)
                EventSrcmod['y2Utm'].append(yTopCenterUtm - km2m * yTopOffset)
                EventSrcmod['z2Utm'].append(km2m * (zTopCenter - zTopOffset))
                EventSrcmod['x3Utm'].append(xTopCenterUtm + km2m * (xTopBottomOffset + xTopOffset))
                EventSrcmod['y3Utm'].append(yTopCenterUtm + km2m * (yTopBottomOffset + yTopOffset))
                EventSrcmod['z3Utm'].append(km2m * (zTopCenter + zTopBottomOffset - zTopOffset))
                EventSrcmod['x4Utm'].append(xTopCenterUtm + km2m * (xTopBottomOffset - xTopOffset))
                EventSrcmod['y4Utm'].append(yTopCenterUtm + km2m * (yTopBottomOffset - yTopOffset))
                EventSrcmod['z4Utm'].append(km2m * (zTopCenter + zTopBottomOffset - zTopOffset))
             
                # Extract patch dip, strike, width, and length
                EventSrcmod['dip'].append(F['seg' + str(iPanel) + 'DipAn'][0][0][0])
                EventSrcmod['strike'].append(F['seg' + str(iPanel) + 'AStke'][0][0][0])
                EventSrcmod['dipMean'] = np.mean(np.array(EventSrcmod['dip'])) # not width weighted
                EventSrcmod['strikeMean'] = np.mean(np.array(EventSrcmod['strike'])) # not length weighted
                EventSrcmod['dipMean'] = np.mean(np.array(EventSrcmod['dip'])) # not width weighted
                EventSrcmod['strikeMean'] = np.mean(np.array(EventSrcmod['strike'])) # not length weighted
                EventSrcmod['rake'].append(F['seg' + str(iPanel) + 'RAKE'][0][iDownDip][iAlongStrike])
                EventSrcmod['angle'].append(angle)
                EventSrcmod['width'].append(km2m * W)
                EventSrcmod['length'].append(km2m * L)

                # Extract fault slip
                EventSrcmod['slip'].append(cm2m*(F['seg' + str(iPanel) + 'SLIP'][0][iDownDip][iAlongStrike]))
                rTemp = np.array([[math.cos(math.radians(EventSrcmod['rake'][-1])), 
                                   -math.sin(math.radians(EventSrcmod['rake'][-1]))],
                                  [math.sin(math.radians(EventSrcmod['rake'][-1])),
                                   math.cos(math.radians(EventSrcmod['rake'][-1]))]])
                xTempOrig = np.array([[EventSrcmod['slip'][-1]], [0]])
                xTempRot = np.dot(rTemp, xTempOrig)
                EventSrcmod['slipStrike'].append(xTempRot[0])
                EventSrcmod['slipDip'].append(xTempRot[1])
    return(EventSrcmod)


def calcCfs(Tensor, nVecNormal, nVecInPlane, coefficientOfFriction):
    cfs = np.zeros(np.shape(Tensor['sxx']))
    for iObs in range(0, len(Tensor['sxx'])):
        deltaTau = (Tensor['sxx'][iObs] * nVecNormal[0] * nVecInPlane[0] + 
                    Tensor['sxy'][iObs] * nVecNormal[1] * nVecInPlane[0] + 
                    Tensor['sxz'][iObs] * nVecNormal[2] * nVecInPlane[0] + 
                    Tensor['sxy'][iObs] * nVecNormal[0] * nVecInPlane[1] + 
                    Tensor['syy'][iObs] * nVecNormal[1] * nVecInPlane[1] + 
                    Tensor['syz'][iObs] * nVecNormal[2] * nVecInPlane[1] + 
                    Tensor['sxz'][iObs] * nVecNormal[0] * nVecInPlane[2] + 
                    Tensor['syz'][iObs] * nVecNormal[1] * nVecInPlane[2] + 
                    Tensor['szz'][iObs] * nVecNormal[2] * nVecInPlane[2])
        deltaSigma = (Tensor['sxx'][iObs] * nVecNormal[0] * nVecNormal[0] + 
                      Tensor['sxy'][iObs] * nVecNormal[1] * nVecNormal[0] + 
                      Tensor['sxz'][iObs] * nVecNormal[2] * nVecNormal[0] + 
                      Tensor['sxy'][iObs] * nVecNormal[0] * nVecNormal[1] + 
                      Tensor['syy'][iObs] * nVecNormal[1] * nVecNormal[1] + 
                      Tensor['syz'][iObs] * nVecNormal[2] * nVecNormal[1] + 
                      Tensor['sxz'][iObs] * nVecNormal[0] * nVecNormal[2] + 
                      Tensor['syz'][iObs] * nVecNormal[1] * nVecNormal[2] + 
                      Tensor['szz'][iObs] * nVecNormal[2] * nVecNormal[2])
        cfs[iObs] = deltaTau - coefficientOfFriction * deltaSigma
    return(cfs)


def calcCfsNormal(Tensor, nVecNormal, nVecInPlane, coefficientOfFriction):
    cfsNormal = np.zeros(np.shape(Tensor['sxx']))
    for iObs in range(0, len(Tensor['sxx'])):
        deltaSigma = (Tensor['sxx'][iObs] * nVecNormal[0] * nVecNormal[0] + 
                      Tensor['sxy'][iObs] * nVecNormal[1] * nVecNormal[0] + 
                      Tensor['sxz'][iObs] * nVecNormal[2] * nVecNormal[0] + 
                      Tensor['sxy'][iObs] * nVecNormal[0] * nVecNormal[1] + 
                      Tensor['syy'][iObs] * nVecNormal[1] * nVecNormal[1] + 
                      Tensor['syz'][iObs] * nVecNormal[2] * nVecNormal[1] + 
                      Tensor['sxz'][iObs] * nVecNormal[0] * nVecNormal[2] + 
                      Tensor['syz'][iObs] * nVecNormal[1] * nVecNormal[2] + 
                      Tensor['szz'][iObs] * nVecNormal[2] * nVecNormal[2])
        cfs[iObs] = -coefficientOfFriction * deltaSigma
    return(cfsNormal)


def calcCfsTotal(Tensor, nVecNormal, nVecInPlane, coefficientOfFriction):
    cfsTotal = np.zeros(np.shape(Tensor['sxx']))
    nVecCross = np.cross(nVecNormal, nVecInPlane)
    for iObs in range(0, len(Tensor['sxx'])):
        deltaTau1 = (Tensor['sxx'][iObs] * nVecNormal[0] * nVecInPlane[0] + 
                     Tensor['sxy'][iObs] * nVecNormal[1] * nVecInPlane[0] + 
                     Tensor['sxz'][iObs] * nVecNormal[2] * nVecInPlane[0] + 
                     Tensor['sxy'][iObs] * nVecNormal[0] * nVecInPlane[1] + 
                     Tensor['syy'][iObs] * nVecNormal[1] * nVecInPlane[1] + 
                     Tensor['syz'][iObs] * nVecNormal[2] * nVecInPlane[1] + 
                     Tensor['sxz'][iObs] * nVecNormal[0] * nVecInPlane[2] + 
                     Tensor['syz'][iObs] * nVecNormal[1] * nVecInPlane[2] + 
                     Tensor['szz'][iObs] * nVecNormal[2] * nVecInPlane[2])
        deltaTau2 = (Tensor['sxx'][iObs] * nVecNormal[0] * nVecCross[0] + 
                     Tensor['sxy'][iObs] * nVecNormal[1] * nVecCross[0] + 
                     Tensor['sxz'][iObs] * nVecNormal[2] * nVecCross[0] + 
                     Tensor['sxy'][iObs] * nVecNormal[0] * nVecCross[1] + 
                     Tensor['syy'][iObs] * nVecNormal[1] * nVecCross[1] + 
                     Tensor['syz'][iObs] * nVecNormal[2] * nVecCross[1] + 
                     Tensor['sxz'][iObs] * nVecNormal[0] * nVecCross[2] + 
                     Tensor['syz'][iObs] * nVecNormal[1] * nVecCross[2] + 
                     Tensor['szz'][iObs] * nVecNormal[2] * nVecCross[2])
        deltaSigma = (Tensor['sxx'][iObs] * nVecNormal[0] * nVecNormal[0] + 
                      Tensor['sxy'][iObs] * nVecNormal[1] * nVecNormal[0] + 
                      Tensor['sxz'][iObs] * nVecNormal[2] * nVecNormal[0] + 
                      Tensor['sxy'][iObs] * nVecNormal[0] * nVecNormal[1] + 
                      Tensor['syy'][iObs] * nVecNormal[1] * nVecNormal[1] + 
                      Tensor['syz'][iObs] * nVecNormal[2] * nVecNormal[1] + 
                      Tensor['sxz'][iObs] * nVecNormal[0] * nVecNormal[2] + 
                      Tensor['syz'][iObs] * nVecNormal[1] * nVecNormal[2] + 
                      Tensor['szz'][iObs] * nVecNormal[2] * nVecNormal[2])
        cfsTotal[iObs] = np.abs(deltaTau1) + np.abs(deltaTau2) - coefficientOfFriction * deltaSigma
    return(cfsTotal)



def calcOkadaDisplacementStress(xVec, yVec, zVec, EventSrcmod, lambdaLame, muLame, useUtm):
    # Define the material parameter that Okada's Greens functions is sensitive too
    alpha = (lambdaLame+muLame) / (lambdaLame+2*muLame)

    # Calculate elastic displacement, strain, and stress fields associated with all fault patches
    DisplacementVector = dict()
    DisplacementVector['ux'] = np.zeros(xVec.size)
    DisplacementVector['uy'] = np.zeros(xVec.size)
    DisplacementVector['uz'] = np.zeros(xVec.size)
    StrainTensor = dict()
    StrainTensor['sxx'] = np.zeros(xVec.size)
    StrainTensor['sxy'] = np.zeros(xVec.size)
    StrainTensor['sxz'] = np.zeros(xVec.size)
    StrainTensor['syy'] = np.zeros(xVec.size)
    StrainTensor['syz'] = np.zeros(xVec.size)
    StrainTensor['szz'] = np.zeros(xVec.size)
    StressTensor = dict()
    StressTensor['sxx'] = np.zeros(xVec.size)
    StressTensor['sxy'] = np.zeros(xVec.size)
    StressTensor['sxz'] = np.zeros(xVec.size)
    StressTensor['syy'] = np.zeros(xVec.size)
    StressTensor['syz'] = np.zeros(xVec.size)
    StressTensor['szz'] = np.zeros(xVec.size)

    # Check if UTM coordinates should be used
    if useUtm == True:
        utmString = 'Utm'
    else:
        utmString = ''

    # Loop over each patch and calculation elastic displacements and strains
    for iPatch in range(0, len(EventSrcmod['x1'])): # Loop over source patches
        # print 'patch ' + str(iPatch+1) + ' of ' + str(len(EventSrcmod['x1']))
        for iObs in range(0, len(xVec)): # Loop over observation coordinates
            # Translate and (un)rotate observation coordinates
            xTemp = xVec[iObs]-EventSrcmod['x1'+utmString][iPatch]
            yTemp = yVec[iObs]-EventSrcmod['y1'+utmString][iPatch]
            rTemp = np.array([[math.cos(math.radians(-EventSrcmod['angle'][iPatch])), 
                              -math.sin(math.radians(-EventSrcmod['angle'][iPatch]))],
                              [math.sin(math.radians(-EventSrcmod['angle'][iPatch])),
                               math.cos(math.radians(EventSrcmod['angle'][iPatch]))]])
            xTempOrig = np.array([xTemp, yTemp])
            xTempRot = np.dot(rTemp, xTempOrig)
            xTemp = xTempRot[0]
            yTemp = xTempRot[1]

            # Calculate elastic deformation using Okada 1992 (BSSA)
            # Seven arguments to DC3DWrapper are required:
            # alpha = (lambda + mu) / (lambda + 2 * mu)
            # xo = 3-vector representing the observation point (x, y, z in the original)
            # depth = the depth of the fault origin
            # dip = the dip-angle of the rectangular dislocation surface
            # strike_width = the along-strike range of the surface (al1,al2 in the original)
            # dip_width = the along-dip range of the surface (aw1, aw2 in the original)
            # dislocation = 3-vector representing the direction of motion on the surface (DISL1, DISL2, DISL3)
            success, u, uGrad = dc3dwrapper(alpha, [xTemp, yTemp, zVec[iObs]],
                                            EventSrcmod['z3'][iPatch], EventSrcmod['dip'][iPatch],
                                            [0.0, EventSrcmod['length'][iPatch]],
                                            [0.0, EventSrcmod['width'][iPatch]],
                                            [EventSrcmod['slipStrike'][iPatch], EventSrcmod['slipDip'][iPatch], 0.0])
            DisplacementVector['ux'][iObs] += u[0]
            DisplacementVector['uy'][iObs] += u[1]
            DisplacementVector['uz'][iObs] += u[2]
            strainSxx = uGrad[0, 0]
            strainSxy = 0.5*(uGrad[0, 1] + uGrad[1, 0])
            strainSxz = 0.5*(uGrad[0, 2] + uGrad[2, 0])
            strainSyy = uGrad[1, 1]
            strainSyz = 0.5*(uGrad[1, 2] + uGrad[2, 1])
            strainSzz = uGrad[2, 2]
            StrainTensor['sxx'][iObs] += strainSxx
            StrainTensor['sxy'][iObs] += strainSxy
            StrainTensor['sxz'][iObs] += strainSxz
            StrainTensor['syy'][iObs] += strainSyy
            StrainTensor['syz'][iObs] += strainSyz
            StrainTensor['szz'][iObs] += strainSzz
            StressTensor['sxx'][iObs] += lambdaLame*(strainSxx+strainSyy+strainSzz) + 2*muLame*strainSxx
            StressTensor['sxy'][iObs] += 2*muLame*strainSxy
            StressTensor['sxz'][iObs] += 2*muLame*strainSxz
            StressTensor['syy'][iObs] += lambdaLame*(strainSxx+strainSyy+strainSzz) + 2*muLame*strainSyy
            StressTensor['syz'][iObs] += 2*muLame*strainSyz
            StressTensor['szz'][iObs] += lambdaLame*(strainSxx+strainSyy+strainSzz) + 2*muLame*strainSzz
    return(DisplacementVector, StrainTensor, StressTensor)


def plotSrcmodStressAndEarthquakes(EventSrcmod, xVec, yVec, Cfs, Catalog, obsDepth):
    # Clip CFS values for plotting purposes
    cfsHighIdx1 = (Cfs['cfs']>Cfs['cfsUpperLimit']).nonzero()
    cfsHighIdx2 = (Cfs['cfs']>0).nonzero()
    cfsHighIdx = np.intersect1d(np.array(cfsHighIdx1), np.array(cfsHighIdx2))
    cfsLowIdx1 = (Cfs['cfs']<Cfs['cfsLowerLimit']).nonzero()
    cfsLowIdx2 = (Cfs['cfs']<0).nonzero()
    cfsLowIdx = np.intersect1d(np.array(cfsLowIdx1), np.array(cfsLowIdx2))
    Cfs['cfs'][cfsHighIdx] = Cfs['cfsUpperLimit']
    Cfs['cfs'][cfsLowIdx] = Cfs['cfsLowerLimit']

    # Generate figure showing fault geometry and CFS field
    fig = plt.figure(facecolor='white', figsize=(10, 6), dpi=100)
    plt.subplot(1, 2, 1)
    ax = fig.gca()
    cs = plt.tricontourf(xVec.flatten(), yVec.flatten(), 1e-6*Cfs['cfs'].flatten(), 10, cmap=cm.bwr, origin='lower', hold='on', extend='both')
    cs2 = plt.tricontour(xVec.flatten(), yVec.flatten(), 1e-6*Cfs['cfs'].flatten(), 0, linewidths=1.0, colors='k', origin='lower', hold='on')

    # Draw a black line around the CFS field
    xCircle = 100e3*np.cos(np.arange(0, 2*np.pi+0.01, 0.01))
    yCircle = 100e3*np.sin(np.arange(0, 2*np.pi+0.01, 0.01))
    ax.plot(xCircle + EventSrcmod['epicenterXUtm'], yCircle + EventSrcmod['epicenterYUtm'], color=[0.0, 0.0, 0.0], linewidth=1.0)

    for iPatch in range(0, len(EventSrcmod['x1'])): # Plot the edges of each fault patch fault patches
        ax.plot([EventSrcmod['x1Utm'][iPatch], EventSrcmod['x2Utm'][iPatch]], [EventSrcmod['y1Utm'][iPatch], EventSrcmod['y2Utm'][iPatch]], color=[0.0, 0.0, 0.0], linewidth=0.5)
        ax.plot([EventSrcmod['x2Utm'][iPatch], EventSrcmod['x4Utm'][iPatch]], [EventSrcmod['y2Utm'][iPatch], EventSrcmod['y4Utm'][iPatch]], color=[0.0, 0.0, 0.0], linewidth=0.5)
        ax.plot([EventSrcmod['x1Utm'][iPatch], EventSrcmod['x3Utm'][iPatch]], [EventSrcmod['y1Utm'][iPatch], EventSrcmod['y3Utm'][iPatch]], color=[0.0, 0.0, 0.0], linewidth=0.5)
        ax.plot([EventSrcmod['x3Utm'][iPatch], EventSrcmod['x4Utm'][iPatch]], [EventSrcmod['y3Utm'][iPatch], EventSrcmod['y4Utm'][iPatch]], color=[0.0, 0.0, 0.0], linewidth=0.5)

    # Plot scale bar
    ax.plot([-100e3 + EventSrcmod['epicenterXUtm'], -50e3 + EventSrcmod['epicenterXUtm']], 
            [-125e3 + EventSrcmod['epicenterYUtm'], -125e3 + EventSrcmod['epicenterYUtm']], color=[0.0, 0.0, 0.0], linewidth=1)
    ax.text(-75e3 + EventSrcmod['epicenterXUtm'], -130e3 + EventSrcmod['epicenterYUtm'], '50 km', fontsize=10, verticalalignment='center', horizontalalignment='center')

    # Plot text for depth, azimuth, and dip
    ax.text(-75e3 + EventSrcmod['epicenterXUtm'], -110e3 + EventSrcmod['epicenterYUtm'], 'depth = ' + str(obsDepth) + ' (km)', fontsize=10, verticalalignment='center', horizontalalignment='center')
    ax.text(75e3 + EventSrcmod['epicenterXUtm'], -110e3 + EventSrcmod['epicenterYUtm'], 'azimuth = ' + '%0.2f' % Cfs['faultAzimuth'] + ' (deg)', fontsize=10, verticalalignment='center', horizontalalignment='center')
    ax.text(75e3 + EventSrcmod['epicenterXUtm'], -120e3 + EventSrcmod['epicenterYUtm'], 'dip = ' +  '%0.2f' % Cfs['faultDip'] + ' (deg)', fontsize=10, verticalalignment='center', horizontalalignment='center')

    # Plot ISC earthquake locations if they are close enough to the epicenter
    for iIsc in range(0, len(Catalog['xUtm'])):
        ax.plot(Catalog['xUtm'][iIsc], Catalog['yUtm'][iIsc], marker='o', color='white', linestyle='none',
                markerfacecoloralt='gray', markersize=5, alpha=1.0)

    # Standard decorations
    plt.title(EventSrcmod['tag'])
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.axis('equal')
    plt.axis('off') # turning off axes labels...replaced with scale bar
    cbar = plt.colorbar(cs, shrink=0.3, pad=0.05, aspect=25.0,
                        orientation='horizontal', ticks=[-1e-1, 0, 1e-1]) # Make a colorbar for the ContourSet returned by the contourf call
    cbar.ax.set_xlabel(r'$\Delta \mathrm{CFS} \, \mathrm{(MPa)}$')
    cbar.ax.tick_params(length=0)
    ax.set_xlim([-110e3 + EventSrcmod['epicenterXUtm'], 110e3 + EventSrcmod['epicenterXUtm']])
    ax.set_ylim([-140e3 + EventSrcmod['epicenterYUtm'], 120e3 + EventSrcmod['epicenterYUtm']])

    # Plot CFS as a function of distance from SRCMOD epicenter
    plt.subplot(3, 2, 2)
    ax = fig.gca()
    ax.plot([3.0, 5.0], [0.0, 0.0], marker=' ', color='k', linestyle='-', alpha=1.0)
    iscCfsPositives = 0
    iscCfsNegatives = 0
    for iIsc in range(0, len(Catalog['xUtm'])):
        if (Catalog['cfs'][iIsc] > 0):
            iscCfsPositives = iscCfsPositives + 1
            ax.plot(np.log10(Catalog['distanceToEpicenter'][iIsc]), 1e-6*Catalog['cfs'][iIsc], marker='o', color='red', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
        else:
            iscCfsNegatives = iscCfsNegatives + 1
            ax.plot(np.log10(Catalog['distanceToEpicenter'][iIsc]), 1e-6*Catalog['cfs'][iIsc], marker='o', color='blue', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
    ax.set_xlim([3.0, 5.0])
    ax.set_ylim([-0.5e2, 0.5e2])
    plt.xlabel(r'$\log_{10}\,d \, \mathrm{(m)}$')
    plt.ylabel(r'$\Delta \mathrm{CFS} \, \mathrm{(MPa)}$')
    plt.title(r'$N(\Delta\mathrm{CFS}>0) = $' + str(iscCfsPositives) + ', $N(\Delta\mathrm{CFS}<0) = $' + str(iscCfsNegatives))

    # CFS as a function of time
    plt.subplot(3, 2, 4)
    ax = fig.gca()
    ax.plot([0.0, 30.0], [0.0, 0.0], marker=' ', color='k', linestyle='-', alpha=1.0)
    for iIsc in range(0, len(Catalog['xUtm'])):
        if (Catalog['cfs'][iIsc] > 0):
            ax.plot((Catalog['datetime'][iIsc]-EventSrcmod['datetime']).days, 1e-6*Catalog['cfs'][iIsc], marker='o', color='red', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
        else:
            ax.plot((Catalog['datetime'][iIsc]-EventSrcmod['datetime']).days, 1e-6*Catalog['cfs'][iIsc], marker='o', color='blue', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
    ax.set_ylim([-0.5e2, 0.5e2])
    plt.xlabel(r'$t \, \mathrm{(days)}$')
    plt.ylabel(r'$\Delta \mathrm{CFS} \, \mathrm{(MPa)}$')

    # CFS as a function of magnitude
    plt.subplot(3, 2, 6)
    ax = fig.gca()
    ax.plot([3.0, 7.0], [0.0, 0.0], marker=' ', color='k', linestyle='-', alpha=1.0)
    for iIsc in range(0, len(Catalog['xUtm'])):
        if (Catalog['cfs'][iIsc] > 0):
            ax.plot(Catalog['magnitude'][iIsc], 1e-6*Catalog['cfs'][iIsc], marker='o', color='red', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
        else:
            ax.plot(Catalog['magnitude'][iIsc], 1e-6*Catalog['cfs'][iIsc], marker='o', color='blue', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
    ax.set_xlim([3, 7])
    ax.set_ylim([-0.5e2, 0.5e2])
    plt.xlabel(r'$\mathrm{M_W}$')
    plt.ylabel(r'$\Delta \mathrm{CFS} \, \mathrm{(MPa)}$')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)


def diskObservationPoints(obsDepth, xOffset, yOffset):
    # Observation coordinates on a circular grid
    # http://matplotlib.org/1.4.3/examples/pylab_examples/tricontour_demo.html
    n_angles = 30
    n_radii = 30
    radii = np.linspace(0, 1, n_radii)
    radii = 100e3 * np.sqrt(radii)
    angles = np.linspace(0, 2*math.pi, n_angles, endpoint=False)
    angles = np.repeat(angles[..., np.newaxis], n_radii, axis=1)
    angles[:, 1::2] += math.pi/n_angles
    xVec = (radii*np.cos(angles)).flatten()
    yVec = (radii*np.sin(angles)).flatten()
    xVec = xVec + xOffset
    yVec = yVec + yOffset
    zVec = obsDepth*np.ones(xVec.size)
    return(xVec, yVec, zVec)


def cfsVectorsFromAzimuth(faultAzimuth, faultDip):
    # Need to check these updates
    nVecInPlaneRef = [0, 1, 0]
    nVecNormalRef = [1, 0, 0]
    rotationAngle = faultDip - 90 # This is the angle through wich we rotation nVecNormalRef
    rTempAzimuth = np.array([[math.cos(math.radians(faultAzimuth)), math.sin(math.radians(faultAzimuth)), 0],
                             [-math.sin(math.radians(faultAzimuth)), math.cos(math.radians(faultAzimuth)), 0],
                             [0, 0, 1]])
    rTempDip = np.array([[math.cos(math.radians(rotationAngle)), math.sin(math.radians(rotationAngle)), 0],
                         [-math.sin(math.radians(rotationAngle)), math.cos(math.radians(rotationAngle)), 0],
                         [0, 0, 1]])
    nVecInPlane = np.dot(rTempAzimuth, nVecInPlaneRef)
    nVecInPlane = np.dot(rTempDip, nVecInPlane)
    nVecNormal = np.dot(rTempAzimuth, nVecNormalRef)
    nVecNormal = np.dot(rTempDip, nVecNormal)
    return(nVecInPlane, nVecNormal)


def getNearFieldIscEvents(Catalog, nearFieldDistance, EventSrcmod):
    # Convert longitude and latitudes to local UTM coordinates
    Catalog['xUtm'], Catalog['yUtm'] = EventSrcmod['projEpicenter'](Catalog['longitude'], Catalog['latitude'])

    # Determine distanance from ISC events to SRCMOD event and delete those more than NNN km away
    deleteIdx = []
    Catalog['distanceToEpicenter'] = []
    for iIsc in range(0, len(Catalog['xUtm'])):
        srcmodIscDistance = np.sqrt((Catalog['xUtm'][iIsc] - EventSrcmod['epicenterXUtm'])**2 + 
                                    (Catalog['yUtm'][iIsc] - EventSrcmod['epicenterYUtm'])**2)
        Catalog['distanceToEpicenter'].append(srcmodIscDistance)
        if (srcmodIscDistance > nearFieldDistance):
            deleteIdx.append(iIsc)

    # Remove all ISC earthquakes that are not in field of interest from lists in dict Catalog
    deleteIdx = sorted(deleteIdx, reverse=True)
    for iKey in Catalog.keys():
        for iDeleteIdx in range(0, len(deleteIdx)):
            del Catalog[iKey][deleteIdx[iDeleteIdx]]
    return(Catalog)


def calcMaximumShear(tensor):
    # Calculate maximum shear (stress or strain) in three-dimensions.
    # Not yet visualized or checked.
    maximumShear = []
    for iObs in range(0, len(tensor['sxx'])): # Compute invariants at each observation coordinate
        matrix = np.array([[tensor['sxx'][iObs], tensor['sxy'][iObs], tensor['sxz'][iObs]],
                           [tensor['sxy'][iObs], tensor['syy'][iObs], tensor['syz'][iObs]],
                           [tensor['sxz'][iObs], tensor['syz'][iObs], tensor['szz'][iObs]]])
        # Compute principal strains/stresses using eigenvalues
        eigenValues = list(np.linalg.eigvalsh(matrix))
        maximumShear.append((max(eigenValues) - min(eigenValues)) / 2.0) # maximum shear stress
    return(maximumShear)


def calcTensorInvariants(tensor):
    # Calculate invariants of tensor (stress or strain) in three-dimensions.
    # Not yet visualized or checked.
    invariant1 = []
    invariant2 = []
    invariant3 = []
    for iObs in range(0, len(tensor['sxx'])): # Compute invariants at each observation coordinate
        invariant1.append(tensor['sxx'][iObs] + tensor['syy'][iObs] + tensor['szz'][iObs])
        invariant2.append(tensor['sxx'][iObs] * tensor['syy'][iObs] + 
                          tensor['sxx'][iObs] * tensor['szz'][iObs] +
                          tensor['syy'][iObs] * tensor['szz'][iObs] -
                          tensor['sxy'][iObs] * tensor['sxy'][iObs] -
                          tensor['sxz'][iObs] * tensor['sxz'][iObs] -
                          tensor['syz'][iObs] * tensor['syz'][iObs])
        invariant3.append(tensor['sxx'][iObs] * tensor['syy'][iObs] * tensor['szz'][iObs] +
                          2 * tensor['sxy'][iObs] * tensor['sxz'][iObs] * tensor['syz'][iObs] -
                          tensor['sxx'][iObs] * tensor['syz'][iObs] * tensor['syz'][iObs] -
                          tensor['syy'][iObs] * tensor['sxz'][iObs] * tensor['sxz'][iObs] -
                          tensor['szz'][iObs] * tensor['sxy'][iObs] * tensor['sxy'][iObs])
    return(invariant1, invariant2, invariant3)


def calcDeviatoricTensor(tensor):
    sxxDeviatoric = []
    syyDeviatoric = []
    szzDeviatoric = []
    for iObs in range(0, len(tensor['sxx'])): # Compute deviatoric tensor at each observation coordinate
        isotropicComponent = 1.0 / 3.0 * (tensor['sxx'][iObs] + tensor['syy'][iObs] + tensor['szz'][iObs])
        sxxDeviatoric.append(tensor['sxx'][iObs] - isotropicComponent)
        syyDeviatoric.append(tensor['syy'][iObs] - isotropicComponent)
        szzDeviatoric.append(tensor['szz'][iObs] - isotropicComponent)
    return(sxxDeviatoric, syyDeviatoric, szzDeviatoric)


def calcFaultBuffer(EventSrcmod, distance):
    # Create buffer around fault with shapely
    circles = []
    for iPatch in range(0, len(EventSrcmod['x1'])): # Plot the edges of each fault patch fault patches
        circles.append(Point(EventSrcmod['x1Utm'][iPatch], EventSrcmod['y1Utm'][iPatch]).buffer(distance))
        circles.append(Point(EventSrcmod['x2Utm'][iPatch], EventSrcmod['y2Utm'][iPatch]).buffer(distance))
        circles.append(Point(EventSrcmod['x3Utm'][iPatch], EventSrcmod['y3Utm'][iPatch]).buffer(distance))
        circles.append(Point(EventSrcmod['x4Utm'][iPatch], EventSrcmod['y4Utm'][iPatch]).buffer(distance))
    polygonBuffer = cascaded_union(circles)
    temp = np.array(polygonBuffer.exterior).flatten()
    xBuffer = temp[0::2]
    yBuffer = temp[1::2]
    # polygonBuffer = Polygon(polygonBuffer)
    return(xBuffer, yBuffer, polygonBuffer)


def calcBufferGridPoints(xBuffer, yBuffer, polygonBuffer, spacingGrid):
    xBufferFillVec = np.arange(np.min(xBuffer), np.max(xBuffer), spacingGrid)
    yBufferFillVec = np.arange(np.min(yBuffer), np.max(yBuffer), spacingGrid)
    xBufferFillGrid, yBufferFillGrid = np.meshgrid(xBufferFillVec, yBufferFillVec)

    # Select only those grid points inside of buffered region
    xBufferFillGridVec = xBufferFillGrid.flatten()
    yBufferFillGridVec = yBufferFillGrid.flatten()
    deleteIdx = []
    for iBufferFillGrid in range(0, len(xBufferFillGridVec)):
       candidatePoint = Point(xBufferFillGridVec[iBufferFillGrid], yBufferFillGridVec[iBufferFillGrid])
       isIn = polygonBuffer.contains(candidatePoint)
       if isIn == False:
           deleteIdx.append(iBufferFillGrid)

    xBufferFillGridVec = np.delete(xBufferFillGridVec, deleteIdx)
    yBufferFillGridVec = np.delete(yBufferFillGridVec, deleteIdx)
    return(xBufferFillGridVec, yBufferFillGridVec)


def plotSrcmodStressAndEarthquakesBuffer(EventSrcmod, xVec, yVec, Cfs, Catalog, obsDepth, xBuffer, yBuffer):
    # Clip CFS values for plotting purposes
    cfsHighIdx1 = (Cfs['cfs']>Cfs['cfsUpperLimit']).nonzero()
    cfsHighIdx2 = (Cfs['cfs']>0).nonzero()
    cfsHighIdx = np.intersect1d(np.array(cfsHighIdx1), np.array(cfsHighIdx2))
    cfsLowIdx1 = (Cfs['cfs']<Cfs['cfsLowerLimit']).nonzero()
    cfsLowIdx2 = (Cfs['cfs']<0).nonzero()
    cfsLowIdx = np.intersect1d(np.array(cfsLowIdx1), np.array(cfsLowIdx2))
    Cfs['cfs'][cfsHighIdx] = Cfs['cfsUpperLimit']
    Cfs['cfs'][cfsLowIdx] = Cfs['cfsLowerLimit']

    # Generate figure showing fault geometry and CFS field
    fig = plt.figure(facecolor='white', figsize=(10, 6), dpi=100)
    plt.subplot(1, 2, 1)
    ax = fig.gca()
    cs = plt.tricontourf(xVec.flatten(), yVec.flatten(), 1e-6*Cfs['cfs'].flatten(), 10, cmap=cm.bwr, origin='lower', hold='on', extend='both')
    cs2 = plt.tricontour(xVec.flatten(), yVec.flatten(), 1e-6*Cfs['cfs'].flatten(), 0, linewidths=1.0, colors='k', origin='lower', hold='on')

    # Draw the buffer perimeter
    ax.plot(xBuffer, yBuffer, color=[0.0, 0.0, 0.0], linewidth=1.0)

    for iPatch in range(0, len(EventSrcmod['x1'])): # Plot the edges of each fault patch fault patches
        ax.plot([EventSrcmod['x1Utm'][iPatch], EventSrcmod['x2Utm'][iPatch]], [EventSrcmod['y1Utm'][iPatch], EventSrcmod['y2Utm'][iPatch]], color=[0.0, 0.0, 0.0], linewidth=0.5)
        ax.plot([EventSrcmod['x2Utm'][iPatch], EventSrcmod['x4Utm'][iPatch]], [EventSrcmod['y2Utm'][iPatch], EventSrcmod['y4Utm'][iPatch]], color=[0.0, 0.0, 0.0], linewidth=0.5)
        ax.plot([EventSrcmod['x1Utm'][iPatch], EventSrcmod['x3Utm'][iPatch]], [EventSrcmod['y1Utm'][iPatch], EventSrcmod['y3Utm'][iPatch]], color=[0.0, 0.0, 0.0], linewidth=0.5)
        ax.plot([EventSrcmod['x3Utm'][iPatch], EventSrcmod['x4Utm'][iPatch]], [EventSrcmod['y3Utm'][iPatch], EventSrcmod['y4Utm'][iPatch]], color=[0.0, 0.0, 0.0], linewidth=0.5)

    # Plot scale bar
    ax.plot([-100e3 + EventSrcmod['epicenterXUtm'], -50e3 + EventSrcmod['epicenterXUtm']], 
            [np.min(yBuffer) - 20e3, np.min(yBuffer) - 20e3], color=[0.0, 0.0, 0.0], linewidth=1)
    ax.text(-75e3 + EventSrcmod['epicenterXUtm'], np.min(yBuffer) - 30e3, '50 km', fontsize=10, verticalalignment='center', horizontalalignment='center')

    # Plot text for depth, azimuth, and dip
    ax.text(-75e3 + EventSrcmod['epicenterXUtm'], np.min(yBuffer) - 10e3, r'$d = $' + str(obsDepth) + ' (km)', fontsize=10, verticalalignment='center', horizontalalignment='center')
    ax.text(75e3 + EventSrcmod['epicenterXUtm'], np.min(yBuffer) - 10e3, r'$\alpha = $' + '%0.2f' % Cfs['faultAzimuth'] + ' (deg)', fontsize=10, verticalalignment='center', horizontalalignment='center')
    ax.text(75e3 + EventSrcmod['epicenterXUtm'], np.min(yBuffer) - 20e3, r'$\delta = $' +  '%0.2f' % Cfs['faultDip'] + ' (deg)', fontsize=10, verticalalignment='center', horizontalalignment='center')

    # Plot ISC earthquake locations if they are close enough to the epicenter
    for iIsc in range(0, len(Catalog['xUtm'])):
        ax.plot(Catalog['xUtm'][iIsc], Catalog['yUtm'][iIsc], marker='o', color='white', linestyle='none',
                markerfacecoloralt='gray', markersize=5, alpha=1.0)

    # Standard decorations
    plt.title(EventSrcmod['tag'])
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.axis('equal')
    plt.axis('off') # turning off axes labels...replaced with scale bar
    cbar = plt.colorbar(cs, shrink=0.3, pad=0.05, aspect=25.0,
                        orientation='horizontal', ticks=[-1e-1, 0, 1e-1]) # Make a colorbar for the ContourSet returned by the contourf call
    cbar.ax.set_xlabel(r'$\Delta \mathrm{CFS} \, \mathrm{(MPa)}$')
    cbar.ax.tick_params(length=0)
    ax.set_xlim([np.min(xBuffer) - 10e3, np.max(xBuffer) + 10e3])
    ax.set_ylim([np.min(yBuffer) - 10e3, np.max(yBuffer) + 10e3])

    # Plot CFS as a function of distance from SRCMOD epicenter
    plt.subplot(3, 2, 2)
    ax = fig.gca()
    ax.plot([3.0, 5.0], [0.0, 0.0], marker=' ', color='k', linestyle='-', alpha=1.0)
    iscCfsPositives = 0
    iscCfsNegatives = 0
    for iIsc in range(0, len(Catalog['xUtm'])):
        if (Catalog['cfs'][iIsc] > 0):
            iscCfsPositives = iscCfsPositives + 1
            ax.plot(np.log10(Catalog['distanceToEpicenter'][iIsc]), 1e-6*Catalog['cfs'][iIsc], marker='o', color='red', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
        else:
            iscCfsNegatives = iscCfsNegatives + 1
            ax.plot(np.log10(Catalog['distanceToEpicenter'][iIsc]), 1e-6*Catalog['cfs'][iIsc], marker='o', color='blue', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
    ax.set_xlim([3.0, 5.0])
    ax.set_ylim([-0.5e2, 0.5e2])
    plt.xlabel(r'$\log_{10}\,d \, \mathrm{(m)}$')
    plt.ylabel(r'$\Delta \mathrm{CFS} \, \mathrm{(MPa)}$')
    plt.title(r'$N(\Delta\mathrm{CFS}>0) = $' + str(iscCfsPositives) + ', $N(\Delta\mathrm{CFS}<0) = $' + str(iscCfsNegatives))

    # CFS as a function of time
    plt.subplot(3, 2, 4)
    ax = fig.gca()
    ax.plot([0.0, 30.0], [0.0, 0.0], marker=' ', color='k', linestyle='-', alpha=1.0)
    for iIsc in range(0, len(Catalog['xUtm'])):
        if (Catalog['cfs'][iIsc] > 0):
            ax.plot((Catalog['datetime'][iIsc]-EventSrcmod['datetime']).days, 1e-6*Catalog['cfs'][iIsc], marker='o', color='red', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
        else:
            ax.plot((Catalog['datetime'][iIsc]-EventSrcmod['datetime']).days, 1e-6*Catalog['cfs'][iIsc], marker='o', color='blue', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
    ax.set_ylim([-0.5e2, 0.5e2])
    plt.xlabel(r'$t \, \mathrm{(days)}$')
    plt.ylabel(r'$\Delta \mathrm{CFS} \, \mathrm{(MPa)}$')

    # CFS as a function of magnitude
    plt.subplot(3, 2, 6)
    ax = fig.gca()
    ax.plot([3.0, 7.0], [0.0, 0.0], marker=' ', color='k', linestyle='-', alpha=1.0)
    for iIsc in range(0, len(Catalog['xUtm'])):
        if (Catalog['cfs'][iIsc] > 0):
            ax.plot(Catalog['magnitude'][iIsc], 1e-6*Catalog['cfs'][iIsc], marker='o', color='red', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
        else:
            ax.plot(Catalog['magnitude'][iIsc], 1e-6*Catalog['cfs'][iIsc], marker='o', color='blue', linestyle='none',
                    markerfacecoloralt='gray', markersize=5, alpha=1.0)
    ax.set_xlim([3, 7])
    ax.set_ylim([-0.5e2, 0.5e2])
    plt.xlabel(r'$\mathrm{M_W}$')
    plt.ylabel(r'$\Delta \mathrm{CFS} \, \mathrm{(MPa)}$')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)


def main():
    # Name of Srcmod file to read
    fileName = 's1999HECTOR01SALI'

    # Parameters for CFS calculation and visualization
    lambdaLame = 3e10 # First Lame parameter (Pascals)
    muLame = 3e10 # Second Lame parameter (shear modulus, Pascals)
    coefficientOfFriction = 0.85 # Coefficient of friction
    obsDepth = -5e3; # depth of observation coordinates just for disc visualization
    Cfs = dict()
    Cfs['cfsUpperLimit'] = 1e5; # for visualziation purposes
    Cfs['cfsLowerLimit'] = -1e5; # for visualization purposes
    useUtm = True
    catalogType = 'REVIEWED' # Options are:
                             # 'COMPREHENSIVE': most earthquakes. Not human reviewed
                             # 'REVIEWED': slightly fewer earthquakes. Some human quality control
                             # 'EHB': many fewer earthquakes. Human quality control and precise relocations
    captureDays = 30 # Consider ISC earthquakes for this many days after day of main shock
    nearFieldDistance = 100e3 # Keep only those ISC earthquakes withing this disance of the SRCMOD epicenter
    spacingGrid = 10e3

    # Read in Srcmod fault geometry and slip distribution for this representation of the event
    EventSrcmod = readSrcmodFile(fileName)

    # Generate regular grid over region inside of fault buffer
    xBuffer, yBuffer, polygonBuffer = calcFaultBuffer(EventSrcmod, nearFieldDistance)
    xBufferFillGridVec, yBufferFillGridVec = calcBufferGridPoints(xBuffer, yBuffer, polygonBuffer, spacingGrid)

    # Calculate displacement vector and stress tensor at observation coordinates
    DisplacementVectorBuffer, StrainTensorBuffer, StressTensorBuffer = calcOkadaDisplacementStress(xBufferFillGridVec, yBufferFillGridVec, obsDepth + 0*xBufferFillGridVec, EventSrcmod, lambdaLame, muLame, useUtm)

    # Generate observation coordinates on a regular grid around epicenter
    obsX, obsY, obsZ = diskObservationPoints(obsDepth, EventSrcmod['epicenterXUtm'], EventSrcmod['epicenterYUtm'])

    # Calculate displacement vector and stress tensor at observation coordinates
    DisplacementVector, StrainTensor, StressTensor = calcOkadaDisplacementStress(obsX, obsY, obsZ, EventSrcmod, lambdaLame, muLame, useUtm)

    # Resolve Coulomb failure stresses on reciever plane
    Cfs['faultAzimuth'] = EventSrcmod['strikeMean']
    Cfs['faultDip'] = EventSrcmod['dipMean']
    Cfs['nVecInPlane'], Cfs['nVecNormal'] = cfsVectorsFromAzimuth(Cfs['faultAzimuth'], Cfs['faultDip'])
    Cfs['cfs'] = calcCfs(StressTensor, Cfs['nVecNormal'], Cfs['nVecInPlane'], coefficientOfFriction)
    Cfs['cfs'] = calcCfs(StressTensorBuffer, Cfs['nVecNormal'], Cfs['nVecInPlane'], coefficientOfFriction)


    # Read in ISC data for capture days after the date of the SRCMOD event
    startDateTime = EventSrcmod['datetime'] + datetime.timedelta(days=1) # Start day after event
    endDateTime = startDateTime + datetime.timedelta(days=captureDays) - datetime.timedelta(seconds=1)
    Catalog = getIscEventCatalog(startDateTime, endDateTime, catalogType) # Read ISC events
    Catalog = getNearFieldIscEvents(Catalog, nearFieldDistance, EventSrcmod)

    # Calculate Coulomb failure stress at ISC event locations
    Catalog['cfs'] = []
    for iIsc in range(0, len(Catalog['xUtm'])):
        DisplacementVectorIsc, StrainTensorIsc, StressTensorIsc = calcOkadaDisplacementStress(np.array([Catalog['xUtm'][iIsc]]), 
                                                                                              np.array([Catalog['yUtm'][iIsc]]), 
                                                                                              np.array([Catalog['depth'][iIsc]]), 
                                                                                              EventSrcmod, lambdaLame, muLame, useUtm)
        Catalog['cfs'].append(calcCfs(StressTensorIsc, Cfs['nVecNormal'], Cfs['nVecInPlane'], coefficientOfFriction)[0])
    
    # Calculate maximum shear magnitude (not yet orientation!)
    maximumShear = calcMaximumShear(StressTensorIsc)

    # Calculate stress invariants
    invariant1, invariant2, invariant3 = calcTensorInvariants(StressTensorIsc)

    # Calculate deviatoric stress tensor
    sxxDeviatoric, syyDeviatoric, szzDeviatoric = calcDeviatoricTensor(StressTensorIsc)

    # Plot CFS with SRCMOD event and ISC events
    #plotSrcmodStressAndEarthquakes(EventSrcmod, obsX, obsY, Cfs, Catalog, obsDepth)
    plotSrcmodStressAndEarthquakesBuffer(EventSrcmod, xBufferFillGridVec, yBufferFillGridVec, Cfs, Catalog, obsDepth, xBuffer, yBuffer)
    plt.show()


    # fig, ax = plt.subplots(figsize=(8, 8))
    # ax.plot(xBuffer, yBuffer, color=[0.0, 0.0, 0.0], linewidth=0.5)
    # ax.plot(xBufferFillGridVec, yBufferFillGridVec, 'ro', linewidth=0.0)
    # ax.plot(xBuffer, yBuffer, 'go', linewidth=0.0)

    # ax.relim()
    # ax.autoscale()
    # plt.xlabel('x (m)')
    # plt.ylabel('y (m)')
    # plt.axis('equal')
    # plt.show()

    # Provide keyboard control to interact with variables
    code.interact(local=locals())


if __name__ == '__main__':
   main()
