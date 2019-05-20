import getopt
import argparse
parser = argparse.ArgumentParser(prog = 'lake_area_monitoring.py',
description = 'Calculate time series of lake surface area given lake name. \
(Example: python lake_area_monitoring.py "Great Lake")')

parser.add_argument('lake_name', help = 'Lake name for which the areas will be calculated', type = str)
parser.add_argument('-ln', '--list_lake_names', help = 'print names of all the lakes', action = 'store_true')
args = parser.parse_args()

# load the required functions
from functions_landsat import *

import ee

ee.Initialize()

jrc = ee.Image("JRC/GSW1_0/GlobalSurfaceWater")
lake1 = ee.FeatureCollection("users/sarinarinab/CitSci_Washington_Lakes")
lake2 = ee.FeatureCollection("users/sntopp/Lakes_ply_20171106")
lake3 = ee.FeatureCollection("users/sntopp/NCSC_Horsepen_Lake")

occ = jrc.select('occurrence')

def calcAreaImage(maskImage):
    return(ee.Image.pixelArea().divide(1000000).mask(maskImage).rename(maskImage.bandNames()))

def convert_to_km2(f):
    # calculate area for the lake polygon in km^2
    return(f.set('polygonArea', f.area(10).divide(1000000)))

def convert_date(f):
    # convert timestamp to date string
    return(f.set('Date', ee.Date(f.get('timestamp'))))

def CalculateAreaGen(aoi, MAXDISTANCE, waterClassification):
    # parameterize a function that calculates lake area

    def CalculateArea(img):
        i = CalculateWaterAddFlagsSR(img, waterClassification) # water classificaton method can be changed here
        aoiBuffered = aoi.buffer(200)

        waterMask = i.select(['waterMask']) # waterMask
        obsMask = i.select(['fmask']).gte(2).rename(['obsMask']) # obstructionMask

        # reconstruct the water surface
        recon = ee.Image(fillWater(waterMask, obsMask, aoiBuffered))

        # print(recon.bandNames().getInfo())

        # calculate lakeMask after reconstructing the water surface
        reconstructedLake = (ExtractChannel(recon.select(['Reconstructed']), aoi.buffer(-50), MAXDISTANCE)
        .rename(['Reconstructed']))

        reconstructedLakeAll = (ExtractChannel(recon.select(['ReconstructedAll']), aoi.buffer(-50), MAXDISTANCE)
        .rename(['ReconstructedAll']))

        rawLake = waterMask.And(reconstructedLake).rename(['lakeMask'])

        recon = recon.select(['obsMask']).addBands(reconstructedLake).addBands(rawLake).addBands(reconstructedLakeAll)

        # construct area image, each pixel value = its area in km^2, each band masked with its corresponding lake mask
        areaImage = (calcAreaImage(recon.select(['lakeMask']))
        .addBands(calcAreaImage(recon.select(['obsMask'])))
        .addBands(calcAreaImage(recon.select(['Reconstructed'])))
        .addBands(calcAreaImage(recon.select(['ReconstructedAll']))))

        # add the pixel areas together to calculate total area
        lakeArea = areaImage.reduceRegion(reducer = ee.Reducer.sum(), geometry = aoiBuffered, scale = 30)

        return(i.addBands(recon).setMulti(lakeArea).copyProperties(recon).set('polygonArea', f.get('polygonArea')))

    return(CalculateArea)


# MAIN SCRIPT STARTS

# set parameters
name = args.lake_name
MAXDISTANCE = 1000 # used in cumulative_cost function to isolate lake area from water area
waterClassification = 'Jones2019_2' # one of ['Jones2019_3', 'Jones2019_2', 'Zou2018']

# load and merge lake polygons and sort them by their polygon area
lake1 = lake1.select(['GNIS_Name', 'AreaSqKm'], ['GNIS_NAME', 'AREASQKM'])
lakes = (lake1.merge(lake2).merge(lake3)
.map(convert_to_km2)
.sort('polygonArea'))

# print out lake names if requested
lakeNames = lakes.aggregate_array('GNIS_NAME')
if args.list_lake_names:
    print('')
    print('Here is the list of all supported lake names:')
    print(lakeNames.getInfo())

# load landsat images
dat = (merge_collections_std_bandnames_collection1tier1_sr()
.filterDate('2017-01-01', '2030-01-01')
.filterMetadata('CLOUD_COVER_LAND', 'less_than', 50))

# specify lake name
f = ee.Feature(lakes.filterMetadata('GNIS_NAME', 'equals', name).first())
aoi = f.geometry()

# parameterize lake area function
CalculateArea = CalculateAreaGen(aoi, MAXDISTANCE, waterClassification)

result = ee.ImageCollection(filterContained(dat, aoi)).map(CalculateArea)

# remove obs where exposed lake area is smaller than half of the lake polygon area
halfArea = ee.Number(f.get('polygonArea')).divide(2)
result = ee.FeatureCollection(result.filterMetadata('lakeMask', 'greater_than', halfArea))

# export the result as a link to a csv file
output = result.map(convert_date).select(propertySelectors = ['Date', 'ReconstructedAll'], newProperties = ['Date', 'ReconstructedArea'], retainGeometry = False)

print('')
print('The csv file for ' + name + ' is located at:')
print(output.getDownloadURL(filetype = 'csv', selectors = ['Date', 'ReconstructedArea'], filename = name))
print('')
