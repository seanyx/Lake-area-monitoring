import ee

def merge_collections_std_bandnames_collection1tier1_sr():
    """merge landsat 5, 7, 8 collection 1 tier 1 SR imageCollections and standardize band names
    """
    ## standardize band names
    bn8 = ['B1', 'B2', 'B3', 'B4', 'B6', 'pixel_qa', 'B5', 'B7']
    bn7 = ['B1', 'B1', 'B2', 'B3', 'B5', 'pixel_qa', 'B4', 'B7']
    bn5 = ['B1', 'B1', 'B2', 'B3', 'B5', 'pixel_qa', 'B4', 'B7']
    bns = ['uBlue', 'Blue', 'Green', 'Red', 'Swir1', 'BQA', 'Nir', 'Swir2']

    # create a merged collection from landsat 5, 7, and 8
    ls5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR").select(bn5, bns)

    ls7 = (ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")
           .filterDate('1999-04-15', '2003-05-30')
           .select(bn7, bns))

    ls8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR").select(bn8, bns)

    merged = ls5.merge(ls7).merge(ls8)

    return(merged)

def id2Img(id):
    return(ee.Image(merge_collections_std_bandnames_collection1tier1_sr()
    .filterMetadata('LANDSAT_ID', 'equals', id)
    .first()))

def Unpack(bitBand, startingBit, bitWidth):
    # unpacking bit bands
    # see: https:#groups.google.com/forum/#!starred/google-earth-engine-developers/iSV4LwzIW7A
    return (ee.Image(bitBand)
            .rightShift(startingBit)
            .bitwiseAnd(ee.Number(2).pow(ee.Number(bitWidth)).subtract(ee.Number(1)).int()))
def UnpackAllSR(bitBand):
    # apply Unpack function for multiple pixel qualities
    bitInfoSR = {
    'Cloud': [5, 1],
    'CloudShadow': [3, 1],
    'SnowIce': [4, 1],
    'Water': [2, 1]
    }
    unpackedImage = ee.Image.cat([Unpack(bitBand, bitInfoSR[key][0], bitInfoSR[key][1]).rename([key]) for key in bitInfoSR])
    return unpackedImage
def AddFmaskSR(image):
    # # add fmask as a separate band to the input image
    temp = UnpackAllSR(image.select(['BQA']))

    fmask = (temp.select(['Water']).rename(['fmask'])
    .where(temp.select(['SnowIce']), ee.Image(3))
    .where(temp.select(['CloudShadow']), ee.Image(2))
    .where(temp.select(['Cloud']), ee.Image(4))
    .mask(temp.select(['Cloud']).gte(0)))

    return image.addBands(fmask)
def CalcHillShadesSR(image):
    # # calculate hill shades
    mergedDEM = ee.Image("users/eeProject/MERIT").select('elevation').clip(image.geometry())

    shiftDistance = 30
    dp1 = ee.Image.cat(ee.Image(shiftDistance), ee.Image(shiftDistance))
    dp2 = ee.Image.cat(ee.Image(-shiftDistance), ee.Image(shiftDistance))
    dp3 = ee.Image.cat(ee.Image(shiftDistance), ee.Image(-shiftDistance))
    dp4 = ee.Image.cat(ee.Image(-shiftDistance), ee.Image(-shiftDistance))

    SOLAR_AZIMUTH_ANGLE = ee.Number(image.get('SOLAR_AZIMUTH_ANGLE'))
    SOLAR_ZENITH_ANGLE = ee.Number(image.get('SOLAR_ZENITH_ANGLE'))

    hillShade1 = (ee.Terrain.hillshade(mergedDEM.displace(dp1),
    SOLAR_AZIMUTH_ANGLE.add(360), SOLAR_ZENITH_ANGLE))

    hillShade2 = (ee.Terrain.hillshade(mergedDEM.displace(dp2),
    SOLAR_AZIMUTH_ANGLE.add(360), SOLAR_ZENITH_ANGLE))

    hillShade3 = (ee.Terrain.hillshade(mergedDEM.displace(dp3),
    SOLAR_AZIMUTH_ANGLE.add(360), SOLAR_ZENITH_ANGLE))

    hillShade4 = (ee.Terrain.hillshade(mergedDEM.displace(dp4),
    SOLAR_AZIMUTH_ANGLE.add(360), SOLAR_ZENITH_ANGLE))

    # # ee.Algorithms.HillShadow
    hillShade = ee.ImageCollection.fromImages([hillShade1, hillShade2, hillShade3, hillShade4]).mean()

    return hillShade.rename(['hillshade'])

# /* functions to classify water (default) */
def ClassifyWater(imgIn, method = 'Jones2019_2'):

    from functions_waterClassification_Jones2019 import ClassifyWaterJones2019
    if method == 'Jones2019_2':
        return(ClassifyWaterJones2019(imgIn, 2))
    if method == 'Jones2019_3':
        return(ClassifyWaterJones2019(imgIn, 3))
    elif method == 'Zou2018':
        from functions_waterClassification_Zou2018 import ClassifyWaterZou2018
        return(ClassifyWaterZou2018(imgIn))

# /* water function */
def CalculateWaterAddFlagsSR(imgIn, waterMethod = 'Jones2019_2'):
    # waterMethod = typeof waterMethod !== 'undefined' ? waterMethod : 'Jones2019'

    fmask = AddFmaskSR(imgIn).select(['fmask'])

    fmaskUnpacked = (fmask.eq(4).rename(['flag_cloud'])
    .addBands(fmask.eq(2).rename(['flag_cldShadow']))
    .addBands(fmask.eq(3).rename(['flag_snowIce']))
    .addBands(fmask.eq(1).rename(['flag_water'])))

    water = ee.Image(ClassifyWater(imgIn, waterMethod)).rename(['waterMask']).where(test = fmask.gte(2), value = 0)
    hillshade = CalcHillShadesSR(imgIn).rename(['flag_hillshade'])

    imgOut = (ee.Image(water.addBands(fmask).addBands(hillshade).addBands(fmaskUnpacked)
    .setMulti({
        'image_id': imgIn.get('LANDSAT_ID'),
        'timestamp': imgIn.get('system:time_start'),
        'scale': imgIn.projection().nominalScale(),
        'crs': imgIn.projection().crs()
    })))

    return(imgOut)

def ExtractChannel(image, centerline, maxDistance):
    # # extract the channel water bodies from the water mask, based on connectivity to the reference centerline.
    connectedToCl = (image.Not().cumulativeCost(
        source = ee.Image().toByte().paint(centerline, 1).And(image), ## only use the centerline that overlaps with the water mask
        maxDistance = maxDistance,
        geodeticDistance = False).eq(0))

    channel = image.updateMask(connectedToCl).unmask(0).updateMask(image.gte(0)).rename(['channelMask'])
    return(channel)

# /* water interpolation function*/
def fillWater(waterMask, obstructionMask, aoi):

    jrc = ee.Image("JRC/GSW1_0/GlobalSurfaceWater")
    occ = jrc.select(['occurrence'])
    woi = occ.clip(aoi).unmask(0)

    kernel = ee.Kernel.circle(radius = 1, units = 'pixels')

    lakeShore = (waterMask
    .subtract(waterMask.focal_min(kernel = kernel, iterations = 1))
    .subtract(obstructionMask.focal_max(kernel = kernel, iterations = 1))
    .gt(0))

    combinedReducer = (ee.Reducer.mean()
    .combine(ee.Reducer.median(), None, True)
    .combine(ee.Reducer.percentile([0, 5, 10, 50, 90, 100]), None, True)
    .combine(ee.Reducer.stdDev(), None, True))

    shoreLineStats = woi.mask(lakeShore).reduceRegion(reducer = combinedReducer, geometry = aoi)

    shoreLineOcc = shoreLineStats.get('occurrence_p5')

    statusFlag = ee.Algorithms.If(condition = shoreLineOcc, trueCase = 'success', falseCase = 'fail')
    shoreLineOcc = ee.Algorithms.If(condition = shoreLineOcc, trueCase = shoreLineOcc, falseCase = 1)


    # reconstruct lake extent
    reconLakeFull = woi.gte(ee.Image.constant(shoreLineOcc)).rename(['ReconstructedAll'])#.Or(waterMask)
    reconObstructedArea = reconLakeFull.And(obstructionMask)
    reconLake = waterMask.Or(reconObstructedArea).rename(['Reconstructed'])


    return(waterMask.addBands(reconLake).addBands(obstructionMask).addBands(reconLakeFull)
    .setMulti({
        'fillStatus': statusFlag,
        'shoreLineOcc': shoreLineOcc,
        'shoreLineOccStd': shoreLineStats.get('occurrence_stdDev')
    }))

def filterContained(collection, geometry):
    col1 = collection.filterBounds(geometry)
    # define a spatial filter based on contain

    spatialFilter = ee.Filter.isContained(leftField = '.geo', rightField = '.geo', maxError = 10)

    # Define a save all join.
    saveAllJoin = ee.Join.saveAll(matchesKey = 'containedBy')

    # Apply the join.
    containJoined = saveAllJoin.apply(geometry, col1, spatialFilter)

    col2 = ee.List(ee.Feature(containJoined
    .first())
    .get('containedBy'))

    return(col2)
