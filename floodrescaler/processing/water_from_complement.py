# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""
__version__ = '2024.05.12'
import pprint, os, datetime, tempfile
#from osgeo import gdal
from qgis import processing
from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsFeatureSink,
                       QgsProcessingException, #still works
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterRasterDestination,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterString,
                       QgsRasterLayer,
                       QgsCoordinateTransformContext,
                       QgsMapLayerStore,
                       QgsProcessingOutputLayerDefinition,
                       QgsProject,
                       QgsProcessingFeatureSourceDefinition,
                       QgsFeatureRequest,
                       QgsVectorLayer,
 
                       )

from qgis.analysis import QgsNativeAlgorithms, QgsRasterCalculatorEntry, QgsRasterCalculator

import pandas as pd
import numpy as np
from osgeo import gdal

 
class WaterFromComp(QgsProcessingAlgorithm):
    """build a WSH or WSE grid from teh DEM and the complement
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT_DEM = 'INPUT_DEM'
    INPUT_WATER = 'INPUT_WATER'
 
    INPUT_TYPE = 'INPUT_TYPE'
    
    #outputs 
    OUTPUT_WATER = 'OUTPUT_WATER'
    
    
 
 
    
    

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return WaterFromComp()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'WaterFromComp'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Water Grid from Complement')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr('FloodRescaling')

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'floodrescaler'

    def shortHelpString(self):

        return self.tr(
            """
            <p>Two tools to generate a Water Surface Elevation (WSE) or Water Surface Height (WSH) grid from its complement and a DEM. 
            Useful for pre-processing other tools.</p>
            <ul>
                <li><strong>WSH to WSE</strong>: Add the WSH to the DEM and mask zeros</li>
                <li><strong>WSE to WSH</strong>: Subtract the DEM from the WSE and set masked values to zero </li>
            </ul>
            
            <h3>Tips and Tricks and Notes</h3>
            Grids need to have identical Geospatial attributes (e.g., extents and shape match)
            Inputs (and outputs) adopt dry=null for WSE and dry=0 for WSH grids
            DEM and Input grids are assumed to be hydraulically consistent (i.e., the same DEM was used to generate the Input)
            
            <h3>Issues and Updates</h3>
            See the <a href="https://github.com/cefect/FloodRescaler">project repository</a> for updates, to post an issue/question, or to request new features.

            <h3>Attribution</h3>
            If you use these tools for your work, please cite <a href="https://doi.org/10.1029/2023WR035100">Bryant et. al., (2023)</a>
            """

        )
 

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """
        
        #=======================================================================
        # control
        #=======================================================================
 

        #=======================================================================
        # INPUTS-----------
        #=======================================================================
        #=======================================================================
        # raster layers
        #=======================================================================
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT_DEM, self.tr('DEM'))
        )
        
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT_WATER, self.tr('Input Water Grid'))
        )
                
 
        
        #=======================================================================
        # pars
        #=======================================================================        
        # add methods parameter with some hints
        param = QgsProcessingParameterString(self.INPUT_TYPE, 
                                             self.tr('Input Water Grid Type'), 
                                             defaultValue='WSH', multiLine=False)
        
        param.setMetadata({'widget_wrapper':{ 'value_hints': ['WSH','WSE'] }})
        self.addParameter(param)
        
        
 
        
        #=======================================================================
        # OUTPUTS
        #======================================================================= 
        obj = QgsProcessingParameterRasterDestination       
   
        self.addParameter(
            obj(self.OUTPUT_WATER, self.tr('Output Water Grid'))
        )
        
        


    def processAlgorithm(self, params, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        #=======================================================================
        # defaults
        #=======================================================================
        feedback.pushInfo('starting w/ \n%s'%(pprint.pformat(params, width=30)))
        
        self._init_algo(params, context, feedback)
        
        #=======================================================================
        # retrieve inputs---
        #=======================================================================
        #=======================================================================
        # input rasters
        #=======================================================================
        def get_rlay(attn):
            """load a raster layer and do some checks"""
            inrlay =  self.parameterAsRasterLayer(params, getattr(self, attn), context)
            if not isinstance(inrlay, QgsRasterLayer):
                raise QgsProcessingException(f'bad type on {attn}')            
 
            return inrlay                
 
 
        
        input_dem = get_rlay('INPUT_DEM')
        
        input_water = get_rlay('INPUT_WATER')
        
        
 
        feedback.pushInfo('finished loading raster layers')
        
        #=======================================================================
        # params
        #=======================================================================
        input_type = self.parameterAsString(params, self.INPUT_TYPE, context)
        kwargs = dict()
        
 
        
        feedback.pushInfo(f'running with \'{input_type}\'\n{kwargs}') 
        
        if feedback.isCanceled():
            return {}
        #=======================================================================
        # exceute specified method----
        #=======================================================================.
 
 
        return self.run_water_grid_convert(input_dem, input_water, input_type, 
                               #context=context, feedback=feedback,
                               **params, **kwargs)
        
    def run_water_grid_convert(self, dem_rlay, water_rlay, water_type,
                               **params):
        
        """convert from WSE or WSH using the DEM"""
        
        
        #=======================================================================
        # defaults
        #=======================================================================
        start = now()

        OUTPUT=self._get_out(self.OUTPUT_WATER)
        assert not OUTPUT==''
        #=======================================================================
        # #check consistency
        #=======================================================================
        assert_spatial_equal(dem_rlay, water_rlay, msg='input grid extent mismatch')  
        
        
        #=======================================================================
        # run based on types
        #=======================================================================
        skwargs = dict(OUTPUT=OUTPUT)
        if water_type == 'WSE':
            func = self.wse_to_wsh            
            
        elif water_type=='WSH':
            func = self.wsh_to_wse
        
        else:
            raise KeyError(water_type)
        
        
        result = func(water_rlay, dem_rlay, **skwargs)
        
        return {self.OUTPUT_WATER:result}
        
    def wse_to_wsh(self, wse_rlay, dem_rlay, 
                   OUTPUT='TEMPORARY_OUTPUT'):
        """convert WSE to WSH
        
         ds['WSHf'] = (ds['WSEf'] - ds['DEM']).fillna(0)
        """
        
        tfp =lambda x:self._tfp(prefix=x, _dir=self.temp_dir)
        
        f=self.feedback
        #cf_kwargs = dict(context=self.context, feedback=f) 
        
        f.pushInfo('\ncomputing WSH from WSE\n=========================\n\n')
        
        #=======================================================================
        # check WSE expectations
        #=======================================================================
        assert getNoDataCount(wse_rlay.source())>0, f'WSE is expected to contain some nulls (dry cells)'
        
        #=======================================================================
        # #get raw delta
        #=======================================================================
        #still has nulls from WSE mask
        delta_fp = self._gdal_calc({ 
            'INPUT_A' : wse_rlay, 'BAND_A' : 1, 'INPUT_B' : dem_rlay, 'BAND_B' : 1, 
            'FORMULA':'A-B',            
            'NO_DATA' : -9999,  'OUTPUT' : tfp('delta_'), 'RTYPE' : 5, 
                   })
        
        f.pushInfo(f'got delta\n    {delta_fp}')
        
        #=======================================================================
        # handle mask
        #=======================================================================
        wsh_fp = processing.run('native:fillnodata', 
                                  { 'BAND' : 1, 'FILL_VALUE' : 0, 
                                   'INPUT' : delta_fp, 'OUTPUT' : OUTPUT }, 
                                  **self.proc_kwargs)['OUTPUT']
                                  
        f.pushInfo(f'nulls filled')
        
        #=======================================================================
        # wrap
        #=======================================================================
        stats_d = self._rasterlayerstatistics(wsh_fp)
        
        f.pushInfo(f'WSH grid computed w/\n%s'%{'%.2f'%stats_d[k] for k in ['MAX', 'MEAN', 'MIN', 'SUM']})
        
        if not stats_d['MIN']==0:            
            f.pushWarning('returned a WSH with min!=0.\nYour DEM and WSE may be inconsistent')
            
        if stats_d['MAX']>99:
            f.pushWarning('WSH result has an unexpectedly high max')
        
        return OUTPUT
            
        
        
    def wsh_to_wse(self, wsh_rlay, dem_rlay, 
                   OUTPUT='TEMPORARY_OUTPUT'):
        """convert WSH to WSE
        
        wse_mar = ma.array(wsh_mar.data+dem_mar.data,mask = np.logical_or(wsh_mar.mask, dem_mar.mask), 
                              fill_value=np.nan)
        """
        
        tfp =lambda x:self._tfp(prefix=x, _dir=self.temp_dir)
        
        f=self.feedback
        #cf_kwargs = dict(context=self.context, feedback=f) 
        
        f.pushInfo('\ncomputing WSE from WSH\n=========================\n\n')
        
        #=======================================================================
        # check WSH meets expectations
        #=======================================================================
        stats_d = self._rasterlayerstatistics(wsh_rlay)
        assert stats_d['MIN']==0, f'WSH layer must have a minimum of zero (negative depths not supported)'
        
        #=======================================================================
        # add
        #=======================================================================
        delta_fp = self._gdal_calc({ 
            'INPUT_A' : wsh_rlay, 'BAND_A' : 1, 'INPUT_B' : dem_rlay, 'BAND_B' : 1, 
            'FORMULA':'A+B',            
            'NO_DATA' : -9999,  'OUTPUT' : tfp('delta_'), 'RTYPE' : 5, 
                   })
        f.pushInfo(f'\n\ncomputed delta\n    {delta_fp}')
        #=======================================================================
        # mask
        #=======================================================================
        mask_fp = self._gdal_calc_mask(wsh_rlay)
        
        wse_fp = self._gdal_calc({ 
            'INPUT_A' : delta_fp, 'BAND_A' : 1, 'INPUT_B' : mask_fp, 'BAND_B' : 1, 
            'FORMULA':'A*B',            
            'NO_DATA' : -9999,  'OUTPUT' : OUTPUT, 'RTYPE' : 5, 
                   })
        
        f.pushInfo(f'\n\napplied mask and got \n    {wse_fp}')
        #=======================================================================
        # check
        #=======================================================================
        stats_d = self._rasterlayerstatistics(wse_fp)
        
        stats_d['null_cnt'] = getNoDataCount(wse_fp)
        
        
        f.pushInfo(f'WSE grid computed w/\n%s'%{'%.2f'%stats_d[k] for k in ['MAX', 'MEAN', 'MIN', 'SUM']})
        
        if not stats_d['null_cnt']>0:
            f.pushWarning(f'WSE result has no nulls (dry cells)')
            
        if stats_d['MIN']==0:
            f.pushWarning('WSE result has a min of zero')
            
        return wse_fp
 
        
    def _init_algo(self, params, context, feedback,
                   temp_dir=None,
                   ):
        """common init for tests"""
        self.proc_kwargs = dict(feedback=feedback, context=context, is_child_algorithm=True)
        self.context, self.feedback, self.params = context, feedback, params
        
        if temp_dir is None:
            temp_dir = os.path.join(tempfile.TemporaryDirectory().name, 'dscale')
            
        if not os.path.exists(temp_dir):os.makedirs(temp_dir)
        self.temp_dir=temp_dir
        
        self.logger = log(feedback=feedback)

        
        
    def _get_out(self, attn):
        return self.parameterAsOutputLayer(self.params, getattr(self, attn), self.context)
        
     
    
    def _gdal_calc(self, pars_d):
        """Raster calculator (gdal:rastercalculator)
        - 0: Byte
        - 1: Int16
        - 2: UInt16
        - 3: UInt32
        - 4: Int32
        - 5: Float32
        - 6: Float64
        """
 
        ofp =  processing.run('gdal:rastercalculator', pars_d, **self.proc_kwargs)['OUTPUT']
        
        if not os.path.exists(ofp):
            raise QgsProcessingException('gdal:rastercalculator failed to get a result for \n%s'%pars_d['FORMULA'])
        
        return ofp
    
    
    def _gdal_calc_mask(self, rlay, OUTPUT='TEMPORARY_OUTPUT'):
        """1=data, null=nodata"""
        return self._gdal_calc({'FORMULA':'A/A', 'INPUT_A':rlay, 'BAND_A':1,'NO_DATA':-9999,'OUTPUT':OUTPUT, 'RTYPE':5})
    
    def _rasterlayerstatistics(self, rlay):
        return processing.run('native:rasterlayerstatistics', {'INPUT':rlay}, **self.proc_kwargs)

    def _tfp(self, prefix=None, suffix='.tif', _dir=None):
        if _dir is None: _dir=self.temp_dir
        return tempfile.NamedTemporaryFile(suffix=suffix,prefix=prefix, dir=_dir).name
        
 

    
#===============================================================================
# HELPERS----------
#===============================================================================
"""cant do imports without creating a provider and a plugin"""


        
def now():
    return datetime.datetime.now()



 

def rlay_get_resolution(rlay):
    
    return (rlay.rasterUnitsPerPixelY() + rlay.rasterUnitsPerPixelX())*0.5

 

class log(object):
    """log emulator"""
    def __init__(self, feedback):
        self.feedback=feedback
    def debug(self, msg):
        self.feedback.pushDebugInfo(msg)
    def info(self, msg):
        self.feedback.pushInfo(msg)


def assert_spatial_equal(left, right, msg='',): 
    """ extents check"""
 
    assert isinstance(left, QgsRasterLayer), f'got bad type :{type(left).__name__}'+ '\n'+msg
    assert isinstance(right, QgsRasterLayer), type(right).__name__+ '\n'+msg
    __tracebackhide__ = True
    
    #===========================================================================
    # crs
    #===========================================================================
    if not left.crs()==right.crs():
        raise QgsProcessingException('crs mismatch')
    #===========================================================================
    # extents
    #===========================================================================
    if not left.extent()==right.extent():
        raise QgsProcessingException('%s != %s extent\n    %s != %s\n    '%(
                left.name(),   right.name(), left.extent(), right.extent()) +msg) 
        
        
    #===========================================================================
    # shape
    #===========================================================================
    #left.height()
    
    lshape = rlay_get_shape(left)
    rshape=rlay_get_shape(right)
    if not lshape==rshape:
        raise QgsProcessingException('%s != %s shape\n    %s != %s\n    '%(
                left.name(),   right.name(), lshape,rshape) +msg)
        

def assert_extent_equal(left, right, msg='',): 
    """ extents check"""
 
    assert isinstance(left, QgsRasterLayer), f'got bad type :{type(left).__name__}'+ '\n'+msg
    assert isinstance(right, QgsRasterLayer), type(right).__name__+ '\n'+msg
    __tracebackhide__ = True
    
    #===========================================================================
    # crs
    #===========================================================================
    if not left.crs()==right.crs():
        raise QgsProcessingException('crs mismatch')
    #===========================================================================
    # extents
    #===========================================================================
    if not left.extent()==right.extent():
        raise QgsProcessingException('%s != %s extent\n    %s != %s\n    '%(
                left.name(),   right.name(), left.extent(), right.extent()) +msg) 


def rlay_get_shape(r):
    return (r.height(), r.width())


def getNoDataCount(fp, **kwargs):
    """2022-05-10: this was returning some nulls
    for rasters where I could not find any nulls"""
    #get raw data
    ar = rlay_to_array(fp, **kwargs)
    
    return np.isnan(ar).astype(int).sum()

def rlay_to_array(rlay_fp, dtype=np.dtype('float32')):
    """context managger?"""
    #get raw data
    ds = gdal.Open(rlay_fp)
    band = ds.GetRasterBand(1)
    
    
    ar_raw = np.array(band.ReadAsArray(), dtype=dtype)
    
    #remove nodata values
    ndval = band.GetNoDataValue()
    
 
    
    del ds
    del band
    
    return np.where(ar_raw==ndval, np.nan, ar_raw)
