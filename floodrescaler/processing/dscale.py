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
import pprint, os, datetime
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
                       QgsRasterLayer ,
                       QgsCoordinateTransformContext,
                       QgsMapLayerStore
 
                       )

from qgis.analysis import QgsNativeAlgorithms, QgsRasterCalculatorEntry, QgsRasterCalculator
 
class Dscale(QgsProcessingAlgorithm):
    """
    Downscaling port to QGIS processing script
        from  L:\09_REPOS\03_TOOLS\FloodDownscaler\fdsc\control.py
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT_DEM = 'INPUT_DEM'
    INPUT_WSE = 'INPUT_WSE'
 
    INPUT_METHOD = 'INPUT_METH'
    
    #outputs
 
    OUTPUT_WSE = 'OUTPUT_WSE'
 
    
    

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return Dscale()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'fres_dscale'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Downscaling')

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
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        return self.tr(
            """Downscaling  
            """)

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """
        
            #=======================================================================
        # control
        #=======================================================================
        self.methods_d = {            
            #'CostGrow':self.run_cg, 
            'Resample':self.run_Resample, 
            #'TerrainFilter':self.run_tf, 
            #'Schumann14':self.run_s14
            }

        #=======================================================================
        # INPUTS-----------
        #=======================================================================
        #=======================================================================
        # raster layers
        #=======================================================================
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT_DEM, self.tr('DEM s1'))
        )
        
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT_WSE, self.tr('WSE s2'))
        )
                
 
        
        #=======================================================================
        # pars
        #======================================================================= 
 
        
        # add methods parameter with some hints
        param = QgsProcessingParameterString(self.INPUT_METHOD, 
                                             self.tr('downscaling method'), 
                                             defaultValue='Resample', multiLine=False)
        
        param.setMetadata({'widget_wrapper':{ 'value_hints': list(self.methods_d.keys()) }})
        self.addParameter(param)
        
        #=======================================================================
        # OUTPUTS
        #======================================================================= 
        obj = QgsProcessingParameterRasterDestination       
   
        self.addParameter(
            obj(self.OUTPUT_WSE, self.tr('WSE s1'))
        )
        
        


    def processAlgorithm(self, params, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        #=======================================================================
        # defaults
        #=======================================================================
        feedback.pushInfo('starting w/ \n%s'%(pprint.pformat(params, width=30)))
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
 
 
        input_wse = get_rlay('INPUT_WSE')
        input_dem = get_rlay('INPUT_DEM')
        
        #check consistency
        assert_extent_equal(input_dem, input_wse, msg='DEM vs. WSE') 
        feedback.pushInfo('finished loading raster layers')
        
        #=======================================================================
        # params
        #=======================================================================
        mName = self.parameterAsString(params, self.INPUT_METHOD, context)
 
        
        feedback.pushInfo(f'running with \'{mName}\'') 
        
        if feedback.isCanceled():
            return {}
        #=======================================================================
        # exceute specified method----
        #=======================================================================.
 
 
        return self.run_dscale(input_dem, input_wse, mName, context=context, feedback=feedback,
                               **params)
    
    def run_dscale(self,   dem1_rlay,wse2_rlay, method,
                   feedback=None,**kwargs):
        """main runner for downscaling using specified method
        
        made separate function for easy testing
        
        Params
        ------------
 
            
        wse2_fp: str
            filepath to WSE raster layer at low-resolution (to be downscaled)
            
        dem1_fp: str
            filepath to DEM raster layer at high-resolution (used to infer downscaled WSE)
            
        method: str
            downsccaling method to apply
            
            
        Note
        -------
        no AOI clipping is performed. raster layers must have the same spatial extents. 
        see p0_clip_rasters to pre-clip the rasters
        """
        
        #=======================================================================
        # defaults
        #=======================================================================
        start = now()
        
        
        #=======================================================================
        # precheck
        #=======================================================================
        assert_extent_equal(dem1_rlay, wse2_rlay)
        #HydTypes('DEM').assert_fp(dem1_fp)
        #HydTypes('WSE').assert_fp(wse2_fp)
        
        #=======================================================================
        # get downscaling ratio
        #=======================================================================
        downscale = get_resolution_ratio(dem1_rlay, wse2_rlay)
        
        feedback.pushInfo(f'started w/ \'{method}\' and downscale={downscale:.2f}')
        
        #=======================================================================
        # executre
        #=======================================================================
        return self.methods_d[method](dem1_rlay,wse2_rlay, downscale, feedback=feedback, **kwargs)
        
        
    def run_Resample(self,
               _,wse2_rlay,downscale,
               feedback=None, context=None,
               OUTPUT_WSE='TEMPORARY_OUTPUT',
               RESAMPLING=1, 
               **kwargs):
        """simple Resamlper"""
        
        target_res = rlay_get_resolution(wse2_rlay)/downscale
        
        feedback.pushInfo(f'gdalwarp on {wse2_rlay}')
        pars_d = { 'DATA_TYPE' : 0, 'EXTRA' : '', 
                  'INPUT' : wse2_rlay, 
                  'OUTPUT' : OUTPUT_WSE,  
                  'MULTITHREADING' : False, 
                  'NODATA' : -9999, 
                  'OPTIONS' : '',                  
                  'RESAMPLING' : RESAMPLING,  #bilinear
                  'SOURCE_CRS' : None, 
                  'TARGET_CRS' : None, 'TARGET_EXTENT' : None, 'TARGET_EXTENT_CRS' : None, 
                  'TARGET_RESOLUTION' : target_res}
        
        #gives a filepath regardless
        ofp =  processing.run('gdal:warpreproject', pars_d, 
                              context=context, feedback=feedback, is_child_algorithm=True)['OUTPUT']
                              
        if not os.path.exists(ofp):
            raise QgsProcessingException('gdal:warpreproject failed to get a result for \n%s'%pars_d['INTPUT'])
        
        return {self.OUTPUT_WSE:ofp} 
        
        
        
        
        
 

    
#===============================================================================
# HELPERS----------
#===============================================================================
"""cant do imports without creating a provider and a plugin"""


        
def now():
    return datetime.datetime.now()

def get_resolution_ratio( 
                             rlay_s1, #fine
                             rlay_s2, #coarse
                             ):
        
        
    s1 = rlay_get_resolution(rlay_s1)
    s2 = rlay_get_resolution(rlay_s2)
    assert (s2/ s1)>=1.0
    return s2 / s1

def rlay_get_resolution(rlay):
    
    return (rlay.rasterUnitsPerPixelY() + rlay.rasterUnitsPerPixelX())*0.5
def assert_extent_equal(left, right, msg='',): 
    """ extents check"""
 
    assert isinstance(left, QgsRasterLayer), type(left).__name__+ '\n'+msg
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


