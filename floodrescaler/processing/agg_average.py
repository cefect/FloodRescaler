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
                       )


 
class AggAverage(QgsProcessingAlgorithm):
    """
    Aggregating (upscalig) flood layers 
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT_DEM = 'INPUT_DEM'
    INPUT_WSE = 'INPUT_WSE'
    INPUT_WSH = 'INPUT_WSH'
    INPUT_UPSCALE = 'INPUT_S'
    INPUT_METHOD = 'INPUT_METH'
    
    #outputs
    OUTPUT_DEM = 'OUTPUT_DEM'
    OUTPUT_WSE = 'OUTPUT_WSE'
    OUTPUT_WSH = 'OUTPUT_WSH'
    
    methods_d = {'WSH Averaging':'direct', 'WSE Averaging':'filter'}

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return AggAverage()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'fres_agg_averaging'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Aggregation with Averaging')

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
            """Aggregating (upscalig) flood layers with methods described in Bryant et. al., (2022).
            Two methods are provided, both use simple averaging to compute the DEM grid.
             \'WSH Averaging\': usees the WSH (water depths) grid as the primary 
            layer, which is aggregated through simple averaging, while the WSE grid is computed via addition
            with the DEM grid. \'WSE Averaging\' first aggregates the WSE grid, before filtering against the DEM, then subtracting to compute the WSH grid.            
            """)

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        #=======================================================================
        # INPUTS-----------
        #=======================================================================
        #=======================================================================
        # raster layers
        #=======================================================================
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT_DEM, self.tr('DEM s1 layer'))
        )
        
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT_WSE, self.tr('WSE s1 layer'))
        )
                
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT_WSH, self.tr('WSH s1 layer'))
        )
        
        #=======================================================================
        # pars
        #======================================================================= 
        self.addParameter(
            QgsProcessingParameterNumber(self.INPUT_UPSCALE, self.tr('upscaling value'), defaultValue=10.0)
        )
        
        # add methods parameter with some hints
        param = QgsProcessingParameterString(self.INPUT_METHOD, self.tr('aggregation method'), defaultValue='WSH Averaging', multiLine=False)
        param.setMetadata({'widget_wrapper':{ 'value_hints': self.methods_d.keys() }})
        self.addParameter(param)
        
        #=======================================================================
        # OUTPUTS
        #======================================================================= 
        obj = QgsProcessingParameterRasterDestination       
        self.addParameter(
            obj(self.OUTPUT_DEM, self.tr('DEM s2 layer'))
        )
        
        self.addParameter(
            obj(self.OUTPUT_WSE, self.tr('WSE s2 layer'))
        )
                
        self.addParameter(
            obj(self.OUTPUT_WSH, self.tr('WSH s2 layer'))
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        #=======================================================================
        # retrieve inputs---
        #=======================================================================
        #=======================================================================
        # rasters
        #=======================================================================
        def get_rlay(attn):
            inrlay =  self.parameterAsRasterLayer(parameters, getattr(self, attn), context)
            if not isinstance(inrlay, QgsRasterLayer):
                raise QgsProcessingException(f'bad type on {attn}')
            
            assert_rlay_simple(inrlay)
            return inrlay
                
 
        input_wsh = get_rlay('INPUT_WSH')
        input_wse = get_rlay('INPUT_WSE')
        input_dem = get_rlay('INPUT_DEM')
        
        #check consistency
        assert_extent_equal(input_dem, input_wse, msg='DEM vs. WSE')
        assert_extent_equal(input_dem, input_wsh, msg='DEM vs. WSH')

        feedback.pushInfo('finished loading raster layers')
        
        #=======================================================================
        # params
        #=======================================================================
        mName = self.parameterAsString(parameters, self.INPUT_METHOD, context)
        scale = self.parameterAsDouble(parameters, self.INPUT_UPSCALE, context)
        
        feedback.pushInfo(f'running with \'{mName}\' and scale={scale:.2f}')
        
        #=======================================================================
        # exceute specified method
        #=======================================================================.
        args = (input_dem, input_wsh,  input_wse, scale)
        if self.methods_d[mName] == 'direct':
            self.agg_direct(*args)
        elif self.methods_d[mName] == 'filter':
            self.agg_filter(*args)
        else:
            raise QgsProcessingException(f'unrecognized method {mName}')

        # To run another Processing algorithm as part of this algorithm, you can use
        # processing.run(...). Make sure you pass the current context and feedback
        # to processing.run to ensure that all temporary layer outputs are available
        # to the executed algorithm, and that the executed algorithm can send feedback
        # reports to the user (and correctly handle cancellation and progress reports!)
        #=======================================================================
        # if False:
        #     buffered_layer = processing.run("native:buffer", {
        #         'INPUT': dest_id,
        #         'DISTANCE': 1.5,
        #         'SEGMENTS': 5,
        #         'END_CAP_STYLE': 0,
        #         'JOIN_STYLE': 0,
        #         'MITER_LIMIT': 2,
        #         'DISSOLVE': False,
        #         'OUTPUT': 'memory:'
        #     }, context=context, feedback=feedback)['OUTPUT']
        #=======================================================================

        # Return the results of the algorithm. In this case our only result is
        # the feature sink which contains the processed features, but some
        # algorithms may return multiple feature sinks, calculated numeric
        # statistics, etc. These should all be included in the returned
        # dictionary, with keys matching the feature corresponding parameter
        # or output names.
        return {self.OUTPUT_DEM: 'output_dem'}
    
    def agg_direct(self, dem1, wsh1, wse1, scale):
        """WSH Averaging method"""
        
    
    def agg_filter(self, dem1, wsh1, wse1, scale):
        """WSH Averaging method"""
    
#===============================================================================
# HELPERS
#===============================================================================
"""cant do imports without creating a provider and a plugin"""
def assert_extent_equal(left, right, msg='',): 
    """ extents check"""
    if not __debug__: # true if Python was not started with an -O option
        return
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
        
def assert_rlay_simple(rlay, msg='',): 
    """square pixels with integer size"""
    if not __debug__: # true if Python was not started with an -O option
        return
 
 
    
    __tracebackhide__ = True  
    
    x = rlay.rasterUnitsPerPixelX()
    y = rlay.rasterUnitsPerPixelY()
    
    if not x==y:
        raise QgsProcessingException('non-square pixels\n' + msg)
 
    if not round(x, 10)==int(x):
        raise QgsProcessingException('non-integer pixel size\n' + msg)
