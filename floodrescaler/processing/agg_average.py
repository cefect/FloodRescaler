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
import pprint
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
 
                       )

from qgis.analysis import QgsNativeAlgorithms, QgsRasterCalculatorEntry, QgsRasterCalculator
 
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
            """Aggregating (upscalig) flood layers with methods described in Bryant et. al., (2022). Two methods are provided, both use simple averaging to compute the DEM grid.
             \'WSH Averaging\': usees the WSH (water depths) grid as the primary layer, which is aggregated through simple averaging, while the WSE grid is computed via addition with the DEM grid.
             \'WSE Averaging\' first aggregates the WSE grid, before filtering against the DEM, then subtracting to compute the WSH grid.            
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
            QgsProcessingParameterNumber(self.INPUT_UPSCALE, self.tr('upscaling value'), defaultValue=2.0)
        )
        
        # add methods parameter with some hints
        param = QgsProcessingParameterString(self.INPUT_METHOD, self.tr('aggregation method'), defaultValue='WSH Averaging', multiLine=False)
        param.setMetadata({'widget_wrapper':{ 'value_hints': list(self.methods_d.keys()) }})
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
            inrlay =  self.parameterAsRasterLayer(params, getattr(self, attn), context)
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
        
        #TODO: check null status

        feedback.pushInfo('finished loading raster layers')
        
        #=======================================================================
        # params
        #=======================================================================
        mName = self.parameterAsString(params, self.INPUT_METHOD, context)
        scale = self.parameterAsDouble(params, self.INPUT_UPSCALE, context)
        
        feedback.pushInfo(f'running with \'{mName}\' and scale={scale:.2f}')
        
        #=======================================================================
        # output rasters
        #=======================================================================
        #=======================================================================
        # OUTPUT_DEM = 'OUTPUT_DEM'
        # OUTPUT_WSE = 'OUTPUT_WSE'
        # OUTPUT_WSH = 'OUTPUT_WSH'
        # def get_orlay(attn):
        #     return self.parameterAsOutputLayer(parameters, getattr(self, attn), context)
        #     
        # output_dem = get_orlay('OUTPUT_DEM')
        # output_wsh = get_orlay('OUTPUT_WSH')
        # output_wse = get_orlay('OUTPUT_WSE')
        #=======================================================================
        
        if feedback.isCanceled():
            return {}
        #=======================================================================
        # exceute specified method
        #=======================================================================.
        args = (params, input_dem, input_wsh,  input_wse, scale, context, feedback)
        if self.methods_d[mName] == 'direct':
            res_d = self.agg_direct(*args)       
        
        elif self.methods_d[mName] == 'filter':
            res_d = self.agg_filter(*args)
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
        return res_d
    
    
    def agg_simple(self, inp, scale, 
                   context=None, feedback=None,
                   output='TEMPORARY_LAYER', 
                   ):
        """apply gdal warp"""
        
        feedback.pushInfo(f'gdalwarp on {inp}')
        pars_d = { 'DATA_TYPE' : 0, 'EXTRA' : '', 
                  'INPUT' : inp, 
                  'OUTPUT' : output, #only matters if a filepath si specified? 
                  'MULTITHREADING' : False, 
                  'NODATA' : -9999, 
                  'OPTIONS' : '',                  
                  'RESAMPLING' : 0,  #nearest
                  'SOURCE_CRS' : None, 'TARGET_CRS' : None, 'TARGET_EXTENT' : None, 'TARGET_EXTENT_CRS' : None, 
                  'TARGET_RESOLUTION' : scale}
        
        #gives a filepath regardless
        return processing.run('gdal:warpreproject', pars_d, 
                              context=context, feedback=feedback, is_child_algorithm=True)['OUTPUT']
                              
    def gdal_calc_1(self, inp, formula, output='TEMPORARY_OUTPUT', **kwargs):
        """run gdal raster calc on a single raster. useful for masking operations"""
        
        pars_d = { 
            'INPUT_A' : inp, 'BAND_A' : 1, 'FORMULA':formula,  
            #'FORMULA' : '(A!=0)/(A!=0)',           
            'NO_DATA' : -9999,  'OUTPUT' : output, 'RTYPE' : 5, 
            #'BAND_B' : None, 'BAND_C' : None, 'BAND_D' : None, 'BAND_E' : None, 'BAND_F' : None, 
            #'INPUT_B' : None, 'INPUT_C' : None, 'INPUT_D' : None, 'INPUT_E' : None, 'INPUT_F' : None, 'EXTRA' : '',  'OPTIONS' : '',
                   }
        
        return self._gdal_calc(pars_d, **kwargs)
                              
    def gdal_calc_add(self, inpA, inpB, output='TEMPORARY_OUTPUT', **kwargs):
        """add two rasters together"""
        
        pars_d = { 
            'INPUT_A' : inpA, 'BAND_A' : 1, 'INPUT_B' : inpB, 'BAND_B' : 1, 
            'FORMULA':'A+B',            
            'NO_DATA' : -9999,  'OUTPUT' : output, 'RTYPE' : 5, 
            #'BAND_B' : None, 'BAND_C' : None, 'BAND_D' : None, 'BAND_E' : None, 'BAND_F' : None, 
            #'INPUT_B' : None, 'INPUT_C' : None, 'INPUT_D' : None, 'INPUT_E' : None, 'INPUT_F' : None, 'EXTRA' : '',  'OPTIONS' : '',
                   }
        
        return self._gdal_calc(pars_d, **kwargs)
                              
    def _gdal_calc(self, pars_d, context=None, feedback=None):
        """
        help(processing.run)
        help(processing)
        
        help(processing.runAndLoadResults)
        """
        return processing.run('gdal:rastercalculator', pars_d, 
                              context=context, feedback=feedback, is_child_algorithm=True)['OUTPUT']
    
    def agg_direct(self, params, dem1, wsh1, wse1, scale, context=None, feedback=None):
        """WSH Averaging method"""
        cf_kwargs = dict(context=context, feedback=feedback)    
        
 
        # simple DEM aggregate
        dem2_fp = self.agg_simple(dem1, scale,  
                                  #output=params['OUTPUT_DEM'],
                                   **cf_kwargs)
        
        #simple WSH aggregate
        wsh2_fp = self.agg_simple(wsh1, scale,  output=params['OUTPUT_DEM'], **cf_kwargs)
        
        if feedback.isCanceled():
            return {} 
        
        #=======================================================================
        # #compute wse
        #=======================================================================
        feedback.pushInfo('computing WSE')
        #mask out zeros
        wsh2_maskd_fp = self.gdal_calc_1(wsh2_fp,'A*(A!=0)/(A!=0)', **cf_kwargs)
        
        #DEm + WSH
        wse2_fp = self.gdal_calc_add(wsh2_maskd_fp, dem2_fp, output=params['OUTPUT_WSE'], **cf_kwargs)
         
        
        
        return {self.OUTPUT_DEM:dem2_fp, self.OUTPUT_WSH:wsh2_fp, self.OUTPUT_WSE:wse2_fp}
        
    
    def agg_filter(self, dem1, wsh1, wse1, scale, context, feedback):
        """WSH Averaging method"""
        
        
    def __enter__(self,*args,**kwargs):
        return self
    
#===============================================================================
# HELPERS
#===============================================================================
"""cant do imports without creating a provider and a plugin"""


class RasterCalc(object):
    
    def __init__(self, ref_lay, feedback=None):
        assert isinstance(ref_lay, QgsRasterLayer)
        self.ref_lay = ref_lay
        
        
        
        self.feedback = feedback
        self.rasterEntries = list()
        
        self.ref_rcentry = self.add_rcentry(ref_lay)
        
    def add_rcentry(self, rlay, bandNumber=1):
        """add a QgsRasterCalculatorEntry"""
 
        #=======================================================================
        # check
        #=======================================================================
        assert isinstance(rlay, QgsRasterLayer)
        
        #=======================================================================
        # build the entry
        #=======================================================================
        rcentry = QgsRasterCalculatorEntry()
        rcentry.raster = rlay  # not accesesible for some reason
        rcentry.ref = '%s@%i' % (rlay.name(), bandNumber)
        rcentry.bandNumber = bandNumber
        
        self.rasterEntries.append(rcentry)
        return rcentry
    
    def rcalc(self, formula, output,
              ref_lay=None, rasterEntries=None):
        #=======================================================================
        # defaults
        #=======================================================================
        if ref_lay is None: ref_lay = self.ref_lay
        if rasterEntries is None: rasterEntries = self.rasterEntries
        
        assert isinstance(rasterEntries, list)
        #=======================================================================
        # assemble parameters
        #=======================================================================

        d = dict(
                formula=formula,
                ofp=output,
                outputExtent=ref_lay.extent(),
                outputFormat='GTiff',
                crs=self.qproj.crs(),
                nOutputColumns=ref_lay.width(),
                nOutputRows=ref_lay.height(),
                crsTrnsf=QgsCoordinateTransformContext(),
                rasterEntries=rasterEntries,
            )
        #=======================================================================
        # execute
        #=======================================================================
        
        rcalc = QgsRasterCalculator(d['formula'], d['ofp'],
                                     d['outputFormat'], d['outputExtent'], d['crs'],
                                     d['nOutputColumns'], d['nOutputRows'], d['rasterEntries'], d['crsTrnsf'])
        
        try:
            result = rcalc.processCalculation(feedback=self.feedback)
        except Exception as e:
            raise QgsProcessingException('failed to processCalculation w/ \n    %s' % e)
    
        #=======================================================================
        # check    
        #=======================================================================
        if not result == 0:
            raise QgsProcessingException('formula=%s failed w/ \n    %s' % (formula, rcalc.lastError()))
        
        self.feedback.pushInfo(f'finished on {output}')
        
        return output
    
    def __enter__(self, *args, **kwargs):
        return self


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
