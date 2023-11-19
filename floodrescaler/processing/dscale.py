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
import pprint, os, datetime, tempfile
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
                       QgsMapLayerStore,
                       QgsProcessingOutputLayerDefinition,
                       QgsProject
 
                       )

from qgis.analysis import QgsNativeAlgorithms, QgsRasterCalculatorEntry, QgsRasterCalculator

import pandas as pd
 
class Dscale(QgsProcessingAlgorithm):
    """Downscaling port to QGIS processing script
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
            """Downscaling coarse WSE to a finer resolution with methods described in Bryant et. al., (2023) [10.5194/hess-2023-156] which apply simple hydraulic assumptions and a fine resolution DEM.
            Three methods are provided with varying treatment of the three resample cases (wet-wet (WW), wet-partial (WP), and dry-partial (DP)):
            \'Resample\': [WW] simple bilinear resampling (e.g., gdal_warp)
            \'TerrainFilter\': [WW+WP]: bilinear resampling followed by a terrain filter
            \'CostGrow\': [WW+WP+DP]: as in \'TerrainFilter\' followed by a cost distance routine to grow the inundation and an isolated flooding filter. 
            
        
        \'CostGrow\' requires the \'WhiteboxTools for QGIS\' plugin to be installed and configured (https://www.whiteboxgeo.com/manual/wbt_book/qgis_plugin.html)
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
            'CostGrow':self.run_CostGrow, 
            'Resample':self.run_Resample, 
            'TerrainFilter':self.run_TerrainFilter, 
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
            obj(self.OUTPUT_WSE, self.tr('WSE (downscaled)'))
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
 
 
        return self.run_dscale(input_dem, input_wse, mName, 
                               #context=context, feedback=feedback,
                               **params)
        
    def _init_algo(self, params, context, feedback):
        """common init for tests"""
        self.proc_kwargs = dict(feedback=feedback, context=context, is_child_algorithm=True)
        self.context, self.feedback, self.params = context, feedback, params
        self.temp_dir=tempfile.TemporaryDirectory()

        
        
    def _get_out(self, attn):
        return self.parameterAsOutputLayer(self.params, getattr(self, attn), self.context)
        
        
    
    def run_dscale(self,   dem1_rlay,wse2_rlay, method,
                   **kwargs):
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
        feedback=self.feedback
        context=self.context
        
        #=======================================================================
        # precheck
        #=======================================================================
        assert_extent_equal(dem1_rlay, wse2_rlay)
        #HydTypes('DEM').assert_fp(dem1_fp)
        #HydTypes('WSE').assert_fp(wse2_fp)
        
        #=======================================================================
        # get downscaling ratio
        #=======================================================================
        self.downscale = get_resolution_ratio(dem1_rlay, wse2_rlay)
        
        feedback.pushInfo(f'started w/ \'{method}\' and downscale={self.downscale:.2f}')
        
        #=======================================================================
        # executre
        #=======================================================================
        result= self.methods_d[method](dem1_rlay,wse2_rlay, **kwargs)
    
        #=======================================================================
        #wrap
        #=======================================================================

        feedback.pushInfo(f'finished with temp_dir\n    {self.temp_dir.name}')


        return result

        
        



    def run_Resample(self,_,wse2_rlay,downscale=None,**kwargs):
        """simple Resamlper. wrapper around gdal_warp        
        
        """
        if downscale is None: downscale=self.downscale
        

                
        ofp = self._gdal_warp(wse2_rlay, downscale, OUTPUT=self._get_out(self.OUTPUT_WSE), **kwargs)        
        return {self.OUTPUT_WSE:ofp}
    

    def _dem_filter(self, dem, wse2_fp, OUTPUT='TEMPORARY_OUTPUT'):
 
        return self._gdal_calc({
                'FORMULA':'A*(A > B)', 
                'INPUT_A':wse2_fp, 'BAND_A':1, 
                'INPUT_B':dem, 'BAND_B':1, 
                #'FORMULA' : '(A!=0)/(A!=0)',
                'NO_DATA':0.0, 
                'OUTPUT':OUTPUT, 'RTYPE':5})

    def run_TerrainFilter(self, dem, wse, downscale=None, 
 
                           **kwargs):
        """resampler with terrain filter"""
         
        if downscale is None: downscale=self.downscale
        feedback=self.feedback
        #=======================================================================
        # resample 
        #======================================================================= 
        feedback.pushInfo(f'resample w/ downscale={downscale:.2f}')       
        wse2_fp = self._gdal_warp(wse, downscale) 
        
        if feedback.isCanceled():
            return {} 
        #=======================================================================
        # filter
        #=======================================================================
        feedback.pushInfo('running filter against dem')  
        ofp = self._dem_filter(dem, wse2_fp, OUTPUT=self._get_out(self.OUTPUT_WSE))
        
        #=======================================================================
        # warp
        #=======================================================================
        return {self.OUTPUT_WSE:ofp}
        
        
        
    def run_CostGrow(self, dem, wse, **kwargs):
        """costGrow"""
        
        downscale=self.downscale
        feedback=self.feedback
        
        #=======================================================================
        # 01resample 
        #======================================================================= 
        feedback.pushInfo(f'\n\nresample w/ downscale={downscale:.2f}==============================================')       
        wse2_fp = self._gdal_warp(wse, downscale) 
        
        if feedback.isCanceled():
            return {} 
        #=======================================================================
        # 02filter
        #=======================================================================
        feedback.pushInfo('\n\nDEM filter 1==============================================')  
        wse3_fp = self._dem_filter(dem, wse2_fp,
                                   OUTPUT=os.path.join(self.temp_dir.name, '02dem_filter2.tif'),
                                   )
        
        #=======================================================================
        # 03costDistanceGrow_wbt
        #=======================================================================
        feedback.pushInfo('\n\ncostdistance extrapolation==============================================') 
        wse4_fp = self._costdistance(wse3_fp,
                                     OUTPUT=os.path.join(self.temp_dir.name, '03costdistance.tif'),
                                     )
        
        #=======================================================================
        # 04DEM filter again
        #=======================================================================
        feedback.pushInfo('\n\nDEM filter 2==============================================') 
        wse5_fp = self._dem_filter(dem, wse4_fp, 
                                   OUTPUT=os.path.join(self.temp_dir.name, '04dem_filter2.tif')
                                   )
        
        #=======================================================================
        # 05isolated filter
        #=======================================================================
        feedback.pushInfo('\n\nisolated filter==============================================')
        ofp= self._filter_isolated(wse5_fp, OUTPUT=self._get_out(self.OUTPUT_WSE))
        
        
        #=======================================================================
        # wrap
        #=======================================================================
        return {self.OUTPUT_WSE:ofp}
        
    def _filter_isolated(self, wse_fp, OUTPUT='TEMPORARY_OUTPUT', **kwargs):
        """remove isolated cells from grid using WBT"""
        feedback=self.feedback
        temp_dir = os.path.join(self.temp_dir.name, 'filter_isolated')
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        tfp =lambda x:self._tfp(prefix=x, dir=temp_dir)
 
        #=======================================================================
        # #convert to mask
        #=======================================================================
        """wbt.clump needs 1s and 0=dry"""
        feedback.pushInfo(f'\nprocessing 1=wet 0=dry WSE mask for WBT')
        mask1_fp = self._gdal_calc_mask(wse_fp) 
        
        mask2_fp = processing.run('native:fillnodata', 
                                  { 'BAND' : 1, 'FILL_VALUE' : 0, 
                                   'INPUT' : mask1_fp, 'OUTPUT' : tfp('01_WSEmask_') }, 
                                  **self.proc_kwargs)['OUTPUT']
        
        #=======================================================================
        # clump
        #=======================================================================
        feedback.pushInfo(f'\nprocessing clump mosaic from WSE mask')
        clump_fp = processing.run('wbt:Clump', 
                      { 'diag' : False, 'input' : mask2_fp, 'output' : tfp('02_clumps_'), 'zero_back' : True }, 
                                  **self.proc_kwargs)['output']
        
        #set the nodata value (needed because of the zero_back option)
        clump_fp2 = processing.run('wbt:SetNodataValue', 
                      { 'back_value' : 0, 'input' : clump_fp, 'output' : tfp('02b_clumps2_') }, 
                                  **self.proc_kwargs)['output']


        #identify main clump 
        feedback.pushInfo(f'\nidentifying main clump')       
        unique_fp = processing.run('native:rasterlayeruniquevaluesreport', 
                      { 'BAND' : 1, 'INPUT' : clump_fp2, 
                       #'OUTPUT_HTML_FILE' : 'TEMPORARY_OUTPUT',
                        'OUTPUT_TABLE' : self._tfp(dir=temp_dir, suffix='.csv') }, **self.proc_kwargs)['OUTPUT_TABLE']

                        
        df = pd.read_csv(unique_fp).sort_values('count', ascending=False)
        assert df['value'].min()>0.0, f'failed to filter zerodata'
        ser = df.iloc[0, :]#take largest
        assert ser['count']>0
        max_clump_id = ser['value']
        
        
        
        #build mask of this (1= main clump, 0=not)
        feedback.pushInfo(f'\nidentified largest clump {max_clump_id}\n{ser.to_dict()}\nbuilding mask')
        clump_mask_fp = self._gdal_calc({'FORMULA':f'A=={int(max_clump_id)}', 
                                         'INPUT_A':clump_fp, 'BAND_A':1,'NO_DATA':0,'RTYPE':5, 
                                         'OUTPUT':tfp('03_clumpMask')})
        
        #apply clump filter
        feedback.pushInfo(f'\napplying clump filter from mask on WSE:\n    mask:{clump_mask_fp}\n    wse:{wse_fp}')
        return self._gdal_calc({'FORMULA':f'A/B', 
                                'INPUT_A':wse_fp, 'BAND_A':1,
                                'INPUT_B':clump_mask_fp, 'BAND_B':1,
                                'NO_DATA':0,'RTYPE':5, 'OUTPUT':OUTPUT})
        
        
        
    
    def _costdistance(self, wse_fp, OUTPUT=None,
                      **kwargs):
        """build costallocation WSE extrapolation using WBT"""
        feedback=self.feedback
        context=self.context
        feedback.pushInfo(f'costDistance + costAllocation on {wse_fp}')

        if OUTPUT is None:
            OUTPUT=self._tfp()
        
        #=======================================================================
        # costDistance
        #=======================================================================
 
        """
        couldnt get the native temp files to work
 
        """
        #fillnodata in wse (for source)
        wse1_fp = processing.run('wbt:ConvertNodataToZero', 
                                 { 'input' : wse_fp, 'output' : self._tfp()}, 
                                 **self.proc_kwargs)['output']
        
        #build cost friction (constant)
        cost_fric_fp = processing.run('wbt:NewRasterFromBase',
                                      { 'base' : wse1_fp,  'data_type' : 1, 'output' : self._tfp(), 'value' : 1.0}, 
                                      **self.proc_kwargs)['output']
        
        #compute backlink raster
        backlink_fp = processing.run('wbt:CostDistance',
                                      { 'cost' : cost_fric_fp, 'out_accum' : self._tfp(), 
                                       'out_backlink' : self._tfp(), 'source' : wse1_fp },
                                      **self.proc_kwargs)['out_backlink']
        
        feedback.pushInfo(f'built costDistance backlink raster \n    {backlink_fp}')
        
        #=======================================================================
        # costAllocation
        #=======================================================================
        costAlloc_fp = processing.run('wbt:CostAllocation',
                                      { 'backlink' : backlink_fp, 
                                       'output' : OUTPUT, 
                                       'source' : wse1_fp},
                                      **self.proc_kwargs)['output']
                                      
        feedback.pushInfo(f'built CostAllocation \n    {costAlloc_fp}')
        
        return costAlloc_fp
        
         
        
    def _gdal_warp(self, wse2_rlay, downscale, OUTPUT='TEMPORARY_OUTPUT', RESAMPLING=1, **kwargs):
        """convenience for gdal warp
        
        
        not sure how to set the layer name"""
        
        target_res = rlay_get_resolution(wse2_rlay) / downscale
        #feedback.pushInfo(f'gdalwarp on {wse2_rlay}')
        pars_d = {'DATA_TYPE':0, 'EXTRA':'', 
            'INPUT':wse2_rlay, 
            'OUTPUT':OUTPUT, 
            'MULTITHREADING':False, 
            'NODATA':-9999, 
            'OPTIONS':'', 
            'RESAMPLING':RESAMPLING, #bilinear
            'SOURCE_CRS':None, 
            'TARGET_CRS':None, 'TARGET_EXTENT':None, 'TARGET_EXTENT_CRS':None, 
            'TARGET_RESOLUTION':target_res}
        
        #gives a filepath regardless
        ofp = processing.run('gdal:warpreproject', pars_d,  **self.proc_kwargs)['OUTPUT']
        if not os.path.exists(ofp):
            raise QgsProcessingException('gdal:warpreproject failed to get a result for \n%s' % pars_d['INPUT'])
        return ofp
    
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
        return self._gdal_calc({'FORMULA':'A/A', 'INPUT_A':rlay, 'BAND_A':1,'NO_DATA':-9999,'OUTPUT':OUTPUT, 'RTYPE':5})

    def _tfp(self, prefix=None, suffix='.tif', dir=None):
        if dir is None: dir=self.temp_dir.name
        return tempfile.NamedTemporaryFile(suffix=suffix,prefix=prefix, dir=dir).name
        
 

    
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


