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
__version__ = '2024.01.19'
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

#descriptions_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'descriptions')
 
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
    
    
    #parameters
    cg_isolated_methods_l = ['area', 'pixel']
    INPUT_CG_ISO_clump_cnt='INPUT_CG_ISO_clump_cnt'
    INPUT_CG_ISO_METHOD='INPUT_CG_ISO_METHOD'
 
    
    

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
        return 'dscale'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr('Downscaling WSE')

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
            <p>Downscaling coarse WSE to a finer resolution with methods described in <a href="https://doi.org/10.5194/hess-2023-156">Bryant et. al., (2023)</a> which apply simple hydraulic assumptions and a fine resolution DEM.</p>
            <p>Three methods are provided with varying treatment of the three resample cases (wet-wet (WW), wet-partial (WP), and dry-partial (DP)):</p>
            <ul>
                <li><strong>Resample</strong>: [WW] simple bilinear resampling (e.g., gdal_warp)</li>
                <li><strong>TerrainFilter</strong>: [WW+WP]: bilinear resampling followed by a terrain filter</li>
                <li><strong>CostGrow</strong>: [WW+WP+DP]: as in 'TerrainFilter' followed by a cost distance routine to grow the inundation and an isolated flooding filter.</li>
            </ul>
            <p>CostGrow requires the <strong>WhiteboxTools for QGIS</strong> plugin to be installed and configured (<a href="https://www.whiteboxgeo.com/manual/wbt_book/qgis_plugin.html">https://www.whiteboxgeo.com/manual/wbt_book/qgis_plugin.html</a>)</p>
            
            <p>CostGrow's 'isolated flooding' filter is used to select the result WSE clump from the cost distance map. Two methods are available for this:</p>
            <ul>
                <li><strong>area</strong>: the n='selection count' clumps with largest area are selected (faster, may be inaccurate for large complex domains)</li>
                <li><strong>pixel</strong>: those clumps overlapping the original coarse WSE are included (slower, more accurate) </li>
            </ul>
            

            <h3>Tips and Tricks</h3>
            The WSE and DEM grids need to have the same extent
            The WSE resolution should be an even multiple of the WSE resolution (e.g., downscale from 4m to 2m... not from 5m to 2m). This may require some pre-processing (gdal_warp) to obtain an even multiple.
            For the CostGrow method, the 'isolated filter' step is often problematic. To debug this (and other steps), compare the output log to the files in the temporary directory to see which step is causing the error.

            <h3>Issues and Updates</h3>
            See the <a href="https://github.com/cefect/FloodRescaler">project repository</a> for updates, to post an issue/question, or to request new features.

            <h3>Attribution</h3>
            If you use these tools for your work, please cite the following:
            <a href="https://doi.org/10.5194/hess-2023-156">Bryant et. al., (2023)</a>
            <a href="https://jblindsay.github.io/ghrg/pubs/LindsayGISRUK2014.pdf">Lindsay (2021)</a>
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
            QgsProcessingParameterRasterLayer(self.INPUT_DEM, self.tr('DEM (fine)'))
        )
        
        self.addParameter(
            QgsProcessingParameterRasterLayer(self.INPUT_WSE, self.tr('WSE (coarse)'))
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
        
        
        #CostGrow._filter_isolated.methods
        param = QgsProcessingParameterString(self.INPUT_CG_ISO_METHOD, 
                                             self.tr('CostGrow: isolated filter: method'), 
                                             defaultValue=self.cg_isolated_methods_l[0], 
                                             optional=True,
                                             multiLine=False)
        
        param.setMetadata({'widget_wrapper':{ 'value_hints':self.cg_isolated_methods_l }})
        self.addParameter(param)
        
        #CostGrow._filter_isolated.clump_cnt
        param = QgsProcessingParameterNumber(self.INPUT_CG_ISO_clump_cnt, 
                                             self.tr('CostGrow: isolated filter: method=\'area\': selection count'), 
                                             defaultValue=1, optional=True,
                                             type=QgsProcessingParameterNumber.Integer)
        
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
        method = self.parameterAsString(params, self.INPUT_METHOD, context)
        kwargs = dict()
        
        if method=='CostGrow':
            kwargs['filter_isolated_kwargs'] = dict(
                clump_cnt = self.parameterAsInt(params, self.INPUT_CG_ISO_clump_cnt, context),
                method = self.parameterAsString(params, self.INPUT_CG_ISO_METHOD, context)
                )

            
            
 
        
        feedback.pushInfo(f'running with \'{method}\'\n{kwargs}') 
        
        if feedback.isCanceled():
            return {}
        #=======================================================================
        # exceute specified method----
        #=======================================================================.
 
 
        return self.run_dscale(input_dem, input_wse, method, 
                               #context=context, feedback=feedback,
                               **params, **kwargs)
        
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
        
        feedback.pushInfo(f'started w/ \'{method}\' and downscale={self.downscale:.2f} w/ verison {__version__}')
        
        #=======================================================================
        # executre
        #=======================================================================
        result= self.methods_d[method](dem1_rlay,wse2_rlay, **kwargs)
    
        #=======================================================================
        #wrap
        #=======================================================================

        feedback.pushInfo(f'finished with temp_dir\n    {self.temp_dir}')


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
        wse2_fp = self._gdal_warp(wse, downscale,
                                  OUTPUT=self._tfp(prefix='01_resample')) 
        
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
        
        
        
    def run_CostGrow(self, dem, wse, 
                     filter_isolated_kwargs=dict(),
                     **kwargs):
        """costGrow"""
        
        downscale=self.downscale
        feedback=self.feedback
        
        #=======================================================================
        # 01resample 
        #======================================================================= 
        feedback.pushInfo(f'\n\nresample w/ downscale={downscale:.2f} \n==============================================')       
        wse2_fp = self._gdal_warp(wse, downscale, OUTPUT=os.path.join(self.temp_dir, '01resample.tif')) 
        
        if feedback.isCanceled():
            return {} 
        #=======================================================================
        # 02filter
        #=======================================================================
        feedback.pushInfo('\n\nDEM filter 1 \n==============================================')  
        wse3_fp = self._dem_filter(dem, wse2_fp,
                                   OUTPUT=os.path.join(self.temp_dir, '02dem_filter2.tif'),
                                   )
        
        #=======================================================================
        # 03costDistanceGrow_wbt
        #=======================================================================
        feedback.pushInfo('\n\ncostdistance extrapolation \n==============================================') 
        wse4_fp = self._costdistance(wse3_fp,
                                     OUTPUT=os.path.join(self.temp_dir, '03costdistance.tif'),
                                     )
        
        #=======================================================================
        # DECAY
        #=======================================================================
        """placeholder"""
        
        #=======================================================================
        # 04DEM filter again
        #=======================================================================
        feedback.pushInfo('\n\nDEM filter 2==============================================') 
        wse5_fp = self._dem_filter(dem, wse4_fp, 
                                   OUTPUT=os.path.join(self.temp_dir, '04dem_filter2.tif')
                                   )
        

        #apply DEM mask
        wse6_fp=self._gdal_calc({
                'FORMULA':'A*(B/B)', 
                'INPUT_A':wse5_fp, 'BAND_A':1, 
                'INPUT_B':dem, 'BAND_B':1, 
                #'FORMULA' : '(A!=0)/(A!=0)',
                'NO_DATA':0.0, 
                'OUTPUT':os.path.join(self.temp_dir, '05bdem_mask.tif'), 
                'RTYPE':5})
        
        #=======================================================================
        # 05isolated filter
        #=======================================================================
        feedback.pushInfo('\n\nisolated filter \n==============================================')
        ofp= self._filter_isolated(wse6_fp, OUTPUT=self._get_out(self.OUTPUT_WSE), wse_raw = wse,
                                   **filter_isolated_kwargs)
        
        assert os.path.exists(ofp)
        
        #=======================================================================
        # wrap
        #=======================================================================
        feedback.pushDebugInfo(f'finished on \n    {ofp}')
        return {self.OUTPUT_WSE:ofp}
        

    def _get_mask_of_clumps(self, tfp, clump_fp, clump_ids):
        """get a mask from specified pixel values"""
        #reclassify to pull out good values
        reclass_vals = ';'.join([f'9999.0;{e}' for e in clump_ids])
        reclass_fp = processing.run('wbt:Reclass', {'assign_mode':True, 'input':clump_fp, 'output':tfp('03_clumpReclassify'), 'reclass_vals':reclass_vals}, **self.proc_kwargs)['output']
        
 
        clump_mask_fp = self._gdal_calc({'FORMULA':f'A==9999.0', 
                'INPUT_A':reclass_fp, 'BAND_A':1, 'NO_DATA':0, 'RTYPE':5, 
                'OUTPUT':tfp('03_clumpMask')})
        
        return clump_mask_fp

    def _filter_isolated(self, wse_fp,
                         clump_cnt=1,
                        method='area',
                         min_pixel_frac=0.01,
                         wse_raw=None,                         
                         OUTPUT='TEMPORARY_OUTPUT', 
                         **kwargs):
        """remove isolated cells from grid using WBT
        
        
        Params
        -------
        wse_fp: str
            filepath to fine resolution WSE on which to apply the isolated filter
            
        clump_cnt: int
            number of clumps to select
        
        method: str
            method for selecting isolated flood groups:
            'area': take the n=clump_cnt largest areas (fast, only one implemented in FloodRescaler)
            'pixel': use polygons and points to select the groups of interest
            
            
        min_pixel_frac: float
            for method='pixel', the minimum fraction of total domain pixels allowed for the clump search
            
        wse_raw: str
            for method='pixel', raw WSE from which to get intersect points
            ideally the raw coarse WSE input
            see note below on resampling
        
        """
        #=======================================================================
        # setup
        #=======================================================================
        feedback=self.feedback
        temp_dir = os.path.join(self.temp_dir, 'filter_isolated')
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        tfp =lambda x:self._tfp(prefix=x, dir=temp_dir) #get temporary file
        
        log = self.logger
                
 
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

                        
        clump_df = pd.read_csv(unique_fp).sort_values('count', ascending=False)
        assert clump_df['value'].min()>0.0, f'failed to filter zerodata'
        #=======================================================================
        # ser = df.iloc[0, :]#take largest
        # assert ser['count']>0
        # max_clump_id = ser['value']
        #=======================================================================
        
        #===================================================================
        # area-based selection-----
        #===================================================================
        feedback.pushInfo(f'clump selection method \'{method}\'')
 
        if method == 'area':
            feedback.pushInfo(f'area-based selecting {clump_cnt}/{len(clump_df)} clumps')            
            clump_ids = clump_df.iloc[0:clump_cnt, 0].values
        
        
        
        #=======================================================================
        # pixel-based selection-------
        #=======================================================================
        elif method=='pixel':
            assert isinstance(wse_raw, QgsRasterLayer), 'for method=pixel must pass wse_raw'
 
            log.debug(f'filtering clumps w/ min_pixel_frac={min_pixel_frac}')
 
            #===================================================================
            # drop some small clumps
            #===================================================================
            """because polygnoizing is very slow... want to drop any super small groups"""
            #select using 100 or some fraction
            bx = clump_df['count']>min(min_pixel_frac*clump_df['count'].sum(), 100)
 
            
            log.debug(f'pre-selected {bx.sum()}/{len(bx)} clumps for pixel selection')
            
            # build a mask of this
            mask_large_clumps = self._get_mask_of_clumps(tfp, clump_fp, clump_df[bx].iloc[:,0].values)
            
            #apply it to the clumps
            clump_major_fp = self._gdal_calc({'FORMULA':f'A*B', 
                                'INPUT_A':clump_fp2, 'BAND_A':1,
                                'INPUT_B':mask_large_clumps, 'BAND_B':1,
                                'NO_DATA':0,'RTYPE':5, 'OUTPUT':tfp('03_clumpsMajor_')})
            
            log.debug(f'wrote major clumps to \n    {clump_major_fp}')
            
            #===================================================================
            # polygonize clumps
            #===================================================================            
            clump_vlay_fp = processing.run('gdal:polygonize', 
                                 { 'BAND' : 1, 'EIGHT_CONNECTEDNESS' : True, 
                                  'EXTRA' : '', 'FIELD' : 'clump_id', 
                                  'INPUT' : clump_major_fp, 
                                  'OUTPUT' : self._tfp(prefix='03_clumpPoly_', dir=temp_dir, suffix='.shp') },
                                 **self.proc_kwargs)['OUTPUT']
            
 
            log.debug(f'vectorized clumps to \n    {clump_vlay_fp}')
 
            #===================================================================
            # intersect of clumps and wse coarse
            #===================================================================
            """NOTE: could speed things up by coarsening this
            could make more accurate by using polygons (poly v poly intersect)
            """
            #get coarse points
            wse_raw_pts = processing.run('wbt:RasterToVectorPoints', 
                      { 'input' : wse_raw, 'output' : self._tfp(prefix='04_rawWsePoints', dir=temp_dir, suffix='.shp') }, 
                                  **self.proc_kwargs)['output']
 
            
            #select intersecting new polygons 
            clump_vlay_sel_fp = processing.run('native:extractbylocation', 
                                 { 'INPUT' : clump_vlay_fp, 'INTERSECT' : wse_raw_pts, 
                                  'OUTPUT' : tfp('05_clumpPolySelected_'), 'PREDICATE' : [0] },
                                  **self.proc_kwargs)['OUTPUT']
            

            
            #list selection
            uqv_s = processing.run('qgis:listuniquevalues', 
                                  { 'FIELDS' : ['clump_id'],
                                    'INPUT' : QgsProcessingFeatureSourceDefinition(clump_vlay_sel_fp, 
                                                                                   selectedFeaturesOnly=False, 
                                                                                   featureLimit=-1, 
                                                                                   geometryCheck=QgsFeatureRequest.GeometryAbortOnInvalid), 
                                   'OUTPUT' : 'TEMPORARY_OUTPUT', 'OUTPUT_HTML_FILE' : 'TEMPORARY_OUTPUT' },
                                  **self.proc_kwargs)['UNIQUE_VALUES']
            
            assert len(uqv_s)>0, f'failed to get any unique polygons during clump selection'
            clump_ids = list(map(int, uqv_s.split(';')))
 
            log.debug(f'selected {len(clump_ids)}/{len(clump_df)} clumps by pixel intersect w/ {wse_raw}')
            
            
        else:
            raise KeyError(method)
        
        #=======================================================================
        # #build mask of this (1= main clump, 0=not)
        #=======================================================================
        feedback.pushInfo(f'\nidentified {len(clump_ids)} clumps\nbuilding mask')
        
        clump_mask_fp = self._get_mask_of_clumps(tfp, clump_fp, clump_ids)
        
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
        if dir is None: dir=self.temp_dir
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

#===============================================================================
# def gdal_get_array(rasterPath):
#  
#     
#     # Open the raster file
#     ds = gdal.Open(rasterPath)
#     
#     # Get the raster band (as an example, getting the first band of the raster)
#     band1 = ds.GetRasterBand(1)
#     
#     # Read the band as a numpy array
#     array1 = band1.ReadAsArray()
#     
#     del ds
#     
#     return array1
#===============================================================================

class log(object):
    """log emulator"""
    def __init__(self, feedback):
        self.feedback=feedback
    def debug(self, msg):
        self.feedback.pushDebugInfo(msg)
    def info(self, msg):
        self.feedback.pushInfo(msg)

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


