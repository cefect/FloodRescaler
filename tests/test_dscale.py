'''
Created on Oct. 17, 2023

@author: cefect
'''


import pytest, copy, os, gc
from qgis.core import (
    QgsRasterLayer, QgsProject,
    QgsProcessingOutputLayerDefinition, QgsApplication
    )

from tests.conftest import get_qrlay
from definitions import wbt_exe, src_dir
from floodrescaler.processing.dscale import Dscale

#NOTE: conftest uses the examples directory
test_data_dir2 = os.path.join(src_dir, 'tests', 'data')


@pytest.fixture(scope='function')
def output_params(qproj, tmpdir):
    """setting up default output parameters for tests"""   
    
    
    def get_out(k):
        """this deletes on exit?""" 
        #return QgsProcessingOutputLayerDefinition(sink='TEMPORARY_OUTPUT', destinationProject=qproj)
        
        """save to pytest dir"""
        return QgsProcessingOutputLayerDefinition(sink=os.path.join(tmpdir, k+'.tif'), destinationProject=qproj)    
    
    
    return {k:get_out(k) for k in ['OUTPUT_WSE']}
 

@pytest.mark.dev
@pytest.mark.parametrize('caseName',['Ahr2021'])
@pytest.mark.parametrize('method, skwargs',[
    ('Resample', dict()),
    ('TerrainFilter', dict()),
    ('CostGrow', dict(filter_isolated_kwargs=dict(method='area'))),
    ('CostGrow', dict(filter_isolated_kwargs=dict(method='pixel'))),   
    ])
def test_runner(dem_layer, wse_layer,   method, skwargs,
                output_params, context, feedback, wbt_init, tmpdir):
    """test the main runner""" 
    assert isinstance(dem_layer, QgsRasterLayer)
    
    #execute
    algo=Dscale()
    algo.initAlgorithm()
    algo._init_algo(output_params, context, feedback, temp_dir=tmpdir)
    
 
    
    res_d = algo.run_dscale( dem_layer, wse_layer, method, **skwargs)
    
    #validate
    assert isinstance(res_d, dict)
    assert set(res_d.keys()).symmetric_difference(output_params.keys())==set()
    
    #todo: add quantiative validation
    

#===============================================================================
# CostGrow-------
#===============================================================================
  
@pytest.mark.parametrize('caseName',['Ahr2021'])
def test_costdistance(wse_fp,  context, feedback, wbt_init):
    """CostDistance WBT extrapolation test""" 
 
    #setup
    algo =  Dscale()
    algo.initAlgorithm()
    algo._init_algo({}, context, feedback) 
    

    
    #execute    
    algo._costdistance(wse_fp )


tdir = lambda x: os.path.join(test_data_dir2, 'test_filter_isolated', x)



@pytest.mark.parametrize('wse_fp, wse_raw_fp',[(tdir('05bdem_mask.tif'), tdir('wse_raw.tif'))])
@pytest.mark.parametrize('method',[
    'area', 
    'pixel',
    ])
@pytest.mark.parametrize('clump_cnt', [3])
def test_filter_isolated(wse_fp, wse_raw_fp, method,  clump_cnt, 
                         context, feedback, wbt_init):
    
    #setup
    algo =  Dscale()
    algo.initAlgorithm()
    algo._init_algo({}, context, feedback) 
    

    
    #execute    
    algo._filter_isolated(wse_fp, wse_raw=QgsRasterLayer(wse_raw_fp), method=method, clump_cnt=clump_cnt)
    


#===============================================================================
# ISSUES----
#===============================================================================
idir = lambda x: os.path.join(src_dir, 'issues', x)


@pytest.mark.parametrize('method, skwargs',[
    ('Resample', dict()),
    ('TerrainFilter', dict()),
    ('CostGrow', dict(filter_isolated_kwargs=dict(method='area'))),
    ('CostGrow', dict(filter_isolated_kwargs=dict(method='pixel'))),   
    ])
@pytest.mark.parametrize('fp_d',
                          [
                              #isusue5
                              #{'dem':idir(r'05\Final_DEM.tif'),'wse':idir(r'05\Final_WSE.tif')},
                              
                              #issue13: filter_isolated
                              {'dem':idir(r'13\input_1.tif'),'wse':idir(r'13\input_2.tif')}
                          ]
                          )
def test_issues(fp_d, method, skwargs,
                output_params, context, feedback, wbt_init, tmpdir):
    """method for debugging tests""" 
    print('testing issues')
    
    #execute
    algo=Dscale()
    algo.initAlgorithm()
    algo._init_algo(output_params, context, feedback, temp_dir=tmpdir)
    res_d = algo.run_dscale( 
        get_qrlay(fp_d['dem']),
        get_qrlay(fp_d['wse']), 
        method, **skwargs)
    
    print(f'temporary directory:\n  {algo.temp_dir}')
    print(f'OUTPUT_WSE:\n    %s'%res_d['OUTPUT_WSE'])
    
    #validate
    assert isinstance(res_d, dict)
    assert set(res_d.keys()).symmetric_difference(output_params.keys())==set()
    
    assert os.path.exists(res_d['OUTPUT_WSE'])
    
    
    