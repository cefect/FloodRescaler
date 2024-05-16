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
from floodrescaler.processing.water_from_complement import WaterFromComp as ProcAlg

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
    
    
    return {k:get_out(k) for k in ['OUTPUT_WATER']}
 

@pytest.mark.dev
@pytest.mark.parametrize('caseName',['Ahr2021'])
@pytest.mark.parametrize('INPUT_TYPE',[
    'WSE', 
    #'WSH'
    ])
def test_runner(dem_coarse_layer, wse_layer,   INPUT_TYPE,
                output_params, context, feedback, tmpdir):
    
    """test the main runner""" 
    assert isinstance(dem_coarse_layer, QgsRasterLayer)
    
    #execute
    algo=ProcAlg()
    algo.initAlgorithm()
    algo._init_algo(output_params, context, feedback, temp_dir=tmpdir)
    
 
    if INPUT_TYPE=='WSE':
        water_rlay = wse_layer
    elif INPUT_TYPE=='WSH':
        water_rlay=wsh_rlay
    else:
        raise KeyError(INPUT_TYPE)
    
    result = algo.run_water_grid_convert( dem_coarse_layer, wse_layer, INPUT_TYPE)
    







   
    