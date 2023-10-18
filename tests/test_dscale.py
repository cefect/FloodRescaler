'''
Created on Oct. 17, 2023

@author: cefect
'''


import pytest, copy, os, gc
from qgis.core import (
    QgsRasterLayer, QgsProject,
    QgsProcessingOutputLayerDefinition, QgsApplication
    )

from definitions import wbt_exe
from floodrescaler.processing.dscale import Dscale




@pytest.fixture(scope='function')
def output_params(qproj):
    """setting up default output parameters for tests"""    
    def get_out():
        return QgsProcessingOutputLayerDefinition(sink='TEMPORARY_OUTPUT', destinationProject=qproj)
    
    return {
        'OUTPUT_WSE':get_out()
        }
    

@pytest.mark.dev 
@pytest.mark.parametrize('caseName',['Ahr2021'])
@pytest.mark.parametrize('method',[
    'Resample',
    'TerrainFilter',
    'CostGrow',    
    ])
def test_runner(dem,wse,   method, output_params, context, feedback, wbt_init):
    """test the main runner""" 
    assert isinstance(dem, QgsRasterLayer)
    
    #execute
    algo=Dscale()
    algo.initAlgorithm()
    algo._init_algo(output_params, context, feedback)
    res_d = algo.run_dscale( dem, wse, method)
    
    #validate
    assert isinstance(res_d, dict)
    assert set(res_d.keys()).symmetric_difference(output_params.keys())==set()
    
    #todo: add quantiative validation
    


  
@pytest.mark.parametrize('caseName',['Ahr2021'])
def test_costdistance(wse_fp,  context, feedback, wbt_init):
    """test the main runner""" 
 
    #setup
    algo =  Dscale()
    algo.initAlgorithm()
    algo._init_algo({}, context, feedback) 
    

    
    #execute    
    algo._costdistance(wse_fp )