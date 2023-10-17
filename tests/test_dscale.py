'''
Created on Oct. 17, 2023

@author: cefect
'''


import pytest, copy, os, gc
from qgis.core import (
    QgsRasterLayer, QgsProject,
    QgsProcessingOutputLayerDefinition
    )


from floodrescaler.processing.dscale import Dscale


@pytest.fixture(scope='function')
def output_params(qproj):
    
    def get_out():
        return QgsProcessingOutputLayerDefinition(sink='TEMPORARY_OUTPUT', destinationProject=qproj)
    
    return {
        'OUTPUT_WSE':get_out()
        }
    
@pytest.mark.dev  
@pytest.mark.parametrize('caseName',['Ahr2021'])
@pytest.mark.parametrize('method',[
    #'CostGrow',
    'Resample'])
def test_runner(dem,wse,   method, output_params, context, feedback):
 
    assert isinstance(dem, QgsRasterLayer)
    #execute
    #with AggAverage() as algo:
    algo=Dscale()
    algo.initAlgorithm()
    res_d = algo.run_dscale( dem, wse, method,
                            context=context, feedback=feedback, **output_params)
    
    #validate
    assert isinstance(res_d, dict)
    assert set(res_d.keys()).symmetric_difference(output_params.keys())==set()
    
    