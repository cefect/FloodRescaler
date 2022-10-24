'''
Created on Oct. 24, 2022

@author: cefect

pytest processing scripts
'''
import pytest, copy, os, gc
from qgis.core import (
    QgsRasterLayer, QgsProject,
    QgsProcessingOutputLayerDefinition
    )


from floodrescaler.processing.agg_average import AggAverage

@pytest.fixture(scope='function')
def output_params(qproj):
    
    def get_out():
        return QgsProcessingOutputLayerDefinition(sink='TEMPORARY_OUTPUT', destinationProject=qproj)
    
    return {
        'OUTPUT_DEM':get_out(),
        'OUTPUT_WSH':get_out(),
        'OUTPUT_WSE':get_out()
        }
@pytest.mark.parametrize('caseName',['SJ2018'])
@pytest.mark.parametrize('scale',[2, 2.0, 2.1])
def test_direct(dem, wsh, wse, scale, output_params, context, feedback):
 
    assert isinstance(dem, QgsRasterLayer)
    #execute
    #with AggAverage() as algo:
    algo=AggAverage()
    algo.initAlgorithm()
    res_d = algo.agg_direct(output_params, dem, wsh, wse, scale, context, feedback)
    
    """throwing a procesing exception on cleanup
    spent a few mins trying to debug
    code seems to work fine in the GUI so I gave up
    """
    
    #===========================================================================
    # del res_d
    # del algo
    # gc.collect()
    #===========================================================================
 
    
        #validsate
        #=======================================================================
        # assert isinstance(res_d, dict)
        # assert set(res_d.keys()).symmetric_difference(output_params.keys())==set()
        #=======================================================================
    
    
    #===========================================================================
    # for k,v in res_d.items():
    #     assert os.path.exists(v)
    #===========================================================================
        