'''
Created on Oct. 24, 2022

@author: cefect
'''
import os, pathlib, pytest
from pytest_qgis.utils import clean_qgis_layer
from qgis.core import (
    QgsRasterLayer, QgsProject, QgsProcessingFeedback, QgsProcessingContext
    )
from definitions import src_dir

 
test_data_dir = os.path.join(src_dir, 'examples')

assert os.path.exists(test_data_dir)

#===============================================================================
# helpers
#===============================================================================
 
def get_rlay(caseName, layName):
    fp = os.path.join(test_data_dir, caseName, layName+'.tif')
    assert os.path.exists(fp), layName
    
    return QgsRasterLayer(fp, f'{caseName}_{layName}')
    
#===============================================================================
# fixtures
#===============================================================================
@pytest.fixture(scope='session')
def qproj(qgis_app, qgis_processing):
    return QgsProject.instance()

@pytest.fixture(scope='session')
def feedback(qproj):
    return QgsProcessingFeedback(False)

@pytest.fixture(scope='session')
def context(qproj):
    return QgsProcessingContext()
 
@pytest.fixture(scope='function')
@clean_qgis_layer
def dem(caseName, qproj):
    return get_rlay(caseName, 'dem')

@pytest.fixture(scope='function')
@clean_qgis_layer
def wsh(caseName, qproj):
    return get_rlay(caseName, 'wsh')


@pytest.fixture(scope='function')
@clean_qgis_layer
def wse(caseName, qproj):
    return get_rlay(caseName, 'wse')
    
