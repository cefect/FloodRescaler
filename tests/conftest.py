'''
Created on Oct. 24, 2022

@author: cefect
'''
import os, pathlib, pytest, logging, sys
from pytest_qgis.utils import clean_qgis_layer
from qgis.core import (
    QgsRasterLayer, QgsProject, QgsProcessingFeedback, QgsProcessingContext, Qgis
    )
from definitions import src_dir

 
print(u'QGIS version: %s, release: %s'%(Qgis.QGIS_VERSION.encode('utf-8'), Qgis.QGIS_RELEASE_NAME.encode('utf-8')))

 
test_data_dir = os.path.join(src_dir, 'examples')

assert os.path.exists(test_data_dir)

#===============================================================================
# helpers
#===============================================================================
 
def get_rlay(caseName, layName):
    fp = os.path.join(test_data_dir, caseName, layName+'.tif')
    assert os.path.exists(fp), layName
    
    return QgsRasterLayer(fp, f'{caseName}_{layName}')


class MyFeedBackQ(QgsProcessingFeedback):
    """special feedback object for testing"""
    
    def __init__(self,logger, *args, **kwargs):        
        self.logger=logger.getChild('FeedBack')        
        super().__init__(*args, **kwargs)
        
    def pushInfo(self, info):
        self.logger.info(info)
        
    def pushDebugInfo(self, info):
        self.logger.debug(info)
    
#===============================================================================
# fixtures
#===============================================================================

@pytest.fixture(scope='session')
def logger():
    logging.basicConfig(
                #filename='xCurve.log', #basicConfig can only do file or stream
                force=True, #overwrite root handlers
                stream=sys.stdout, #send to stdout (supports colors)
                level=logging.INFO, #lowest level to display
                )
    
    return logging.getLogger('r')
    


@pytest.fixture(scope='session')
def qproj(qgis_app, qgis_processing):
    return QgsProject.instance()

@pytest.fixture(scope='session')
def feedback(qproj, logger):
    return MyFeedBackQ(logger, False)

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
    
