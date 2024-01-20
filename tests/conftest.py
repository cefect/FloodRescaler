'''
Created on Oct. 24, 2022

@author: cefect
'''
#===============================================================================
# IMPORTS------
#===============================================================================
import os, pathlib, pytest, logging, sys
from pytest_qgis.utils import clean_qgis_layer
from qgis.core import (
    QgsRasterLayer, QgsProject, QgsProcessingFeedback, QgsProcessingContext, Qgis, QgsSettings, QgsApplication
    )

#pandas settings
import pandas as pd

# Set the maximum number of rows to 5
pd.set_option('display.max_rows', 7)
pd.options.mode.chained_assignment = None   #setting with copy warning handling
 
from definitions import src_dir, wbt_exe

 
#print(u'QGIS version: %s, release: %s'%(Qgis.QGIS_VERSION.encode('utf-8'), Qgis.QGIS_RELEASE_NAME.encode('utf-8')))

 
test_data_dir = os.path.join(src_dir, 'examples')

assert os.path.exists(test_data_dir)






#===============================================================================
# helpers
#===============================================================================
 
def get_rlay(caseName, layName):
    fp = os.path.join(test_data_dir, caseName, layName+'.tif')
    assert os.path.exists(fp), layName
    
    return QgsRasterLayer(fp, f'{caseName}_{layName}')

def get_qrlay(fp):
    return QgsRasterLayer(fp, os.path.basename(fp).replace('.tif',''))


class MyFeedBackQ(QgsProcessingFeedback):
    """special feedback object for testing"""
    
    def __init__(self,logger, *args, **kwargs):        
        self.logger=logger.getChild('FeedBack')        
        super().__init__(*args, **kwargs)
        
    def pushInfo(self, msg):
        #print(msg)
        self.logger.info(msg)
        
    def pushDebugInfo(self, msg):
        self.logger.debug(msg)
    
#===============================================================================
# fixtures
#===============================================================================

@pytest.fixture(scope='session')
def logger():
    logging.basicConfig(
                #filename='xCurve.log', #basicConfig can only do file or stream
                force=True, #overwrite root handlers
                stream=sys.stdout, #send to stdout (supports colors)
                level=logging.DEBUG, #lowest level to display
                )
    
    return logging.getLogger('r')
    


@pytest.fixture(scope='session')
def qproj(qgis_app, qgis_processing):

    #===========================================================================
    # searchTerm='wbt'
    # for alg in qgis_app.processingRegistry().algorithms():
    #     if searchTerm in alg.id() or searchTerm in alg.displayName():
    #         print(alg.id(), "->", alg.displayName())
    #===========================================================================
    """
    from qgis import processing
    
    processing.run('wbt:ConvertNodataToZero', { 'input' : r'L:\09_REPOS\03_TOOLS\FloodRescaler\examples\Ahr2021\wse.tif', 'output' : 'TEMPORARY_OUTPUT' })
    """
 
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
    """loads the wse from teh test_data_dir (examples)"""
    return get_rlay(caseName, 'wse')


@pytest.fixture(scope='function')
def wse_fp(caseName):
    layName='wse'
    fp = os.path.join(test_data_dir, caseName, layName+'.tif')
    assert os.path.exists(fp), layName
    
    return fp


@pytest.fixture(scope='session')
def wbt_init(qgis_app, qgis_processing):
    """initilze the WhiteBoxTools processing provider"""
    
    #load the provider
    from wbt_for_qgis.wbtprovider import WbtProvider #be sure to add the profile folder to sys.path
    whitebox_provider = WbtProvider()
    
    #add to the registry
    assert qgis_app.processingRegistry().addProvider(whitebox_provider)
    
    #add the exe setting
    from processing.core.ProcessingConfig import ProcessingConfig
    wbt_exe = r'l:\06_SOFT\whitebox\v2.2.0\whitebox_tools.exe'    
    ProcessingConfig.setSettingValue('WBT_EXECUTABLE', wbt_exe)
    
    return whitebox_provider #needs to be held somewhere
    
 
    
