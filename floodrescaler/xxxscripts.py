'''
Created on Oct. 24, 2022

@author: cefect

supporting functions
'''

from qgis.core import (QgsRasterLayer,
                       QgsProcessingException,
                       )

def assert_extent_equal(left, right, msg='',): 
    """ extents check"""
    if not __debug__: # true if Python was not started with an -O option
        return
    assert isinstance(left, QgsRasterLayer), type(left).__name__+ '\n'+msg
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
        
def assert_rlay_simple(rlay, msg='',): 
    """square pixels with integer size"""
    if not __debug__: # true if Python was not started with an -O option
        return
 
 
    
    __tracebackhide__ = True  
    
    x = rlay.rasterUnitsPerPixelX()
    y = rlay.rasterUnitsPerPixelY()
    
    if not x==y:
        raise QgsProcessingException('non-square pixels\n' + msg)
 
    if not round(x, 10)==int(x):
        raise QgsProcessingException('non-integer pixel size\n' + msg)