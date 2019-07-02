# ability to read/write fits files
from astropy.io import fits

import numpy as np


def read_image(infile,mosaic=True):
    '''
    read_image
    ----------
    routine to read in an image and either return
    
    mosaic==True
        single mosaicked, non-overscan frame
        
    mosaic==False
        dictionary of numbered quadrants for overscan operations
        
        (can later be combined with arrange_quadrants)
        
        
    inputs
    ----------
    infile         : (string) filename to be read in
    mosaic         : (boolean, default=True) if True, returns a single data frame,
                     if False, returns a dictionary of the four quadrants
    
    outputs
    ----------
    data_quad      : (dictionary or array) if dictionary, keys are [0,1,2,3], each 2056x2048
                     corresponding to each quadrant. if array, single 4122x4096 frame
    
    overscan_quad  : (dictionary) keys are [0,1,2,3], each 2056x2048, corresponding to each quadrant


    dependents
    ----------
    arrange_quadrants : definition to place quadrants in the correct configuration, below

    '''
    
    ofile = fits.open(infile)

    # retreive header for the purposes of sizing the array
    phdr = fits.getheader(infile,0)

    saxis_x = int(phdr['CNAXIS2'])
    saxis_y = int(phdr['CNAXIS1'])+int(phdr['OVRSCAN1'])
        
    # overscan array is solely the overscan region of each quadrant
    overscan_quad = {}
        
    # median_array is the entire array, overscan included
    data_quad = {}
    
    for ext in range(1,5):
        overscan_quad[ext-1] = ofile[ext].data[0:saxis_x,phdr['CNAXIS1']:saxis_y]
        data_quad[ext-1] = ofile[ext].data[0:saxis_x,0:phdr['CNAXIS1']]
        
        
    if mosaic:
        data_mosaic = arrange_quadrants(data_quad)
        
        return data_mosaic,overscan_quad
    
    else:
        return data_quad,overscan_quad
        

        

def arrange_quadrants(quadrants):
    '''
    arrange_quadrants
    -----------------
    rearrange HDI quadrants to be in the proper configuration
    
    can be done with or without overscan, in theory.
    
    inputs
    --------
    quadrants   : (dictionary) dictionary of the four quadrants, with keys [0,1,2,3]
    
    outputs
    --------
    data_array  : (matrix)
    
    '''
    
    saxis_x,saxis_y = quadrants[0].shape
    data_array = np.zeros([2*saxis_x,2*saxis_y])

    # reposition different segments
    data_array[0:saxis_x        ,0:saxis_y]         = quadrants[0]            # lower left
    data_array[0:saxis_x        ,saxis_y:2*saxis_y] = quadrants[1][:,::-1]    # lower right
    data_array[saxis_x:2*saxis_x,0:saxis_y]         = quadrants[2][::-1,:]    # upper left
    data_array[saxis_x:2*saxis_x,saxis_y:2*saxis_y] = quadrants[3][::-1,::-1] # upper right

    # include the left-right flip
    return data_array[:,::-1]



def find_files(directory,img_type,filternum=-1):
    '''
    #
    # definition to go through HDI headers and get the images desired of a particular type
    #    and filter
    #
    
    #
    # limitations:
    #    - not set up to sort by specific targets
    #    - will not enforce one of the filter wheel slots being empty
    # 
    '''
    # grab all HDI files from the specified directory
    files = [infile for infile in glob.glob(directory+'c*t*fits') if not os.path.isdir(infile)]
    
    out_files = []

    for f in files:
        #print(f)
        phdr = fits.getheader(f,0) # the 0 is needed to get the correct extension
    
        # filter by desired_type
        
        # if biases, don't care about filter
        if    (img_type == phdr['OBSTYPE']) \
            & (img_type == 'BIAS'):
                    
            out_files.append(f)
            
        # if flats or objects, the filter matters
        if    (img_type == phdr['OBSTYPE']) \
            & (img_type != 'BIAS')\
            & ( (str(filternum) == phdr['FILTER1']) | (str(filternum) == phdr['FILTER2']) ):

            out_files.append(f)


    return out_files


