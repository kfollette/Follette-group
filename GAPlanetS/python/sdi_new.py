##module with functions to calculate SDI scale factors from Line/Cont cubes

def radarr(xdim,ydim):
    '''
    make an image of size xdim x ydim where every pixel has a value equal to its distance from the center of the array
    '''
    im = np.zeros((xdim, ydim))
    xcen = (xdim - 1) / 2.
    ycen = (ydim - 1) / 2.
    for index1 in range(xdim):
        for index2 in range(ydim):
            im[index1, index2]=np.sqrt((xcen-index1)**2. + (ycen-index2)**2.)        
    return(im)

def annularmask(xdim, ydim, rin, rout):
    '''
    makes a xdim x ydim mask where pixels with distances from center between rin and
    rout have value 1, and everything else has been replaced by NaNs
    '''
    arr = radarr(xdim, ydim)
    new = np.ones((xdim, ydim))
    for y in range(0,ydim):
        for x in range(0,xdim):
            if arr[y][x] <= rin:
                new[y][x] = np.nan
            elif arr[y][x] >= rout:
                new[y][x] = np.nan
    #hdu=fits.PrimaryHDU(new)
    #hduList = fits.HDUList([hdu])
    #hduList.writeto('ctrlmask_test.fits', overwrite=True)
    return(new)

#calculates total starlight in one image
def totalStarlight(image):
    '''
    masks an image and sums non-nan pixels
    image: 2D array
    requires innerRadius and outerRadius to be defined
    outputs float of total object counts between inner and outer radius (IWA and control radius)
    '''
    dim = image.shape[1]
    mask = annularmask(dim, dim, innerRadius, outerRadius)
    maskedim = image*mask
    return np.nansum(maskedim)

#calculates median starlight in one image
def medianStarlight(image):
    '''
    masks an image and finds the (nan)median of the unmasked pixels
    image: 2D array
    requires innerRadius and outerRadius to be defined
    outputs float of median object counts between inner and outer radius
    '''
    dim = image.shape[1]
    mask = annularmask(dim, dim, innerRadius, outerRadius)
    maskedim = image*mask
    return np.nanmedian(maskedim)

#uses total starlight to find SDI ratio for two images
#returns the best ratio, the remaining starlight, and the SDI image
def SFtotal(HA,cont):
    '''
    finds scale factor for a single pair of images using the ratio of the total starlight in the first
    to the total starlight in the second
    HA: 2D array corresponding to a Hydrogen-Alpha image
    cont: 2D array corresponding to a Continuum image
    outputs the optimal SDI ratio, the total residual starlight when that ratio is used to scale and subtract the 
    continuum image from the HA, and the SDI image as a 2D array
    '''
    ratio = totalStarlight(HA)/totalStarlight(cont)
    image = HA - ratio*cont
    starlight = totalStarlight(image)
    return(ratio, starlight, image)

#uses median starlight to find SDI ratio for two images
#returns the best ratio, the remaining median starlight, and the SDI image
def SFmedian(HA,cont):
    '''
    finds scale factor for a single pair of images using the ratio of the median starlight in the first
    to the median starlight in the second
    HA: 2D array corresponding to a Hydrogen-Alpha image
    cont: 2D array corresponding to a Continuum image
    outputs the optimal SDI ratio, the median residual starlight when that ratio is used to scale and subtract the 
    continuum image from the HA, and the SDI image as a 2D array
    '''
    ratio = medianStarlight(HA)/medianStarlight(cont)
    image = HA - ratio*cont
    starlight = medianStarlight(image)
    return(ratio, starlight, image)


#uses total starlight to find SDI ratio for two cubes of images
#returns the best ratio, the remaining starlight, and the SDI cube
def SFtotalCube(HA, cont):
    '''
    finds the scale factor for each pair of images from two data cubes using the ratios of the total 
    starlight in each image of one cube to the total starlight in each image of the second
    HA: 3D array corresponding to a Hydrogen-Alpha data cube
    cont: 3D array corresponding to a Continuum data cube
    outputs an array of the optimal SDI ratios, the total residual starlight in each resulting SDI image,
    and the SDI cube as a 3D array
    '''
    imNum=HA.shape[0]
    ratios=np.zeros(imNum)
    starlight=np.zeros(imNum)
    cube=np.zeros((imNum,HA.shape[1],HA.shape[1])) 
    for i in range(imNum):
        tempRatio = totalStarlight(HA[i])/totalStarlight(cont[i])
        ratios[i] = tempRatio
        starlight[i] = totalStarlight(HA[i] - tempRatio*cont[i])
        cube[i] = HA[i] - ratios[i]*cont[i]
    return(ratios, starlight, cube)

#uses median starlight to find SDI ratio for two cubes of images
#returns the best ratio, the remaining median starlight, and the SDI cube
def SFmedianCube(HA, cont):
    '''
    finds the scale factor for each pair of images from two data cubes using the ratios of the median 
    starlight in each image of one cube to the median starlight in each image of the second
    HA: 3D array corresponding to a Hydrogen-Alpha data cube
    cont: 3D array corresponding to a Continuum data cube
    outputs an array of the optimal SDI ratios, the median residual starlight in each resulting SDI image,
    and the SDI cube as a 3D array
    '''
    imNum=HA.shape[0]
    ratios=np.zeros(imNum)
    starlight=np.zeros(imNum)
    cube=np.zeros((imNum,HA.shape[1],HA.shape[1])) 
    for i in range(imNum):
        tempRatio = medianStarlight(HA[i])/medianStarlight(cont[i])
        ratios[i] = tempRatio
        starlight[i] = medianStarlight(HA[i] - tempRatio*cont[i])
        cube[i] = HA[i] - ratios[i]*cont[i]
    return(ratios, starlight, cube)

def ratioPlot(medRatio, ratios, dataset, total):
    '''
    plots the scale factor ratio vs. image number (which is a proxy for time) and a horizontal line representing
    the scale factor for the median image
    medRatio: scale factor for the median image
    ratios: array of scale factors for each image
    dataset: name of the dataset (String)
    total: boolean; True if scale factors based on total starlight, False if median starlight
    '''
    x = np.arange(len(ratios))
    medians = np.zeros(len(ratios))
    for i in range(len(ratios)):
        medians[i]=medRatio
    if total:
        plt.plot(x,ratios,color='cornflowerblue',label='total individual ratios')
    else:
        plt.plot(x,ratios,color='cornflowerblue',label='median individual ratios')
    plt.plot(x,medians,color='palevioletred',label='median ratio')
    plt.ylabel('Scale Factor')
    plt.xlabel('Image Number')
    plt.title(str(dataset)+' Ratio Plot')
    plt.legend()

def starlightPlot(starlights, dataset, total):
    '''
    plots the residual starlight after subtraction vs. image number (which is a proxy for time) 
    starlights: array of residual starlight values for each image
    dataset: name of the dataset (String)
    total: boolean; True if scale factors based on total starlight, False if median starlight
    '''
    x = np.arange(len(starlights))
    if total:
        plt.plot(x,starlights,color='salmon',label='total starlight remaining')
    else:
        plt.plot(x,starlights,color='salmon',label='median starlight remaining')
    plt.ylabel('Starlight Remaining (counts)')
    plt.xlabel('Image Number')
    plt.title(str(dataset)+' Starlight Plot')
    plt.legend()

def SDI(subdir, innerRadius, outerRadius):
    '''
    outputs plots of scale factor ratios and residual starlight left over for both the total subtraction and median
    subtraction methods. outputs 6 different SDI fits files:
    "SDI_totSF_med.fits" -- 2D median image after total starlight subtraction
    "SDI_medSF_med.fits" -- 2D median image after median starlight subtraction
    "SDI_totSF_indiv.fits" -- 3D cube that contains the individual SDI subtractions of each image using the 
    total starlight subtraction method
    "SDI_medSF_indiv.fits" -- 3D cube that contains the individual SDI subtractions of each image using the 
    median starlight subtraction method
    "tot_scalefactors.fits" -- 1D output of individual scale factors computed for each image using the total starlight
    subtraction method
    "med_scalefactors.fits" -- 1D output of individual scale factors computed for each image using the median starlight
    subtraction method
    '''
    HA = 'Line_clip451_flat_reg.fits'
    cont = 'Cont_clip451_flat_reg.fits'
    
    
    #assign variables
    dataHA = np.array(fits.getdata(main+subdir+HA))
    dataCONT = np.array(fits.getdata(main+subdir+cont))
    head = fits.getheader(main+subdir+HA)
    
    #calculate medians and totals
    medianHA = np.nanmedian(dataHA,axis=0)
    medianCONT = np.nanmedian(dataCONT,axis=0)
    totRatio, totStarlight, totSDI = SFtotal(medianHA, medianCONT)
    medRatio, medStarlight, medSDI = SFmedian(medianHA, medianCONT)
    
    totRatioCube, totStarlightCube, totSDICube = SFtotalCube(dataHA, dataCONT)
    medRatioCube, medStarlightCube, medSDICube = SFmedianCube(dataHA, dataCONT)
    
    #plot
    plt.figure(figsize=(13,10))
    
    plt.subplot(2,2,1)
    ratioPlot(totRatio, totRatioCube, subdir, True)
    
    plt.subplot(2,2,2)
    starlightPlot(totStarlightCube, subdir, True)
    
    plt.subplot(2,2,3)
    ratioPlot(medRatio, medRatioCube, subdir, False)
    
    plt.subplot(2,2,4)
    starlightPlot(medStarlightCube, subdir, False)
    
    #write out
    #head = "header"
    
    writestr = subdir.replace("/","_")
    fits.writeto(writestr+"SDI_totSF_med.fits", totSDI, head, overwrite=True)
    
    fits.writeto(writestr+"SDI_medSF_med.fits", medSDI, head, overwrite=True)
    
    fits.writeto(writestr+"SDI_totSF_indiv.fits", totSDICube, head, overwrite=True)

    fits.writeto(writestr+"SDI_medSF_indiv.fits", medSDICube, head, overwrite=True)
    
    fits.writeto(writestr+"tot_scalefactors.fits", totRatioCube, head, overwrite=True)
    
    fits.writeto(writestr+"med_scalefactors.fits", medRatioCube, head, overwrite=True)