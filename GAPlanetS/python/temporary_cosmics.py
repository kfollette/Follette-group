from astropy.io import fits
import numpy as np

def remove_cosmics(line,cont,cosmics,rotoff):
    hdulist = fits.open(cont)
    Cont_cube = hdulist[0].data
    header_data = fits.getheader(cont)
    hdulist.close()
    hdulist = fits.open(line)
    Line_cube = hdulist[0].data
    hdulist.close()
    hdulist = fits.open(cosmics)
    line_cosmics = hdulist[0].data
    hdulist.close()
    hdulist = fits.open(cosmics)
    cont_cosmics = hdulist[0].data
    hdulist.close()
    hdulist = fits.open(rotoff)
    rotoff_line = hdulist[0].data
    rotoff_cont = rotoff_line
    hdulist.close()
    
    marker = 0
    count = 0
    first = 0
    for a in range(Cont_cube.shape[0]):
        if not a in cont_cosmics:
            rotoff_cont = np.delete(rotoff_cont,a-count)
            count = count+1
        if a in cont_cosmics:
            print("beginning cube with layer " + str(a))
            new_cont = Cont_cube[a]
            marker = a
            break;
    for a in range(marker+1,Cont_cube.shape[0]):
        if not a in cont_cosmics:
            rotoff_cont = np.delete(rotoff_cont,a-count)
            count = count+1
        if a in cont_cosmics:
            print("Adding layer " +str(a))
            if(new_cont.shape[0]==451 and first==0):
                new_cont = np.vstack(([new_cont],[Cont_cube[a]]))
                first = 1
            else:
                new_cont = np.vstack((new_cont,[Cont_cube[a]]))
                
    marker = 0
    count = 0
    first = 0
    for a in range(Line_cube.shape[0]):
        if not a in line_cosmics:
            rotoff_line = np.delete(rotoff_line,a-count)
            count = count+1
        if a in line_cosmics:
            print("beginning cube with layer " + str(a))
            new_line = Line_cube[a]
            marker = a
            break;
    for a in range(marker+1,Line_cube.shape[0]):
        if not a in line_cosmics:
            rotoff_line = np.delete(rotoff_line,a-count)
            count = count+1
        if a in line_cosmics:
            print("Adding layer " +str(a))
            if(new_line.shape[0]==451 and first==0):
                new_line = np.vstack(([new_line],[Line_cube[a]]))
                first = 1
            else:
                new_line = np.vstack((new_line,[Line_cube[a]]))

    lineName = line[0:-5] + "_nocosmics.fits"
    contName = cont[0:-5] + "_nocosmics.fits"
            
    print("Writing to fits files.")    
    fits.writeto(contName, new_cont, header=header_data, overwrite = True)
    fits.writeto(lineName, new_line, header=header_data, overwrite = True)
    fits.writeto('rotoff_nocosmics.fits', rotoff_line, overwrite = True)
    #fits.writeto('rotoff_nocosmics_cont.fits', rotoff_cont, overwrite = True)
    #fits.writeto('rotoff_nocosmics_line.fits', rotoff_line, overwrite = True)
