__author__ = 'JB'


import numpy as np
from copy import copy
from  glob import glob
import csv
import os
import pyklip.kpp.utils.GPIimage as gpiim

def mask_known_objects(cube,prihdr,exthdr,GOI_list_folder = None, mask_radius = 7,
                          include_speckles = False):

    cube_cpy = copy(cube)

    if np.size(cube_cpy.shape) == 3:
        nl,ny,nx = cube_cpy.shape
    elif np.size(cube_cpy.shape) == 2:
        ny,nx = cube_cpy.shape
        cube_cpy = cube_cpy[None,:]
        nl = 1

    width = 2*mask_radius+1
    stamp_x_grid, stamp_y_grid = np.meshgrid(np.arange(0,width,1)-width/2,np.arange(0,width,1)-width/2)
    stamp_mask = np.ones((width,width))
    r_stamp = abs((stamp_x_grid) +(stamp_y_grid)*1j)
    stamp_mask[np.where(r_stamp < mask_radius)] = np.nan

    try:
        # OBJECT: keyword in the primary header with the name of the star.
        object_name = prihdr['OBJECT'].strip().replace (" ", "_")
    except:
        object_name = "UNKNOWN_OBJECT"

    # Get center of the image (star position)
    try:
        # Retrieve the center of the image from the fits headers.
        center = [exthdr['PSFCENTX'], exthdr['PSFCENTY']]
    except:
        # If the keywords could not be found the center is defined as the middle of the image
        print("Couldn't find PSFCENTX and PSFCENTY keywords.")
        center = [(nx-1)/2,(ny-1)/2]

    #Julian Day OBServation
    MJDOBS_fits = prihdr['MJD-OBS']

    row_m = int(np.floor(width/2.0))
    row_p = int(np.ceil(width/2.0))
    col_m = int(np.floor(width/2.0))
    col_p = int(np.ceil(width/2.0))

    if GOI_list_folder is not None:
        object_GOI_filename = GOI_list_folder+os.path.sep+object_name+'_GOI.csv'

        if len(glob(object_GOI_filename)) != 0:
            with open(object_GOI_filename, 'rb') as csvfile_GOI_list:
                GOI_list_reader = csv.reader(csvfile_GOI_list, delimiter=';')
                GOI_csv_as_list = list(GOI_list_reader)
                N_objects = len(GOI_csv_as_list)-1
                attrib_name = GOI_csv_as_list[0]
                GOI_list = np.array(GOI_csv_as_list[1:len(GOI_csv_as_list)])

                pa_id = attrib_name.index("PA")
                sep_id = attrib_name.index("SEP")
                MJDOBS_id = attrib_name.index("MJDOBS")
                STATUS_id = attrib_name.index("STATUS")

                MJDOBS_arr = np.array([ float(it) for it in GOI_list[:,MJDOBS_id]])
                MJDOBS_unique = np.unique(MJDOBS_arr)
                MJDOBS_closest_id = np.argmin(np.abs(MJDOBS_unique-MJDOBS_fits))
                MJDOBS_closest = MJDOBS_unique[MJDOBS_closest_id]
                #Check that the closest MJDOBS is closer than 2 hours
                if abs(MJDOBS_closest-MJDOBS_fits) > 2./24.:
                    # Skip if we couldn't find a matching date.
                    return np.squeeze(cube_cpy)

                for obj_id in np.where(MJDOBS_arr == MJDOBS_closest)[0]:
                    try:
                        pa = float(GOI_list[obj_id,pa_id])
                        radius = gpiim.as2pix(float(GOI_list[obj_id,sep_id]))
                        status = str(GOI_list[obj_id,STATUS_id])
                        if include_speckles or (status in ["Planet","Background","Candidate","Unknown","Brown Dwarf"]):
                            x_max_pos = float(radius)*np.cos(np.radians(90+pa))
                            y_max_pos = float(radius)*np.sin(np.radians(90+pa))
                            col_centroid = x_max_pos+center[0]
                            row_centroid = y_max_pos+center[1]
                            k = int(round(row_centroid))
                            l = int(round(col_centroid))

                            cube_cpy[:,(k-row_m):(k+row_p), (l-col_m):(l+col_p)] = np.tile(stamp_mask,(nl,1,1)) * cube_cpy[:,(k-row_m):(k+row_p), (l-col_m):(l+col_p)]

                    except:
                        print("Missing data in GOI database for {0}".format(object_name))

    for fake_id in range(100):
        try:
            pa = exthdr["FKPA{0:02d}".format(fake_id)]
            radius = exthdr["FKSEP{0:02d}".format(fake_id)]
        except:
            continue

        x_max_pos = float(radius)*np.cos(np.radians(90+pa))
        y_max_pos = float(radius)*np.sin(np.radians(90+pa))
        col_centroid = x_max_pos+center[0]
        row_centroid = y_max_pos+center[1]
        k = int(round(row_centroid))
        l = int(round(col_centroid))

        cube_cpy[:,(k-row_m):(k+row_p), (l-col_m):(l+col_p)] = np.tile(stamp_mask,(nl,1,1)) * cube_cpy[:,(k-row_m):(k+row_p), (l-col_m):(l+col_p)]


    return np.squeeze(cube_cpy)


def get_pos_known_objects(prihdr,exthdr,GOI_list_folder=None,xy = False,pa_sep = False,ignore_fakes = False,fakes_only = False,
                          include_speckles = False,IWA=None,OWA=None):


    try:
        # OBJECT: keyword in the primary header with the name of the star.
        object_name = prihdr['OBJECT'].strip().replace (" ", "_")
    except:
        object_name = "UNKNOWN_OBJECT"

    # Get center of the image (star position)
    try:
        # Retrieve the center of the image from the fits headers.
        center = [exthdr['PSFCENTX'], exthdr['PSFCENTY']]
    except:
        center = [np.nan,np.nan]

    #Julian Day OBServation
    MJDOBS_fits = prihdr['MJD-OBS']



    x_vec = []
    y_vec = []
    col_vec = []
    row_vec = []
    pa_vec = []
    sep_vec = []

    if not fakes_only:
        if GOI_list_folder is not None:
            object_GOI_filename = GOI_list_folder+os.path.sep+object_name+'_GOI.csv'
            if len(glob(object_GOI_filename)) != 0:
                with open(object_GOI_filename, 'rb') as csvfile_GOI_list:
                    GOI_list_reader = csv.reader(csvfile_GOI_list, delimiter=';')
                    GOI_csv_as_list = list(GOI_list_reader)
                    attrib_name = GOI_csv_as_list[0]
                    GOI_list = np.array(GOI_csv_as_list[1:len(GOI_csv_as_list)])

                    pa_id = attrib_name.index("PA")
                    sep_id = attrib_name.index("SEP")
                    MJDOBS_id = attrib_name.index("MJDOBS")
                    STATUS_id = attrib_name.index("STATUS")

                    MJDOBS_arr = np.array([ float(it) for it in GOI_list[:,MJDOBS_id]])
                    MJDOBS_unique = np.unique(MJDOBS_arr)
                    MJDOBS_closest_id = np.argmin(np.abs(MJDOBS_unique-MJDOBS_fits))
                    MJDOBS_closest = MJDOBS_unique[MJDOBS_closest_id]
                    #Check that the closest MJDOBS is closer than 2 hours
                    if abs(MJDOBS_closest-MJDOBS_fits) > 2./24.:
                        # Skip if we couldn't find a matching date.
                        return [],[]

                    for obj_id in np.where(MJDOBS_arr == MJDOBS_closest)[0]:
                        try:
                            pa = float(GOI_list[obj_id,pa_id])
                            radius = float(GOI_list[obj_id,sep_id])
                            if IWA is not None:
                                if gpiim.as2pix(radius) < IWA:
                                    continue
                            if OWA is not None:
                                if gpiim.as2pix(radius) > OWA:
                                    continue
                            status = str(GOI_list[obj_id,STATUS_id])
                            # print(status, include_speckles or (status in ["Planet","Background","Candidate","Unknown","Brown Dwarf"]))
                            if include_speckles or (status in ["Planet","Background","Candidate","Unknown","Brown Dwarf"]):
                                pa_vec.append(pa)
                                sep_vec.append(radius)
                                x_max_pos = float(gpiim.as2pix(radius))*np.cos(np.radians(90+pa))
                                y_max_pos = float(gpiim.as2pix(radius))*np.sin(np.radians(90+pa))
                                x_vec.append(x_max_pos)
                                y_vec.append(y_max_pos)
                                row_vec.append(y_max_pos+center[1])
                                col_vec.append(x_max_pos+center[0])
                        except:
                            print("Missing data in GOI database for {0}".format(object_name))

    if not ignore_fakes:
        for fake_id in range(100):
            try:
                pa = exthdr["FKPA{0:02d}".format(fake_id)]
                radius = gpiim.pix2as(exthdr["FKSEP{0:02d}".format(fake_id)])
                if IWA is not None:
                    if gpiim.as2pix(radius) < IWA:
                        continue
                if OWA is not None:
                    if gpiim.as2pix(radius) > OWA:
                        continue
                pa_vec.append(pa)
                sep_vec.append(radius)
                x_max_pos = float(gpiim.as2pix(radius))*np.cos(np.radians(90+pa))
                y_max_pos = float(gpiim.as2pix(radius))*np.sin(np.radians(90+pa))
                x_vec.append(x_max_pos)
                y_vec.append(y_max_pos)
                row_vec.append(y_max_pos+center[1])
                col_vec.append(x_max_pos+center[0])
            except:
                continue


    if pa_sep:
        return sep_vec,pa_vec
    elif xy:
        return x_vec,y_vec
    else:
        return row_vec,col_vec


def make_GOI_list(outputDir,GOI_list_csv,GPI_TID_csv):
    """
    Generate the GOI files from the GOI table and the TID table (queried from the database).

    :param outputDir: Output directory in which to save the GOI files.
    :param GOI_list_csv: Table with the list of GOIs (including separation, PA...)
    :param GPI_TID_csv: Table giving the TID code for a given object name.
    :return: One .csv file per target for which at list one GOI exists.
            The filename follows: [object]_GOI.csv. For e.g. c_Eri_GOI.csv.
    """
    with open(GOI_list_csv, 'rb') as csvfile_GOI_list:
        GOI_list_reader = csv.reader(csvfile_GOI_list, delimiter=';')
        GOI_csv_as_list = list(GOI_list_reader)
        attrib_name = GOI_csv_as_list[0]
        GOI_list = np.array(GOI_csv_as_list[1:len(GOI_csv_as_list)])

        TID_id = attrib_name.index("TID")
        TID_GOI = GOI_list[:,TID_id]
        TID_unique = np.unique(TID_GOI)

        with open(GPI_TID_csv, 'rb') as csvfile_TID:
            TID_reader = csv.reader(csvfile_TID, delimiter=';')
            TID_csv_as_list = list(TID_reader)
            TID_csv_as_nparr = np.array(TID_csv_as_list)[1:len(TID_csv_as_list),:]
            TID_campaign = np.ndarray.tolist(TID_csv_as_nparr[:,0])
            star_campaign = np.ndarray.tolist(TID_csv_as_nparr[:,1])
            #print(TID_campaign)
            #print(star_campaign)

            dict_matching_TID_to_name = {}
            for TID_it in TID_unique:
                star_raw_name = star_campaign[TID_campaign.index(TID_it)]
                star_name = star_raw_name.replace(" ","_")
                dict_matching_TID_to_name[TID_it] = star_name
            #print(dict_matching_TID_to_name)

        for TID_it in TID_unique:
            where_same_star = np.where(TID_GOI == TID_it)

            with open(outputDir+os.path.sep+dict_matching_TID_to_name[TID_it]+'_GOI.csv', 'w+') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=';')
                table_to_csv = [attrib_name]+np.ndarray.tolist(GOI_list[where_same_star[0],:])#.insert(0,attrib_name)
                #print(table_to_csv)
                csvwriter.writerows(table_to_csv)