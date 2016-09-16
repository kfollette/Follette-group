#Written by Elijah Spiro
#Last edited 8/18/2016

#Syntax to run: [./Program name] [directory location relative from code location] ["cube" for long version, empty for short version] [H-Alpha data cube] [Continuum data cube]
#Example: ./GenerateScale.sh "../HD100546/12Apr14/" "cube" "Line_clip450_reg_circsym.fits" "Cont_clip450_reg_circsym.fits"
# --> Argument 2 allows you to specify whether to run in "cube" mode. This will take ~5x longer, but will produce a unique ratio for every image. Otherwise, you will get a quick and dirty (median) answer

#######################################################################################                                                                                                                    
#                                      ABOUT                                          #                                                                                                                    
#######################################################################################
#This script is designed to produce the perfect scaling ratio for an arbitrary number of files in a FITS image cube. 
#Requirements: Files "Line_clip450_reg_circsym.fits" and "Cont_clip450_reg_circsym.fits" must be in the directory specified by command line argument (format to navigate to this directory from coding directory). Additionally, python files "GenerateMedianScale.py" and "GenerateCubeScale.py" must be in the same directory as this bash script.
#Output: This script will produce 8 FITS files ("SDI_MEDIAN_SR.fits", "SDI_CUBE_SR.fits", "SDI_CUBE_MR.fits", "scale_factors.fits", "MEDIAN_MASK.fits", "CUBE_MASK.fits", "MASKED_DATA_HA.fits", "MASKED_DATA_CONT.fits") and create 2 plots near completion. 
#Expect long runtime (~20 seconds per image in data cube)
#If the script must be halted prior to completion, and any of the eight output files have already been created, please delete them before running the script again.


#######################################################################################                                                                                                                    
#                               ORDER OF EXECUTION                                    #                                                                                                                    
####################################################################################### 
#Due to the long nature of this script, progress will be displayed at every interval. Order of program execution:
#1. Calculate median of H-Alpha and Continuous spectrum data cubes (short)
#2. Scan median image to determine saturation radius, produce "MEDIAN_MASK.fits" (very short)
#3. Apply mask to both median images (very short)
#4. Scan image cube to determine saturation radius, produce "CUBE_MASK.fits" (medium)
#5. Apply mask to each file in both image cubes, produce "MASKED_DATA_HA.fits" and "MASKED_DATA_CONT.fits" for later use (long)
#6. Recursively generate scaling ratio for median image that will yield lowest remaining starlight counts [repeats with alternate method, averages results] (short)
#7. Multiply median images from step 1 by ratio, subtract corresponding pixels, produce "SDI_MEDIAN_SR.fits" (short)
#8. Multiply each image in Continuous spectrum data cube by ratio, subtract corresponding pixels, produce "SDI_CUBE_SR.fits" (long)
#9. Read in masked image cubes created in step 5
#10. Recursively generate scaling ratio for every image in Continuous spectrum data cube (very long)
#11. Generate 1D Fits file with list of these ratios for later use --> "scale_factors.fits"
#12. Multiply each image in continuous spectrum data cube by unique ratio, subtract corresponding pixels, produce "SDI_CUBE_MR.fits" (long)
#13. Produce plot of unique ratios vs. median ratio (short)
#14. Produce plot of starlight remaining post-subtraction in SDI_CUBE_SR vs. SDI_CUBE_MR (medium)


#######################################################################################
#                               BEGIN PROGRAM FLOW                                    #
#######################################################################################
clear
echo "---------------------------------------------------------------------------------"
echo "|          Now running total scaling script. This might take a while.           |"           
echo "---------------------------------------------------------------------------------"
echo ""

##shell scripts do not by default expand aliases, so this is to make it work with gpicruncher. shouldn't affect how it is executed on your machine
shopt -s expand_aliases
source ~/.bashrc

python GenerateMedianScale.py $1 $3 $4

if [ $2 == "cube" ]; then 
    python GenerateCubeScale.py $1 $3 $4
fi

for i in `seq 1 10`;
do
    echo ""
done

echo "---------------------------------------------------------------------------------"
echo "|            Program ran successfully. See FITS files for results.              |"
echo "---------------------------------------------------------------------------------"
