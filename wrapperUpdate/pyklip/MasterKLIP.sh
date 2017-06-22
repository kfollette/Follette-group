#Elijah Spiro
#Version 1.0 - 9/18/16
#Version 1.1 - 3/7/2017
#Version 2.0 - 3/12/2017



#Increase limit on maximum open files
ulimit -n 4096
echo "[INCREASED MEMORY LIMIT]"

#Launch GUI prompting for information
java GetParameters
echo "[CREATED PARAMETER FILE]"

#Read in parameters from 'parameters.txt'
while IFS='' read -r line || [[ -n "$line" ]]; do
    ((counter+=1))
    if [ $counter = 1 ]; then
	FILEPATH=$line
    elif [ $counter = 2 ]; then
	ANNULI=$line
    elif [ $counter = 3 ]; then
	IWA=$line
    elif [ $counter = 4 ]; then
	MOVEMENT=$line
    elif [ $counter = 5 ]; then
	OUTPUTNAME=$line
    elif [ $counter = 6 ]; then
	KLMODES=$line
    else 
	SUBSECTIONS=$line
    fi
done < "parameters.txt"
echo "[SET INTERNAL PARAMETER VARIABLES]"

#Access IWA information
echo $IWA > "iwa.txt"

#Run StartKLIP.py with command line arguments from these variables
/home/anaconda3/bin/python3 StartKLIP.py $FILEPATH $ANNULI $IWA $MOVEMENT $OUTPUTNAME $KLMODES $SUBSECTIONS

#Ensure the correct file exists
rm $OUTPUTNAME"-KLmodes-all.fits"
cp "_"$OUTPUTNAME"-KLmodes-all.fits" $OUTPUTNAME".fits"
rm "_"$OUTPUTNAME"-KLmodes-all.fits"
mv $OUTPUTNAME".fits" $FILEPATH"/.."
rm iwa.txt
rm parameters.txt

