#Elijah Spiro
#Version 2.1 - 3/26/17

#Increase limit on maximum open files
ulimit -n 4096

#Launch GUI prompting for information
javac GetParameters.java
java GetParameters
rm GetParameters.class

#Determine whether to run single reduction of automated reduction
if [ -f single_reduction_parameters.txt ]; then
    echo "Performing single reduction"
    RUNTYPE=1
fi
if [ -f automation_parameters.txt ]; then
    echo "Performing automation"
    RUNTYPE=2
fi

#Read in parameters from one of the two parameter files
if [ $RUNTYPE = 1 ]; then
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
	elif [ $counter = 7 ]; then
	    SUBSECTIONS=$line
	else
	    SNR=$line
	fi
    done < "single_reduction_parameters.txt"
else
    while IFS='' read -r line || [[ -n "$line" ]]; do
        ((counter+=1))
        if [ $counter = 1 ]; then
            FILEPATH=$line
        elif [ $counter = 2 ]; then
            IWA=$line
        elif [ $counter = 3 ]; then
            KLMODES=$line
        elif [ $counter = 4 ]; then
            OUTPUTNAME=$line
        elif [ $counter = 5 ]; then
            A1=$line
        elif [ $counter = 6 ]; then
            A2=$line
        elif [ $counter = 7 ]; then
            A3=$line
        elif [ $counter = 8 ]; then
            M1=$line
        elif [ $counter = 9 ]; then
            M2=$line
        elif [ $counter = 10 ]; then
            M3=$line
	elif [ $counter = 11 ]; then
            S1=$line
        elif [ $counter = 12 ]; then
            S2=$line
	else
            S3=$line
        fi
    done < "automation_parameters.txt"
fi

#Access IWA information
echo $IWA > "iwa.txt"

#Run StartKLIP.py with command line arguments from these variables (if single reduction)
if [ $RUNTYPE = 1 ]; then
    python StartKLIP.py $FILEPATH $ANNULI $IWA $MOVEMENT $OUTPUTNAME $KLMODES $SUBSECTIONS $SNR
    #/home/anaconda3/bin/python3 StartKLIP.py $FILEPATH $ANNULI $IWA $MOVEMENT $OUTPUTNAME $KLMODES $SUBSECTIONS $SNR
    rm $OUTPUTNAME"-KLmodes-all.fits"
    cp "_"$OUTPUTNAME"-KLmodes-all.fits" $OUTPUTNAME".fits"
    rm "_"$OUTPUTNAME"-KLmodes-all.fits"
    cp $OUTPUTNAME".fits" "temp_klip.fits"
    if [ $SNR ]; then
	python SNRMap.py
	#/home/anaconda3/bin/python3 SNRMap
	python AppendFiles.py $OUTPUTNAME
	#/home/anaconda3/bin/python3 AppendFiles
	rm custom_mask.fits
	rm mask_cube.fits
	rm multiplied_cube.fits
	rm nose_map.fits
	rm radial_profile.fits
	rm snr_map.fits
    fi
    rm temp_klip.fits
else 
    python AutomateKLIP.py $FILEPATH $IWA $KLMODES $OUTPUTNAME $A1 $A2 $A3 $M1 $M2 $M3 $S1 $S2 $S3
    #/home/anaconda3/bin/python3 AutomateKLIP.py $FILEPATH $IWA $KLMODES $OUTPUTNAME $A1 $A2 $A3 $M1 $M2 $M3 $S1 $S2 $S3
fi

#Ensure the correct file exists
mv $OUTPUTNAME".fits" $FILEPATH"/.."
rm iwa.txt     
if [ -f "single_reduction_parameters.txt" ]; then
   rm single_reduction_parameters.txt
   rm single_reduction_parameters.txt~
fi
if [ -f "automation_parameters" ]; then
   rm automation_parameters.txt
   rm automation_parametetrs.txt~
fi

rm GetParameters.class