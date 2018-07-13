#Elijah Spiro
#Version 2.1 - 3/26/17
#edited by Clare Leonard - 8/1/2017

#Increase limit on maximum open files
ulimit -n 4096

#Launch GUI prompting for information
if  [ ! -f GetParameters.class ]; then
  javac GetParameters.java
fi

#javac GetParameters.java
java GetParameters
rm *.class

#Determine whether to run single reduction of automated reduction
if [ -f single_reduction_parameters.txt ]; then
    echo "Performing single reduction"
    echo ""
    RUNTYPE=1
fi
if [ -f automation_parameters.txt ]; then
    echo "Performing automation"
    echo ""
    RUNTYPE=2
fi

#Read in parameters from one of the two parameter files
if [[ $RUNTYPE = 1 ]]; then
    while IFS='' read -r line || [[ -n "$line" ]]; do
	((counter+=1))
	if [ $counter = 1 ]; then
	    FILEPATH=$line
    elif [ $counter = 2 ]; then
	    IWA=$line
    elif [ $counter = 3 ]; then
	    KLMODES=$line
	elif [ $counter = 4 ]; then
	    ANNULI=$line
	elif [ $counter = 5 ]; then
	    MOVEMENT=$line
	elif [ $counter = 6 ]; then
	    SUBSECTIONS=$line
    elif [ $counter = 7 ]; then
	    OUTPUTNAME=$line
	elif [ $counter = 8 ]; then
        SNR=$line
	elif [ $counter = 9 ]; then
            SAVE1=$line
	elif [ $counter = 10 ]; then
            FWHM=$line
    elif [ $counter = 11 ]; then
            SMOOTH=$line
	elif [ $counter = 12 ]; then
            RA2=$line
	elif [ $counter = 13 ]; then
            PA2=$line
	elif [ $counter = 14 ]; then
            WID2=$line
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
	elif [ $counter = 13 ]; then
            S3=$line
	elif [ $counter = 14 ]; then
            FWHM=$line
    elif [ $counter = 15 ]; then
            SMOOTH=$line
	elif [ $counter = 16 ]; then
            ra=$line
	elif [ $counter = 17 ]; then
            pa=$line
	elif [ $counter = 18 ]; then
            wi=$line
	elif [ $counter = 19 ]; then
            SAVE2=$line
        fi
    done < "automation_parameters.txt"
fi


#Run StartKLIP.py with command line arguments from these variables (if single reduction)
if [[ $RUNTYPE = 1 ]]; then
    python RunKLIP.py $FILEPATH $IWA $KLMODES $ANNULI $MOVEMENT $SUBSECTIONS $OUTPUTNAME $SNR $SAVE1 $FWHM $SMOOTH $RA2 $PA2 $WID2
    #/home/anaconda3/bin/python3 StartKLIP.py $FILEPATH $ANNULI $IWA $MOVEMENT $OUTPUTNAME $KLMODES $SUBSECTIONS $SNR
    #rm $OUTPUTNAME"-KLmodes-all.fits"
    #cp "_"$OUTPUTNAME"-KLmodes-all.fits" $OUTPUTNAME".fits"
    #rm "_"$OUTPUTNAME"-KLmodes-all.fits"
    #cp $OUTPUTNAME".fits" "temp_klip.fits"
    #rm temp_klip.fits
else 
    python AutomateKLIP.py $FILEPATH $IWA $KLMODES $OUTPUTNAME $A1 $A2 $A3 $M1 $M2 $M3 $S1 $S2 $S3 $FWHM $SMOOTH $ra $pa $wi $SAVE2
fi


     
if [ -f "single_reduction_parameters.txt" ]; then
   rm single_reduction_parameters.txt
fi

if [ -f "automation_parameters.txt" ]; then
   rm automation_parameters.txt
fi

if [ -f "single_reduction_parameters.txt~" ]; then
   rm single_reduction_parameters.txt~
fi

if [ -f "automation_parameters.txt~" ]; then
   rm automation_parameters.txt~
fi

