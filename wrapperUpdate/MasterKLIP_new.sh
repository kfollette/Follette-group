#Elijah Spiro
#Version 2.1 - 3/26/17
#edited by Clare Leonard - 8/1/2017
#Version 3.1 - 7/12/2018

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
    python RunKLIP.py `< single_reduction_parameters.txt`
    rm single_reduction_parameters.txt
fi
if [ -f automation_parameters.txt ]; then
    echo "Performing automation"
    echo ""
    python AutomateKLIP.py `< automation_parameters.txt`
    rm automation_parameters.txt
fi


if [ -f "single_reduction_parameters.txt~" ]; then
   rm single_reduction_parameters.txt~
fi

if [ -f "automation_parameters.txt~" ]; then
   rm automation_parameters.txt~
fi


