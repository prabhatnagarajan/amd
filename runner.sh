#!/bin/bash
TIMESTAMP="$(date +"%T")"
#File where data should be outputted by Program
DATAFILE="data$TIMESTAMP.txt"
#file where times should be outputted
TIMEFILE="times$TIMESTAMP.txt"
#rm $DATAFILE
#rm $TIMEFILE

javac -Xlint -classpath ".:commons-math3-3.6.1/commons-math3-3.6.1.jar:.:/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/cplex.jar" *.java

CORRELATION=0.5
#number of samples from "ground truth" distribution
NUMACTUALSAMPLES=50000
NUMDIRICHLETSAMPLES=10000
CONFIDENCE=0.95

echo "Num Agent 1 Types, Num Agent 2 Types, Number of Actual Samples, Number of Dirichlet Samples (per cond dist), Correlation, Confidence Level, Number of Violated IR Constraints, Number of Violated IC Constraints, Full Surplus, Ex-Post Surplus, Robust Surplus, Learn Time (milliseconds), Ex-Post Mechanism Time (millseconds), Total Time" >> $TIMEFILE

#outer loop, loop through variable values
for z in `seq 1 1`
do
#loop for as many iterations as you would like
for w in `seq 1 1`
do

#echo "NUMACTUALSAMPLES = $NUMACTUALSAMPLES, NUMDIRICHETSAMPLES = $NUMDIRICHLETSAMPLES, CONFIDENCE = $CONFIDENCE, CORRELATION = $CORRELATION" >> $TIMEFILE

#number of types for each bidder
for i in `seq 4 4`;
do

for k in `seq 4 $i`;
do

#parameters to R script plot.r (in order): number bidders, number of types of bidder 1, number of types of bidder 2, correlation, FILENAME
#FILENAME should be set to the name of the file that the discrete dist will be written to                                                                   
#file name should have no spaces
FILENAME="mymatrix.txt"
rm $FILENAME

Rscript plot.r 2 $i $k $CORRELATION $FILENAME

#file name should have no spaces
TYPESFILE="types.txt"
rm $TYPESFILE
#The bidders' types should be separated by spaces
#NOTE: The number of types specified here should be consistent with the number of types specified for R script
SPACE=" "
BIDDERONETYPES=""
BIDDERTWOTYPES=""
for x in `seq 1 $i`;
do
	BIDDERONETYPES=$BIDDERONETYPES$x$SPACE
done
for y in `seq 1 $k`;
do
	BIDDERTWOTYPES=$BIDDERTWOTYPES$y$SPACE
done

echo "$BIDDERONETYPES" >> "$TYPESFILE"
echo "$BIDDERTWOTYPES" >> "$TYPESFILE"
#echo "NumTypesBidderOne = $i, NumTypesBidderTwo = $k" >> $TIMEFILE

#java -classpath ".:commons-math3-3.6.1/commons-math3-3.6.1.jar:.:/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/cplex.jar" -Djava.library.path=/opt/ibm/ILOG/CPLEX_Studio1263/cplex/bin/x86-64_linux -Xmx4g AutomatedMechanismDesign $TYPESFILE $FILENAME $NUMACTUALSAMPLES $NUMDIRICHLETSAMPLES $CONFIDENCE $CORRELATION 

TIME=`(time (java -classpath ".:commons-math3-3.6.1/commons-math3-3.6.1.jar:.:/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/cplex.jar" -Djava.library.path=/opt/ibm/ILOG/CPLEX_Studio1263/cplex/bin/x86-64_linux -Xmx4g AutomatedMechanismDesign $TYPESFILE $FILENAME $NUMACTUALSAMPLES $NUMDIRICHLETSAMPLES $CONFIDENCE $CORRELATION >> $TIMEFILE)) 2>&1 | grep real | awk '{print $2}'`

echo $TIME >> $TIMEFILE

#WARNING: Redirect to TIMEFILE time is accounted for in bash time command (I think)
#time (java -classpath ".:commons-math3-3.6.1/commons-math3-3.6.1.jar:.:/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/cplex.jar" -Djava.library.path=/opt/ibm/ILOG/CPLEX_Studio1263/cplex/bin/x86-64_linux -Xmx4g AutomatedMechanismDesign $TYPESFILE $FILENAME $NUMACTUALSAMPLES $NUMDIRICHLETSAMPLES $CONFIDENCE >> $DATAFILE) 2>> $TIMEFILE

done

done

done
#CONFIDENCE=$(bc <<< "scale=3; $CONFIDENCE + 0.01")
#let NUMACTUALSAMPLES=$NUMACTUALSAMPLES*2
#CORRELATION=$(bc <<< "scale=4; $CORRELATION + 0.1")
done
