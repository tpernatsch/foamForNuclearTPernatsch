#!/bin/bash
# Routine to parse the log.GeN-Foam file to extract the reactivity contribution
# at each time step.
# Author: Thomas Guilbaud, EPFL/Transmutex SA, 2022/04/06
# Last update: 2022/04/07

# Test correct user input
if [ $# -ne 1 ]; then
    echo "Usage:  ./parseFoam.sh path/to/log.GeN-Foam"
    exit 1
fi

# Recover user inputs
genfoamFile=$1

folderName=$(echo "$genfoamFile" | rev | cut -d"/" -f2 | rev)

# Extract reactivity values in list
echo "Parse time"
time=$(grep "^Time = " $genfoamFile | cut -d" " -f3)

echo "Parse reactivity"
totalReactivity=$(grep "totalReactivity   = " $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
Doppler=$(grep -P '^(?=.*Doppler)(?=.*pcm)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
TFuel=$(grep -P '^(?=.*TFuel)(?=.*pcm)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2 | tr -d ',')
TClad=$(grep -P '^(?=.*TClad)(?=.*pcm)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
TCool=$(grep -P '^(?=.*TCool)(?=.*pcm)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
rhoCool=$(grep -P '^(?=.*rhoCool)(?=.*pcm)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
TStruct=$(grep -P '^(?=.*TStruct)(?=.*pcm)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
driveline=$(grep -P '^(?=.*driveline)(?=.*pcm)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
Boron=$(grep -P '^(?=.*Boron)(?=.*pcm)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)

# Set output file
if [ -f $genfoamFile.reactivity ]; then
    rm $genfoamFile.reactivity
fi

echo -e "Time[s]\tTotal[pcm]\tDoppler[pcm]\tTFuel[pcm]\tTClad[pcm]\tTCool[pcm]\trhoCool[pcm]\tTStruct[pcm]\tdriveline[pcm]\tBoron[pcm]" >> $genfoamFile.reactivity

# Write in the output file
paste <(
    echo "${time[*]}" ) <(
    echo "${totalReactivity[*]}" ) <(
    echo "${Doppler[*]}" ) <(
    echo "${TFuel[*]}" ) <(
    echo "${TClad[*]}" ) <(
    echo "${TCool[*]}" ) <(
    echo "${rhoCool[*]}" ) <(
    echo "${TStruct[*]}" ) <(
    echo "${driveline[*]}" ) <(
    echo "${Boron[*]}"
) >> $genfoamFile.reactivity

# ------------------------------------------------------------------------------*
echo "Parse external source modulation"
ExtSrcModulation=$(grep 'extSrcModulation' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)

# Set output file
if [ -f $genfoamFile.mod ]; then
    rm $genfoamFile.mod
fi

echo -e "Time[s]\tExtSrcMod " >> $genfoamFile.mod

# Write in the output file
paste <(
    echo "${time[*]}" ) <(
    echo "${ExtSrcModulation[*]}"
) >> $genfoamFile.mod

# ------------------------------------------------------------------------------*
echo "Parse power"
PowerFission=$(grep -P '^(?=.*fission)(?=.*W)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
PowerDecay=$(grep -P '^(?=.*decay)(?=.*W)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
PowerExtSrcBeam=$(grep -P '^(?=.*extSrcBeamPower)(?=.*W)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)

# Set output file
if [ -f $genfoamFile.power ]; then
    rm $genfoamFile.power
fi

echo -e "Time[s]\tPwrFiss[W]\tPwrDecay[W]\tPwrSrcBeam[W]" >> $genfoamFile.power

# Write in the output file
paste <(
    echo "${time[*]}" ) <(
    echo "${PowerFission[*]}" ) <(
    echo "${PowerDecay[*]}" ) <(
    echo "${PowerExtSrcBeam[*]}"
) >> $genfoamFile.power

# ------------------------------------------------------------------------------*
echo "Parse temperature"
TFuelTemp=$(grep -P '^(?=.*TFuel)(?=.*K)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
TCladTemp=$(grep -P '^(?=.*TClad)(?=.*K)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)
TCoolTemp=$(grep -P '^(?=.*TCool)(?=.*K)' $genfoamFile | cut -d"=" -f2 | cut -d" " -f2)

# Set output file
if [ -f $genfoamFile.temperature ]; then
    rm $genfoamFile.temperature
fi

echo -e "Time[s]\tTFuel[K]\tTClad[K]\tTCool[K]" >> $genfoamFile.temperature

# Write in the output file
paste <(
    echo "${time[*]}" ) <(
    echo "${TFuelTemp[*]}" ) <(
    echo "${TCladTemp[*]}" ) <(
    echo "${TCoolTemp[*]}"
) >> $genfoamFile.temperature

# ------------------------------------------------------------------------------*

echo "Extraction completed"

# ------------------------------------------------------------------------------*
