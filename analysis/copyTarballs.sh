#!/bin/sh

DIR='/store/user/ntran/SUSY/theory_JPM/training'
trainingDir='training'

for file in /eos/uscms$DIR/op_*.tar.gz
do
    filebase=${file##*/}
    echo $file, $filebase
    xrdcp root://cmseos.fnal.gov/$DIR/$filebase $trainingDir/$filebase
done
