anasubs=(res_gran0p005 nores_nogran res_nogran res_gran0p05 res_gran0p005 resx5_nogran resx100_nogran)


for ana in ${anasubs[*]} ; do
    cp anaSubstructure_${ana}.cpp anaSubstructure.cpp
    make
    mv anaSubstructure anaSubstructure_${ana}
    eval `eos root://cmseos.fnal.gov rm /store/user/${USER}/TrackObservablesStudy/anaSubstructure_${ana}`
    xrdcp anaSubstructure_${ana} root://cmseos.fnal.gov:///store/user/${USER}/TrackObservablesStudy/
done
