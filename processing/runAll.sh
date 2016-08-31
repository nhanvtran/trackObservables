anasubs=(nores_nogran res_nogran resx5_nogran resx100_nogran res_gran0p05 res_gran0p005)

for ana in ${anasubs[*]} ; do
    eval `eos root://cmseos.fnal.gov mkdir /store/user/${USER}/${ana}`

    for tfile in $(xrdfs root://cmseos.fnal.gov ls -u /store/user/${USER}/${ana}); do
        filehandle=(${tfile//\// })
        roothandle=${filehandle[${#filehandle[@]}-1]}
        eval `eos root://cmseos.fnal.gov rm /store/user/${USER}/${ana}/${roothandle}`
    done

    rm ./${ana}/*
    python SubmitCondor.py --indir /store/user/${USER}/TrackObservablesStudy/fromMarat/ \
        --outdir ./${ana}/ \
        --maxEvents 50000 \
        --evPerJob 10000 \
        --anaSubLoc root://cmseos.fnal.gov:///store/user/${USER}/TrackObservablesStudy/anaSubstructure_${ana} \
        --outdirname ${ana} 
done
