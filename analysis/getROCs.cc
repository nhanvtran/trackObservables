/****
 * An example macro showcasing the plotROC.C code
 */

#include "./plotROC.C"

void getROCs(TString inputDir,TString outDir="./") {
    
    // analysis information 
    std::vector<TString> type = {"t", "W", "Z", "q", "g" };
    std::vector<TString> trees = { "allpar", "tracks", "tragam" };
    std::vector<TString> treeNames = { "all par.", "tracks", "tra+#gamma" };
    std::vector<TString> pts = { "pt1", "pt5" };
    std::vector<TString> ptNames = { "p_{T} 1 TeV", "p_{T} 5 TeV" };
    std::vector<TString> trainings = { "shapesonly", "massonly", "all" };
    std::vector<TString> trainingsNames = { "shapes", "mass", "all var." };
    outputDir=outDir+"/";

    
    // setup vectors of settings
    std::vector<int> colorSettings;
    std::vector<std::vector<TString>> settingsMap;


    // mass-shapes-all for fixed tree, sig/bkg, pt
    for(size_t i=0; i<trees.size(); i++) {
    for(size_t j=0; j<type.size();  j++) {
    for(size_t k=0; k<type.size();  k++) {
        if(j==k) continue;
    for(size_t l=0; l<pts.size();   l++) {
        TString tree=trees.at(i);
        TString sig = type.at(j);
        TString bkg = type.at(k);
        TString pt  = pts.at(l);

        std::vector<TString> settings={};
        
        for(size_t itrain=0; itrain < trainings.size(); itrain++) {
           TString training = trainings.at(itrain);

           settings.push_back(inputDir+"/"+training+"/MVA_bdtg_"
                                      +sig+sig+"-"+pt+"_"+bkg+bkg
                                      +"-"+pt+"_t_"+tree+".root");
           settings.push_back(sig+" vs "+bkg+" "+treeNames.at(i)+
                                            ", "+trainingsNames.at(itrain));
        }

        settings.push_back(sig+"_v_"+bkg+"_cmptrainings_"+pt+"_"+tree);

        colorSettings.push_back(trainings.size());
        settingsMap.push_back(settings);
    }}}}


    // all-tree for fixed training, sig/bkg, pt
    for(size_t i=0; i<trainings.size(); i++) {
    for(size_t j=0; j<type.size();  j++) {
    for(size_t k=0; k<type.size();  k++) {
        if(j==k) continue;
    for(size_t l=0; l<pts.size();   l++) {
        TString training=trainings.at(i);
        TString sig = type.at(j);
        TString bkg = type.at(k);
        TString pt  = pts.at(l);

        std::vector<TString> settings={};
        
        for(size_t itree=0; itree < trees.size(); itree++) {
           TString tree = trees.at(itree);

           settings.push_back(inputDir+"/"+training+"/MVA_bdtg_"
                                      +sig+sig+"-"+pt+"_"+bkg+bkg
                                      +"-"+pt+"_t_"+tree+".root");
           settings.push_back(sig+" vs "+bkg+" "+tree+", "+training);
        }

        settings.push_back(sig+"_v_"+bkg+"_cmptrees_"+pt+"_"+trainingsNames.at(i));

        colorSettings.push_back(trees.size());
        settingsMap.push_back(settings);
    }}}}


    // all-tree && pt1,pt5 for fixed training, sig/bkg
    for(size_t i=0; i<trainings.size(); i++) {
    for(size_t j=0; j<type.size();  j++) {
    for(size_t k=0; k<type.size();  k++) {
        if(j==k) continue;
        TString training=trainings.at(i);
        TString sig = type.at(j);
        TString bkg = type.at(k);

        std::vector<TString> settings={};
        
        for(size_t l=0; l<pts.size();   l++) {
            TString pt  = pts.at(l);
            for(size_t itree=0; itree < trees.size(); itree++) {
               TString tree = trees.at(itree);

               settings.push_back(inputDir+"/"+training+"/MVA_bdtg_"
                                          +sig+sig+"-"+pt+"_"+bkg+bkg
                                          +"-"+pt+"_t_"+tree+".root");
               settings.push_back(sig+" vs "+bkg+" "+treeNames.at(itree)
                                            +", "+trainingsNames.at(i));
            }
        }

        settings.push_back(sig+"_v_"+bkg+"_cmptrees_ptall_"+training);

        colorSettings.push_back(trees.size());
        settingsMap.push_back(settings);
    }}}
    

    // all-procs && pt1,pt5 for fixed training, sig/bkg
    for(size_t i=0; i<trainings.size(); i++) {
    for(size_t j=0; j<type.size();  j++) {
    for(size_t itree=0; itree < trees.size(); itree++) {
        TString tree = trees.at(itree);
        TString training=trainings.at(i);
        TString sig = type.at(j);

        std::vector<TString> settings={};
        
        for(size_t l=0; l<pts.size();   l++) {
            TString pt  = pts.at(l);
            for(size_t k=0; k<type.size();  k++) {
               if(j==k) continue;
               TString bkg = type.at(k);

               settings.push_back(inputDir+"/"+training+"/MVA_bdtg_"
                                          +sig+sig+"-"+pt+"_"+bkg+bkg
                                          +"-"+pt+"_t_"+tree+".root");
               settings.push_back(sig+" vs "+bkg+" "+treeNames.at(itree)
                                                    +", "+trainingsNames.at(i));
            }
        }

        settings.push_back(sig+"_v_all_"+tree+"_"+training+"_ptall");

        colorSettings.push_back(type.size()-1);
        settingsMap.push_back(settings);
    }}}


    // all-procs for fixed training, sig/bkg, pt
    for(size_t i=0; i<trainings.size(); i++) {
    for(size_t j=0; j<type.size();  j++) {
    for(size_t itree=0; itree < trees.size(); itree++) {
    for(size_t l=0; l<pts.size();   l++) {
        TString tree = trees.at(itree);
        TString training=trainings.at(i);
        TString sig = type.at(j);
        TString pt  = pts.at(l);

        std::vector<TString> settings={};
        
        for(size_t k=0; k<type.size();  k++) {
            if(j==k) continue;
            TString bkg = type.at(k);

            settings.push_back(inputDir+"/"+training+"/MVA_bdtg_"
                    +sig+sig+"-"+pt+"_"+bkg+bkg
                    +"-"+pt+"_t_"+tree+".root");
            settings.push_back(sig+" vs "+bkg+" "+treeNames.at(itree)+", "+trainingsNames.at(i));
        }

        settings.push_back(sig+"_v_all_"+tree+"_"+training+"_"+pt);

        colorSettings.push_back(type.size()-1);
        settingsMap.push_back(settings);
    }}}}


    // Generate plots for all input settings
    for(size_t i=0; i<settingsMap.size(); i++) {
        //for(size_t j=0; j < settingsMap.at(i).size(); j++) {
        //    std::cout<<settingsMap.at(i).at(j)<<std::endl;
        //}
        iCMax = colorSettings.at(i); 
        //std::cout<<iCMax<<std::endl;
        plotROC(settingsMap.at(i));
        gROOT->CloseFiles();
    }
}
