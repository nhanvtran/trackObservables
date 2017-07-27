/****
 * An example macro showcasing the plotROC.C code
 */

#include "./plotROC.C"

void getROCs(TString inputDir,TString outDir="./",TString anasub="r0_h0_e0") {
    
    // analysis information 
    std::vector<TString> type = {"t", "W", "Z", "q", "g" };
    std::vector<TString> trees = { "allpar", "tracks", "tragam" };
    std::vector<TString> treeNames = { "All particles", "Tracks only", "Tracks+#gamma" };
    std::vector<TString> pts = { "pt1", "pt5" };
    std::vector<TString> ptNames = { "p_{T} 1 TeV", "p_{T} 5 TeV" };
    std::vector<TString> trainings = { "shapesonly", "massonly", "all" };
    std::vector<TString> trainingsNames = { "Shapes only", "Mass only", "All variables" };
    std::map<TString,TString> smearNames = { 
        {"r0_h0_e0",        "Perfect" }, 
        {"r05_h05_e005",    "HCAL0.05"}, 
        {"r05_h01_e005",    "HCAL0.01"}, 
        {"r05_h01_e005_t500",  "F1"}, 
        {"r05_h002_e001_t500", "F2"},
        {"r05_h01_e005_t220",  "F1"}, 
        {"r05_h002_e001_t220", "F2"},
        {"r05_h01_e005_t",  "F1"}, 
        {"r05_h002_e001_t", "F2"},
        {"r05_h002_e005_t500", "Extreme-gran."} ,
        {"r05_h005_e005_t500", "High-gran."},
        { "r1_h022_e050_t110", "CMS-like (EB)"},
        { "r1_h022_e0175_t220","CMS-like"}
        };
    TString smearName=smearNames.at(anasub);

    std::cout<<smearName.Data()<<std::endl;

    outputDir=outDir+"/";

    
    // setup vectors of settings
    std::vector<TString> extraTexts;
    std::vector<int> colorSettings;
    std::vector<std::vector<TString>> settingsMap;

    // mass-shapes-all x pt for fixed tree, sig/bkg 
    for(size_t i=0; i<trees.size(); i++) {
    for(size_t j=0; j<type.size();  j++) {
    for(size_t k=0; k<type.size();  k++) {
        if(j==k) continue;
        TString tree=trees.at(i);
        TString sig = type.at(j);
        TString bkg = type.at(k);

        std::vector<TString> settings={};
        
        for(size_t l=0; l<pts.size();   l++) {
        for(size_t itrain=0; itrain < trainings.size(); itrain++) {
           TString training = trainings.at(itrain);
           TString pt  = pts.at(l);

           settings.push_back(inputDir+"/"+training+"/MVA_bdtg_"
                                      +sig+sig+"-"+pt+"_"+bkg+bkg
                                      +"-"+pt+"_t_"+tree+".root");
           settings.push_back(trainingsNames.at(itrain)
                              +(itrain==0 ? ", "+ptNames.at(l) : ""));
        }}

        settings.push_back(sig+"_v_"+bkg+"_cmptrainings_ptall_"+tree);

        TString extra = "#scale[0.8]{#splitline{"+treeNames.at(i)+"}{"+smearName+"}}";
        extra = "#splitline{"+sig+" vs. "+bkg+"}{"+extra+"}";
        extraTexts.push_back(extra);

        colorSettings.push_back(trainings.size());
        settingsMap.push_back(settings);
    }}}


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
           settings.push_back(trainingsNames.at(itrain));
        }

        settings.push_back(sig+"_v_"+bkg+"_cmptrainings_"+pt+"_"+tree);

        TString extra = "#scale[0.8]{#splitline{"+treeNames.at(i)+"}{"+ptNames.at(l)+"}}";
        extra = "#splitline{"+sig+" vs. "+bkg+"}{"+extra+"}";
        extra = "#splitline{"+smearName+"}{"+extra+"}";
        extraTexts.push_back(extra);

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
           settings.push_back(treeNames.at(itree));
        }

        settings.push_back(sig+"_v_"+bkg+"_cmptrees_"+pt+"_"+training);

        TString extra = "#scale[0.8]{#splitline{"+trainingsNames.at(i)+"}{"+ptNames.at(l)+"}}";
        extra = "#splitline{"+sig+" vs. "+bkg+"}{"+extra+"}";
        extra = "#splitline{"+smearName+"}{"+extra+"}";
        extraTexts.push_back(extra);

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
               settings.push_back(treeNames.at(itree)+", "+ptNames.at(l));
            }
        }

        settings.push_back(sig+"_v_"+bkg+"_cmptrees_ptall_"+training);

        TString extra = "#splitline{"+sig+" vs. "+bkg+"}{#scale[0.8]{"+trainingsNames.at(i)+"}}";
        extra = "#splitline{"+smearName+"}{"+extra+"}";
        extraTexts.push_back(extra);

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
               settings.push_back(sig+" vs "+bkg+", "+ptNames.at(l));
            }
        }

        settings.push_back(sig+"_v_all_"+tree+"_"+training+"_ptall");

        TString extra = "#splitline{"+treeNames.at(itree)+"}{#scale[0.8]{"+trainingsNames.at(i)+"}}";
        extra = "#splitline{"+smearName+"}{"+extra+"}";
        extraTexts.push_back(extra);

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
            settings.push_back(sig+" vs "+bkg);
        }

        settings.push_back(sig+"_v_all_"+tree+"_"+training+"_"+pt);

        TString extra = "#scale[0.8]{#splitline{"+trainingsNames.at(i)+"}{"+ptNames.at(l)+"}}";
        extra = "#splitline{"+treeNames.at(itree)+"}{"+extra+"}";
        extra = "#splitline{"+smearName+"}{"+extra+"}";
        extraTexts.push_back(extra);

        colorSettings.push_back(type.size()-1);
        settingsMap.push_back(settings);
    }}}}


    // all trainings for fixed tree, sig/bkg, pt
    for(size_t j=0; j<type.size();  j++) {
    for(size_t itree=0; itree < trees.size(); itree++) {
    for(size_t l=0; l<pts.size();   l++) {
    for(size_t k=0; k<type.size();  k++) {
        if(j==k) continue;
        TString tree = trees.at(itree);
        TString sig = type.at(j);
        TString bkg = type.at(k);
        TString pt  = pts.at(l);

        std::vector<TString> settings={};

        for(size_t i=0; i<trainings.size(); i++) {
            TString training=trainings.at(i);

            settings.push_back(inputDir+"/"+training+"/MVA_bdtg_"
                    +sig+sig+"-"+pt+"_"+bkg+bkg
                    +"-"+pt+"_t_"+tree+".root");
            settings.push_back(trainingsNames.at(i));
        }

        settings.push_back(sig+"_v_"+bkg+"_"+tree+"_"+pt);

        TString extra = "#scale[0.8]{#splitline{"+treeNames.at(itree)+"}{"+ptNames.at(l)+"}}";
        extra = "#splitline{"+sig+" vs. "+bkg+"}{"+extra+"}";
        extra = "#splitline{"+smearName+"}{"+extra+"}";
        extraTexts.push_back(extra);

        colorSettings.push_back(trainings.size()-1);
        settingsMap.push_back(settings);
    }}}}


    // Generate plots for all input settings
    for(size_t i=0; i<settingsMap.size(); i++) {
        //for(size_t j=0; j < settingsMap.at(i).size(); j++) {
        //    std::cout<<settingsMap.at(i).at(j)<<std::endl;
        //}
        iCMax = colorSettings.at(i); 
        extraText=extraTexts.at(i);
        //std::cout<<iCMax<<std::endl;
        plotROC(settingsMap.at(i));
        gROOT->CloseFiles();
    }
}
