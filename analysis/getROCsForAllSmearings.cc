/****
 * An example macro showcasing the plotROC.C code
 */

#include "./plotROC.C"

void getROCsForAllSmearings(TString inputDir,TString outDir="./") {
    
    // analysis information 
    std::vector<TString> type = {"t", "W", "Z", "q", "g" };
    std::vector<TString> trees = { "allpar", "tracks", "tragam" };
    std::vector<TString> treeNames = { "All particles", "Tracks only", "Tracks+#gamma" };
    std::vector<TString> pts = { "pt1", "pt5" };
    std::vector<TString> ptNames = { "p_{T} 1 TeV", "p_{T} 5 TeV" };
    std::vector<TString> trainings = { "shapesonly", "massonly", "all" };
    std::vector<TString> trainingsNames = { "Shapes only", "Mass only", "All variables" };
    std::vector<TString> smears = { "r0_h0_e0", 
                                    "r05_h05_e005", 
                                    "r05_h01_e005", 
                                    "r05_h01_e005_t", 
                                    "r05_h002_e005_t" };
    std::vector<TString> smearNames = { "Perfect", 
                                        "HCAL0.05", 
                                        "HCAL0.01", 
                                        "Tracker degradation", 
                                        "Extreme granularity" };
    outputDir=outDir+"/";

    
    // setup vectors of settings
    std::vector<TString> extraTexts;
    std::vector<int> colorSettings;
    std::vector<std::vector<TString>> settingsMap;


    // all-smears for fixed tree, training, sig/bkg, pt
    for(size_t i=0; i<trainings.size(); i++) {
    for(size_t j=0; j<type.size();  j++) {
    for(size_t k=0; k<type.size();  k++) {
        if(j==k) continue;
    for(size_t l=0; l<pts.size();   l++) {
    for(size_t itree=0; itree < trees.size(); itree++) {
        TString training=trainings.at(i);
        TString sig = type.at(j);
        TString bkg = type.at(k);
        TString pt  = pts.at(l);
        TString tree = trees.at(itree);

        std::vector<TString> settings={};
        
        for(size_t ismear=0; ismear<smears.size(); ismear++) {
           TString smear=smears.at(ismear);

           settings.push_back(inputDir+"/"+smear+"/eosrootfiles/"+training+"/MVA_bdtg_"
                                      +sig+sig+"-"+pt+"_"+bkg+bkg
                                      +"-"+pt+"_t_"+tree+".root");
           settings.push_back(smearNames.at(ismear));
        }

        settings.push_back("cmpsmear_"+sig+"_v_"+bkg+"_"+tree+"_"+pt+"_"+training);

        TString extra = "#splitline{"+trainingsNames.at(i)+"}{"+ptNames.at(l)+"}";
        extra = "#scale[0.8]{#splitline{"+treeNames.at(itree)+"}{"+extra+"}}";
        extra = "#splitline{"+sig+" vs. "+bkg+"}{"+extra+"}";
        extraTexts.push_back(extra);

        colorSettings.push_back(smears.size());
        settingsMap.push_back(settings);
    }}}}}


    // Generate plots for all input settings
    for(size_t i=0; i<settingsMap.size(); i++) {
        iCMax = colorSettings.at(i); 
        extraText=extraTexts.at(i);
        plotROC(settingsMap.at(i));
        gROOT->CloseFiles();
    }
}
