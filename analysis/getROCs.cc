/****
 * An example macro showcasing the plotROC.C code
 * NOTE: not runnable "out of the box" 
 */

#include "./plotROC.C"

void getROCs() {
    
    // analysis information 
    std::vector<TString> type = {"t", "W", "Z", "q", "g" };
    std::vector<TString> trees = { "allpar", "tracks", "tragam" };
    std::vector<TString> trainings = { "_shapesonly", "_massonly", "_all" };
    TString plotDir = "trainings_07152016_all";

    for(size_t i=0; i<type.size(); i++) {
        TString ttype  = type.at(i);
    
        // repeat colorscheme after 3 entries
        iCMax=3;

        for(size_t j=0; j<type.size();j++) {
            if(i==j) continue;
            TString ttype2=type.at(j);
           
            // compare all trainings for each tree 
            for(size_t k=0; k < trees.size(); k++) {
                TString ttree = trees.at(k);
                plotROC("./eosrootfiles/trainings_07142016_massonly/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_"+ttree+".root", 
                            ttype+" vs "+ttype2+" "+ttree+", mass only", 
                        "./eosrootfiles/trainings_07142016_shapesonly/MVA_bdtg_"  +ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_"+ttree+".root", 
                            ttype+" vs "+ttype2+" "+ttree+", shapes only",
                        "./eosrootfiles/trainings_07142016_all/MVA_bdtg_"  +ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_"+ttree+".root", 
                            ttype+" vs "+ttype2+" "+ttree+", shapes + mass",
                       // "./eosrootfiles/trainings_07142016_massonly/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_"+ttree+".root", 
                       //     ttype+" vs "+ttype2+" "+ttree+", mass only 5 TeV", 
                       // "./eosrootfiles/trainings_07142016_shapesonly/MVA_bdtg_"  +ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_"+ttree+".root", 
                       //     ttype+" vs "+ttype2+" "+ttree+", shapes only 5 TeV",
                       // "./eosrootfiles/trainings_07142016_all/MVA_bdtg_"  +ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_"+ttree+".root", 
                       //     ttype+" vs "+ttype2+" "+ttree+", shapes + mass 5 TeV",
                        ttype+"_v_"+ttype2+"_cmppt1_"+ttree);
            }

            // compare all trees for each training
            for(size_t k=0; k < trainings.size(); k++) {
                TString ttrain = trainings.at(k);
                plotROC("./eosrootfiles/trainings_07142016"+ttrain+"/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_allpar.root", 
                            ttype+" vs "+ttype2+" allpar, ", 
                        "./eosrootfiles/trainings_07142016"+ttrain+"/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tragam.root", 
                            ttype+" vs "+ttype2+" tragam, ",
                        "./eosrootfiles/trainings_07142016"+ttrain+"/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tracks.root", 
                            ttype+" vs "+ttype2+" tracks, ",
                       // "./eosrootfiles/trainings_07142016"+ttrain+"/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_allpar.root", 
                       //     ttype+" vs "+ttype2+" allpar, 5 TeV", 
                       // "./eosrootfiles/trainings_07142016"+ttrain+"/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tragam.root", 
                       //     ttype+" vs "+ttype2+" tragam, 5 TeV",
                       // "./eosrootfiles/trainings_07142016"+ttrain+"/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tracks.root", 
                       //     ttype+" vs "+ttype2+" tracks, 5 TeV",
                        ttype+"_v_"+ttype2+"_trainingpt1"+ttrain);
            }
       
        }

        // repeat colorscheme after 4 entries
        iCMax=4;

        // 1 allpar tvall ptall
        plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_allpar.root", ttype+" vs g allpar", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_allpar.root", ttype+" vs W allpar",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_allpar.root", ttype+" vs Z allpar",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_allpar.root", ttype+" vs t allpar",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_allpar.root", ttype+" vs q allpar",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_allpar.root", ttype+" vs g allpar 5 TeV", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_allpar.root", ttype+" vs W allpar 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_allpar.root", ttype+" vs Z allpar 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_allpar.root", ttype+" vs t allpar 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_allpar.root", ttype+" vs q allpar 5 TeV",  ttype+"_v_all_allpar_ptall");

        // 1 tragam tvall ptall
        plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_tragam.root", ttype+" vs g tragam", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_tragam.root", ttype+" vs W tragam",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_tragam.root", ttype+" vs Z tragam",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_tragam.root", ttype+" vs t tragam",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_tragam.root", ttype+" vs q tragam",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_tragam.root", ttype+" vs g tragam 5 TeV", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_tragam.root", ttype+" vs W tragam 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_tragam.root", ttype+" vs Z tragam 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_tragam.root", ttype+" vs t tragam 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_tragam.root", ttype+" vs q tragam 5 TeV", ttype+"_v_all_tragam_ptall");

        // 1 tracks tvall ptall
        plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_tracks.root", ttype+" vs g tracks", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_tracks.root", ttype+" vs W tracks",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_tracks.root", ttype+" vs Z tracks",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_tracks.root", ttype+" vs t tracks",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_tracks.root", ttype+" vs q tracks", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_tracks.root", ttype+" vs g tracks 5 TeV", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_tracks.root", ttype+" vs W tracks 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_tracks.root", ttype+" vs Z tracks 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_tracks.root", ttype+" vs t tracks 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_tracks.root", ttype+" vs q tracks 5 TeV",
                ttype+"_v_all_tracks_ptall");


        iCMax=colors.size();
        // 1 allpar tvall
        plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_allpar.root", ttype+" vs g allpar", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_allpar.root", ttype+" vs W allpar",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_allpar.root", ttype+" vs Z allpar",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_allpar.root", ttype+" vs t allpar",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_allpar.root", ttype+" vs q allpar", 
                ttype+"_v_all_allpar");

        // 1 tragam tvall
        plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_tragam.root", ttype+" vs g tragam", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_tragam.root", ttype+" vs W tragam",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_tragam.root", ttype+" vs Z tragam",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_tragam.root", ttype+" vs t allpar",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_tragam.root", ttype+" vs q tragam", 
                ttype+"_v_all_tragam");

        // 1 tracks tvall
        plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_tracks.root", ttype+" vs g tracks", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_tracks.root", ttype+" vs W tracks",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_tracks.root", ttype+" vs Z tracks",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_tracks.root", ttype+" vs t tracks",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_tracks.root", ttype+" vs q tracks", 
                ttype+"_v_all_tracks");

        // 5 allpar tvall
        plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_allpar.root", ttype+" vs g allpar 5 TeV", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_allpar.root", ttype+" vs W allpar 5 TeV", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_allpar.root", ttype+" vs Z allpar 5 TeV", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_allpar.root", ttype+" vs t allpar 5 TeV", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_allpar.root", ttype+" vs q allpar 5 TeV", 
                ttype+"_v_all_allpar_5");

        // 5 tragam tvall
        plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_tragam.root", ttype+" vs g tragam 5 TeV",  
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_tragam.root", ttype+" vs W tragam 5 TeV", 
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_tragam.root", ttype+" vs Z tragam 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_tragam.root", ttype+" vs t tragam 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_tragam.root", ttype+" vs q tragam 5 TeV", 
                ttype+"_v_all_tragam_5");

        // 5 tracks tvall
        plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_tracks.root", ttype+" vs g tracks 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_tracks.root", ttype+" vs W tracks 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_tracks.root", ttype+" vs Z tracks 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_tracks.root", ttype+" vs t tracks 5 TeV",
                "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_tracks.root", ttype+" vs q tracks 5 TeV", 
                ttype+"_v_all_tracks_5");

        iCMax = 3;
        for(size_t j=0; j<type.size(); j++) {
            if(i==j) continue;
            TString ttype2 = type.at(j);

            // t v W
            plotROC("./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tracks.root", ttype+" vs "+ttype2+" tracks", 
                    "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tragam.root", ttype+" vs "+ttype2+" tragam",
                    "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_allpar.root", ttype+" vs "+ttype2+" allpar",
                    "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tracks.root", ttype+" vs "+ttype2+" tracks 5 TeV",
                    "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tragam.root", ttype+" vs "+ttype2+" tragam 5 TeV",
                    "./eosrootfiles/"+plotDir+"/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_allpar.root", ttype+" vs "+ttype2+" allpar 5 TeV", 
                    ttype+"_v_"+ttype2);


        }
    }
}
