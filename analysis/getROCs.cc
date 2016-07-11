#include "./plotROC.C"

void getROCs() {

    std::vector<TString> type = { "q", "t", "g", "W", "Z" };


    for(size_t i=0; i<type.size(); i++) {
        TString ttype  = type.at(i);

        iCMax=2;

        for(size_t j=0; j<type.size();j++) {
            if(i==j) continue;
            TString ttype2=type.at(j);
            
            //plotROC("./eosrootfiles/shapesonly_tauratio/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_allpar.root", ttype+" vs "+ttype2+" allpar, tau ratios", 
            //        "./eosrootfiles/shapesonly/MVA_bdtg_"  +ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_allpar.root", ttype+" vs "+ttype2+" allpar, all tau",
            //        "./eosrootfiles/shapesonly_tauratio/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_allpar.root", ttype+" vs "+ttype2+" allpar 5 TeV, tau ratios", 
            //        "./eosrootfiles/shapesonly/MVA_bdtg_"  +ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_allpar.root", ttype+" vs "+ttype2+" allpar 5 TeV, all tau",  
            //        ttype+"_v_"+ttype2+"_shapes_tau_cmp_allpar");
       
            //plotROC("./eosrootfiles/shapesonly_tauratio/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tragam.root", ttype+" vs "+ttype2+" tragam, tau ratios", 
            //        "./eosrootfiles/shapesonly/MVA_bdtg_"  +ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tragam.root", ttype+" vs "+ttype2+" tragam, all tau",
            //        "./eosrootfiles/shapesonly_tauratio/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tragam.root", ttype+" vs "+ttype2+" tragam 5 TeV, tau ratios", 
            //        "./eosrootfiles/shapesonly/MVA_bdtg_"  +ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tragam.root", ttype+" vs "+ttype2+" tragam 5 TeV, all tau",  
            //        ttype+"_v_"+ttype2+"_shapes_tau_cmp_tragam");
       
            plotROC("./eosrootfiles/shapesonly_tauratio/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tracks.root", ttype+" vs "+ttype2+" tracks, tau ratios", 
                    "./eosrootfiles/shapesonly/MVA_bdtg_"  +ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tracks.root", ttype+" vs "+ttype2+" tracks, all tau",
                    "./eosrootfiles/shapesonly_tauratio/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tracks.root", ttype+" vs "+ttype2+" tracks 5 TeV, tau ratios", 
                    "./eosrootfiles/shapesonly/MVA_bdtg_"  +ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tracks.root", ttype+" vs "+ttype2+" tracks 5 TeV, all tau",  
                    ttype+"_v_"+ttype2+"_shapes_tau_cmp_tracks");
            break;
        }
        // iCMax=4;
       // // 1 allpar tvall ptall
       // plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_allpar.root", ttype+" vs g allpar", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_allpar.root", ttype+" vs W allpar",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_allpar.root", ttype+" vs Z allpar",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_allpar.root", ttype+" vs t allpar",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_allpar.root", ttype+" vs q allpar",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_allpar.root", ttype+" vs g allpar 5 TeV", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_allpar.root", ttype+" vs W allpar 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_allpar.root", ttype+" vs Z allpar 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_allpar.root", ttype+" vs t allpar 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_allpar.root", ttype+" vs q allpar 5 TeV",  ttype+"_v_all_allpar_ptall");

       // // 1 tragam tvall ptall
       // plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_tragam.root", ttype+" vs g tragam", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_tragam.root", ttype+" vs W tragam",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_tragam.root", ttype+" vs Z tragam",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_tragam.root", ttype+" vs t tragam",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_tragam.root", ttype+" vs q tragam",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_tragam.root", ttype+" vs g tragam 5 TeV", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_tragam.root", ttype+" vs W tragam 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_tragam.root", ttype+" vs Z tragam 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_tragam.root", ttype+" vs t tragam 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_tragam.root", ttype+" vs q tragam 5 TeV", ttype+"_v_all_tragam_ptall");

       // // 1 tracks tvall ptall
       // plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_tracks.root", ttype+" vs g tracks", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_tracks.root", ttype+" vs W tracks",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_tracks.root", ttype+" vs Z tracks",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_tracks.root", ttype+" vs t tracks",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_tracks.root", ttype+" vs q tracks", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_tracks.root", ttype+" vs g tracks 5 TeV", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_tracks.root", ttype+" vs W tracks 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_tracks.root", ttype+" vs Z tracks 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_tracks.root", ttype+" vs t tracks 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_tracks.root", ttype+" vs q tracks 5 TeV",ttype+"_v_all_tracks_ptall");


       // iCMax=colors.size();
       // // 1 allpar tvall
       // plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_allpar.root", ttype+" vs g allpar", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_allpar.root", ttype+" vs W allpar",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_allpar.root", ttype+" vs Z allpar",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_allpar.root", ttype+" vs t allpar",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_allpar.root", ttype+" vs q allpar", ttype+"_v_all_allpar");

       // // 1 tragam tvall
       // plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_tragam.root", ttype+" vs g tragam", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_tragam.root", ttype+" vs W tragam",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_tragam.root", ttype+" vs Z tragam",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_tragam.root", ttype+" vs t allpar",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_tragam.root", ttype+" vs q tragam", ttype+"_v_all_tragam");

       // // 1 tracks tvall
       // plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_gg-pt1_t_tracks.root", ttype+" vs g tracks", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_WW-pt1_t_tracks.root", ttype+" vs W tracks",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_ZZ-pt1_t_tracks.root", ttype+" vs Z tracks",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_tt-pt1_t_tracks.root", ttype+" vs t tracks",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_qq-pt1_t_tracks.root", ttype+" vs q tracks", ttype+"_v_all_tracks");

       // // 5 allpar tvall
       // plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_allpar.root", ttype+" vs g allpar 5 TeV", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_allpar.root", ttype+" vs W allpar 5 TeV", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_allpar.root", ttype+" vs Z allpar 5 TeV", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_allpar.root", ttype+" vs t allpar 5 TeV", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_allpar.root", ttype+" vs q allpar 5 TeV", ttype+"_v_all_allpar_5");

       // // 5 tragam tvall
       // plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_tragam.root", ttype+" vs g tragam 5 TeV",  
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_tragam.root", ttype+" vs W tragam 5 TeV", 
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_tragam.root", ttype+" vs Z tragam 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_tragam.root", ttype+" vs t tragam 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_tragam.root", ttype+" vs q tragam 5 TeV", ttype+"_v_all_tragam_5");

       // // 5 tracks tvall
       // plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_gg-pt5_t_tracks.root", ttype+" vs g tracks 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_WW-pt5_t_tracks.root", ttype+" vs W tracks 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_ZZ-pt5_t_tracks.root", ttype+" vs Z tracks 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_tt-pt5_t_tracks.root", ttype+" vs t tracks 5 TeV",
       //         "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_qq-pt5_t_tracks.root", ttype+" vs q tracks 5 TeV", ttype+"_v_all_tracks_5");

       // iCMax = 3;
       // for(size_t j=0; j<type.size(); j++) {
       //     if(i==j) continue;
       //     TString ttype2 = type.at(j);

       //     // t v W
       //     plotROC("./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tracks.root", ttype+" vs "+ttype2+" tracks", 
       //             "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_tragam.root", ttype+" vs "+ttype2+" tragam",
       //             "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt1_"+ttype2+ttype2+"-pt1_t_allpar.root", ttype+" vs "+ttype2+" allpar",
       //             "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tracks.root", ttype+" vs "+ttype2+" tracks 5 TeV",
       //             "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_tragam.root", ttype+" vs "+ttype2+" tragam 5 TeV",
       //             "./eosrootfiles/MVA_bdtg_"+ttype+ttype+"-pt5_"+ttype2+ttype2+"-pt5_t_allpar.root", ttype+" vs "+ttype2+" allpar 5 TeV", ttype+"_v_"+ttype2);


       // }
    }
}
