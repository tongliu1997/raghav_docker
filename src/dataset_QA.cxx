//! Run 14 AuAu Dataset QA code
//! takes low/mid/high lumi and using HT or MB trigger
//!
//! Using Nick Elsey's jetreader class
//! Raghav Kunnawalkam Elayavalli
//! June 12th 2019
//! Wayne State University, Detroit
//! raghavke@wayne.edu 
//!

#include "jetreader/lib/assert.h"
//#include "jetreader/lib/test_data.h"
#include "jetreader/reader/reader.h"
#include "jetreader/lib/parse_csv.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>     // for getenv, atof, atoi
#include <vector>
#include <assert.h>
#include <cmath>
#include <climits>
#include <sstream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"

#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StPicoEvent/StPicoTrack.h"

#include "jetreader/reader/centrality_def.h"
#include "jetreader/reader/centrality.h"
#include "jetreader/reader/bemc_helper.h"

using namespace std;
using namespace fastjet;
// using namespace contrib;
using namespace jetreader;

#include <random>
#include <algorithm>

int main( int argc, const char** argv ){

  Long64_t Ntot = -1;//! total number of events to run over
  std::string lumi = "low"; //! '', low, mid, high (nothing stands for pre-split)
  std::string picofile = "/home/tl543/raghav_docker/testfile.list"; //! input file name
//  std::string picofile = "/home/tl543/Isobar_Data/st_physics_19124006_raw_4000043.picoDst.root"; //! input file name
  std::string outfile = "outputfile.root"; //! output file name 
  double trackDCACut = 3.0; //! primary track DCA cut
  double trackNhitCut = 15; //! primary track number of hits in tpc
  double trackNhitFracCut = 0.52; //! primary track nhits fraction to remove split tracks
  std::string system="Ru";  
  bool printDebug = false;
  bool removebadruns = false;
  bool removebadtowers = false;
  // Parse arguments
  // ---------------
  vector<string> arguments(argv + 1, argv + argc);
  bool argsokay=true;
  for ( auto parg = arguments.begin() ; parg!=arguments.end() ; ++parg){
    string arg=*parg;
    if ( arg == "-N" ){      
      if ( ++parg == arguments.end() ){ argsokay=false; break; }
      Ntot = atoi( parg->data());
    } else if ( arg == "-lumi" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      lumi = *parg;
    } else if ( arg == "-system" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      system = *parg;
    } else if ( arg == "-picofile" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      picofile = *parg;
    } else if ( arg == "-outfile" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      outfile = *parg;
    } else if ( arg == "-trackDCACut" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      trackDCACut = atof( parg->data());
    } else if ( arg == "-trackNhitCut" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      trackNhitCut = atof( parg->data());
    } else if ( arg == "-trackNhitFracCut" ){      
      if (++parg==arguments.end() ){ argsokay=false; break; }
      trackNhitFracCut = atof( parg->data());
    } else {
      argsokay=false;
      break;
    }
  }
 
  if ( !argsokay ) {
    cerr << "usage: " << argv[0] << endl
      	 << " [-N Nevents (<0 for all)]" << endl
	 << " [-lumi low, mid, high or combined]" << endl
	 << " [-trig HT or MB]" << endl
	 << " [-picofile input file name including path ]" << endl
	 << " [-outfile output file name including path ]" << endl
	 << " [-trackDCACut Primary track DCA cut, nominal value = 1.0]" << endl
	 << " [-trackNhitCut Primary track number of hits, nominal value = 15]" << endl
	 << " [-trackNhitFracCut Primary track Nhit fraction, nominal value = 0.52]" << endl;
    throw std::runtime_error("Not a valid list of options");
  }


  //! declare the output file and histograms
  TFile out(Form("%s", outfile.c_str()), "RECREATE");

  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  TH2D * htrackpT_eta = new TH2D("htrackpT_eta", "", 600, 0, 30, 10, -1, 1);
  TH2D * htrackpT_phi = new TH2D("htrackpT_phi", "", 600, 0, 30, 630, -6.3, 6.3);
  TH2D * htrackpT_cent = new TH2D("htrackpT_cent","",600,0,30,18,-1.5,16.5);
  TH2D * hntrk_cent = new TH2D("hntrk_cent","",300,0,300,18,-1.5,16.5);
//  TTree * evtVars = new TTree ("evtVars", "Event variables to build correlations");
//  double val_vz;
//  evtVars->Branch("vz",&val_vz,"vz/D");  
//  long val_zdcx;
//  evtVars->Branch("zdcx",&val_zdcx,"zdcx/I");  
//  long val_refmult;
//  evtVars->Branch("refmult",&val_refmult,"refmult/I");  
//  vector<double> track_pT;
//  evtVars->Branch("track_pT",&track_pT);
//  vector<double> track_eta;
//  evtVars->Branch("track_eta",&track_eta);
//  vector<double> track_phi;
//  evtVars->Branch("track_phi",&track_phi);
//  vector<double> track_dca;
//  evtVars->Branch("track_dca",&track_dca);
//  vector<int> track_Nhit;
//  evtVars->Branch("track_Nhit",&track_Nhit);

  long event = 0;
  long read_event = 0;
  double weight_event = 0;
  string infile;
  ifstream list;
  list.open(picofile);
  while(getline(list,infile)){
//  if(1){
    // creating reader with test StPicoDst
    cout<<"input file = "<<infile<<endl;
    jetreader::Reader reader(infile);

    //! setup events
    reader.eventSelector()->setVzRange(-35, 25);
    reader.eventSelector()->setVrMax(2);
    reader.eventSelector()->setdVzMax(25);
    reader.eventSelector()->addTriggerId(600001);
    reader.eventSelector()->addTriggerId(600011);
    reader.eventSelector()->addTriggerId(600021);
    reader.eventSelector()->addTriggerId(600031);
    std::string badrunlist = "";
    if(removebadruns){
      badrunlist = "test.csv";
      reader.eventSelector()->addBadRuns(badrunlist);
    }
    
    //! setup tracks 
    reader.trackSelector()->setPtMax(30.0);
    reader.trackSelector()->setPtMin(0.2);
    reader.trackSelector()->rejectEventOnPtFailure(true);
    reader.trackSelector()->setDcaMax(trackDCACut);
    reader.trackSelector()->setNHitsMin(trackNhitCut);
    reader.trackSelector()->setNHitsFracMin(trackNhitFracCut);
    
    //! setup towers
    reader.useHadronicCorrection(true, 1.0);
    reader.useApproximateTrackTowerMatching(true);    
    std::string badtowerlist = "";
    if(removebadtowers){
      badtowerlist = "test.csv";
      reader.towerSelector()->addBadTowers(badtowerlist);
    }
    reader.towerSelector()->setEtMax(30.0);
    reader.towerSelector()->rejectEventOnEtFailure(false);

    //! setup reader centrality 
    if(system=="Ru") 
      reader.centrality().loadCentralityDef(jetreader::CentDefId::Run18Ru);
    else if(system=="Zr")  
      reader.centrality().loadCentralityDef(jetreader::CentDefId::Run18Zr);
    // turn off branches that don't exist in this tree (not necessary)
    reader.SetStatus("BEmcSmdEHit", false);
    reader.SetStatus("BEmcSmdPHit", false);
    
    // initialize the reader
    reader.Init();

    std::cout << "number of events in chain: " << reader.tree()->GetEntries()
	      << std::endl;

    
    while (reader.next()) {
      if (event % 1000 == 0)
	std::cout << "event: " << event << std::endl;
      event++;
      if(Ntot > 0 && event >= Ntot)
	break;

      

      int centbin = -1;
      double cent_weight = 0.0;
      double refm = reader.picoDst()->event()->refMult();
      
      if(printDebug)
       	cout<<"event centrality = "<<reader.centrality16()<<endl;
      centbin = reader.centrality16();
      cent_weight = reader.centrality().weight();
      weight_event += cent_weight; 
      // int cbin = -1;
      // for(int ic = 0; ic<ncentbins; ic++){
      // 	if(centbin >= centbins[ic])
      // 	  cbin = ic;
      // }
      
      // if(cbin < 0)
      // 	continue;

      long zdcx_val = reader.picoDst()->event()->ZDCx();

      double event_weight = 1.0;
      event_weight = cent_weight;
      TVector3 vertex = reader.picoDst()->event()->primaryVertex();
      bool use_primary_tracks_ = true;
      int ntracks = 0;
      for (int track_id = 0; track_id < reader.picoDst()->numberOfTracks(); ++track_id){
	StPicoTrack *track = (StPicoTrack*)reader.picoDst()->track(track_id);
	jetreader::TrackStatus track_status = reader.trackSelector()->select(track, vertex, use_primary_tracks_);
	if (track_status == jetreader::TrackStatus::acceptTrack){
	  ntracks++;
	  TVector3 momentum=track->pMom();
	  short track_charge=track->charge();
	  double track_pt=momentum.Perp();
	  double track_eta=momentum.Eta();
	  double track_phi=momentum.Phi();
	  
	  htrackpT_eta->Fill(track_pt, track_eta, event_weight/track_pt);
	  htrackpT_phi->Fill(track_pt, track_phi, event_weight/track_pt);
  	  htrackpT_cent->Fill(track_pt, centbin, event_weight/track_pt);

	}
      }
      
      // process tracks and towers now through the selected pseudojets
//      std::vector<fastjet::PseudoJet> container = reader.pseudojets();

//      int ntracks = 0, ntowers = 0;

//      track_pT.clear();
//      track_phi.clear();
//      track_eta.clear();
//      track_dca.clear();
//      track_Nhit.clear();
      // track_charge.clear();

/*    
      for (auto &track : container) {
	jetreader::VectorInfo track_info =
	  track.user_info<jetreader::VectorInfo>();
	if (track_info.isPrimary()) {
	  if(track.pt() >= 0.2){

	    int trackid=track_info.trackId();
	    cout<<ntracks<<"\t"<<trackid<<"\t"<<reader.picoDst()->numberOfTracks()<<"\t";
	    StPicoTrack *pico_track = (StPicoTrack*)reader.picoDst()->track(trackid);
	    int pico_id=pico_track->id();
 	    cout<<pico_id<<endl;
	    double pico_pt=pico_track->pMom().Perp();
	    if (pico_pt != track.pt()) cout<<"PT different! "<<pico_pt<<"\t"<<track.pt()<<endl;
	    ntracks++;
	    htrackpT_eta->Fill(track.pt(), track.eta(), event_weight);
	    htrackpT_phi->Fill(track.pt(), track.phi(), event_weight);
  	    htrackpT_cent->Fill(track.pt(),centbin,event_weight);
//	    track_pT.push_back(track.pt());
//	    track_eta.push_back(track.eta());
//	    track_phi.push_back(track.phi());
//	    track_dca.push_back(track_info.dca());
//	    track_Nhit.push_back(track_info.nhits());
	    // track_charge.push_back(track_info.charge());
	    
	    // unsigned matTow = track_info.matchedTower();
	    // if(matTow > 0 && matTow <= 4800){
	    //   double towerE = reader.picoDst()->btowHit(matTow-1)->energy();	    
	    //   hTrackpT_matTowerE[cbin]->Fill(track.pt(), towerE, event_weight);
	    // }
	    // hTrackpT_matTowerCorrE->Fill(track.pt(), towerEt-track.pt(), event_weight);
	  }
	}
	else if (track_info.isBemcTower()) {
	  if(track.Et() >= 0.2) {
	  }
	} else if (track_info.isGlobal()) {
	  std::cout << "global track?" << std::endl;
	} else {
	  JETREADER_THROW("Should not be here; vector labeled as not primary, "
			  "not global, and not bemc tower");
	}
      }
*/
      hntrk_cent->Fill(ntracks,centbin,event_weight);
//    evtVars->Fill();      
      
    }
  read_event+=reader.tree()->GetReadEntry()+1;  
  }//! running over all the files 
    std::cout << event << " good events out of " << read_event << "\n";
    std::cout << "accepted : " << (double) 100.* event / read_event<< "% of events\n";
    
    std::cout <<"Weighted events "<<weight_event<<endl; 

  out.Write();
  out.Close();

  return 0;
}
