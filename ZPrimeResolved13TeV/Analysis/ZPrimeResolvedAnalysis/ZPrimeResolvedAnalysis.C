/////////////////////////////////////////////////////////////
//// ZPrimeResolvedAnalysis code
//// Author: ATLAS Collaboration (2021)
////
////
//// DISCLAIMER:
//// This Software is intended for educational use only!
//// Under no circumstances does it qualify to reproduce actual ATLAS analysis results or produce publishable results!
/////////////////////////////////////////////////////////////

#define ZPrimeResolvedAnalysis_cxx
#include "ZPrimeResolvedAnalysis.h"
#include "ZPrimeResolvedAnalysisHistograms.h"
#include <iostream>
#include <cstring>
#include <string>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>  ////!!!!!!!!!!!!!!!!!!!!

string name;

void ZPrimeResolvedAnalysis::Begin(TTree * )
{

  nEvents=0;

}

void ZPrimeResolvedAnalysis::SlaveBegin(TTree * )
{
  TString option = GetOption();
  printf("Starting analysis with process option: %s \n", option.Data());
  
  name=option;
  
  define_histograms();
  
  FillOutputList();
    
}

Bool_t ZPrimeResolvedAnalysis::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
    
  fChain->GetTree()->GetEntry(entry);
  //  int cut1_mc = 0;

  if(fChain->GetTree()->GetEntries()>0)
  {
    //Do analysis

    //SF
    Float_t scaleFactor = scaleFactor_ELE*scaleFactor_MUON*scaleFactor_TRIGGER;
    //EventW
    Float_t eventWeight = mcWeight*scaleFactor_PILEUP*scaleFactor_ZVERTEX;
    //weight = SF * EventW
    Float_t weight = scaleFactor*eventWeight;

    // Make difference between data and MC
    if (weight == 0. && channelNumber != 110090 && channelNumber != 110091) weight = 1.;

    // Missing Et of the event in GeV ----------------------------------------!!!!!
    Float_t missingEt = met_et/1000.;


    // Preselection cut for electron/muon trigger, Good Run List, and good vertex
    if(trigE || trigM)
    {
      if(passGRL)
      {
        if(hasGoodVertex)
        {
          //Find the good leptons
          int goodlep_index = 0;
          int goodlep_n = 0;
          int lep_index = 0;

          for(int i=0; i<lep_n; i++)
          {
            if(lep_pt->at(i)>25000. && (lep_ptcone30->at(i)/lep_pt->at(i)) < 0.15 && (lep_etcone20->at(i)/lep_pt->at(i)) < 0.15 )
            {
                // electron selection in fiducial region excluding candidates in the transition region between the barrel and endcap
                // electromagnetic calorimeters
                if( lep_type->at(i)==11 && TMath::Abs(lep_eta->at(i)) < 2.47 && ( TMath::Abs(lep_eta->at(i)) < 1.37 || TMath::Abs(lep_eta->at(i)) > 1.52))
                {
                  goodlep_n++;
                  goodlep_index = i;  
                  lep_index++; 
                }
                
                if( lep_type->at(i) == 13 && TMath::Abs(lep_eta->at(i)) < 2.5 ) 
                {
                  goodlep_n++;
                  goodlep_index = i;  
                  lep_index++; 
                }
              }
          }

            
        //Zero cut
        if(goodlep_n==1)
        {
            /* I think this loop is not necessary because there was find the goodlepton_index above ----------------------------------
            for(int i=0; i<lep_n; i++)
            {
                if(lep_pt->at(i)>25000. && (lep_ptcone30[i]/lep_pt[i]) < 0.15 && (lep_etcone20[i]/lep_pt[i]) < 0.15 )
                {
                    goodlep_index = i;
                }
            }
            */
            
            
            // TLorentzVector definitions
            TLorentzVector Lepton_1  = TLorentzVector();
            TLorentzVector      MeT  = TLorentzVector();
            TLorentzVector   Lepton1_MeT = TLorentzVector();

            Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep_index), lep_eta->at(goodlep_index), lep_phi->at(goodlep_index),lep_E->at(goodlep_index));
            MeT.SetPtEtaPhiE(met_et, 0, met_phi , met_et);


            //Calculation of the Invariant Mass using TLorentz vectors (First Lepton + MeT)-------Transverse Mass? 
            Lepton1_MeT = Lepton_1 + MeT;
            float InvMass1       = Lepton1_MeT.Mag();
            float InvMass1_inGeV = InvMass1/1000.;
            float Lepton1_MeT_MT = sqrt(2*Lepton_1.Pt()*MeT.Et()*(1-cos(Lepton_1.DeltaPhi(MeT))));

            /* I think it is not necessary because when we applied cuts for goodlepton, we used these selections
            //First cut : Exactly one good leptons with pT>25GeV
            if(goodlep_n ==1 && lep_pt[goodlep_index] >25000. && (lep_ptcone30[goodlep_index]/lep_pt[goodlep_index]) < 0.15 && (lep_etcone20[goodlep_index]/lep_pt[goodlep_index]) < 0.15)
            {
            */
            
            //Preselection of good jets
            int goodjet_n = 0;
            int goodbjet_n = 0;
            int goodjet_index[jet_n];
            int jet_index = 0;    
            int goodbjet_index[jet_n];
            int bjet_index = 0;

            if(jet_n >= 4) //At least 4 jets -----------------------------HERE!!!!!!!!!!!!!!!!!!!!!!!!!---------------------------
            {
                std::vector<Float_t*> jet_variable={jet_pt,jet_eta,jet_phi,jet_E,jet_m,jet_jvf,jet_SV0,jet_MV1};
                int pass_jet = 0;
                for(int i = 0;i < jet_n;i++)
                {
                  if(jet_pt[i] > 25000. && TMath::Abs(jet_eta[i]) < 2.5)
                  {
                      // JVF cleaning
                      bool jvf_pass=true;
                      if (jet_pt[i] < 50000. && TMath::Abs(jet_eta[i]) < 2.4 && jet_jvf[i] < 0.5) jvf_pass=false;
                      if (jvf_pass) 
                      {
                          goodjet_n++;
                          goodjet_index[jet_index] = i;
                          jet_index++;

			  // cut on 0.7982 is 70% WP
			  if (jet_MV1[i] >0.7982 )
			    {
			      goodbjet_n++;
			      goodbjet_index[bjet_index] = i;
			      bjet_index++;
			    }
			}
		    }
		}
		
		
		if(goodjet_n >= 4)
		  {
		    //At least two b-tagged jets
		    if(goodbjet_n >= 2)
		      {
			
			if(met_et > 30000.)
			  {
			    
			    if(Lepton1_MeT_MT > 30000.)
			      {
                    // Finding 3 hardest pT sum
				int fijet = 0;
				int sejet = 0;
				int thjet = 0;
				int fojet = 0;
				float hardest_pt_sum = 0;
				float hardest3mass = 0;


				//Select these jets using TLorentzVector
				for(int i = 0; i < goodjet_n; i++){
				  for(int j = i+1; j < goodjet_n; j++){
				    for(int k = j+1; k < goodjet_n; k++){
				      for(int h = k+1; h < goodjet_n; h++){
					TLorentzVector jet1 = TLorentzVector ();
					TLorentzVector jet2 = TLorentzVector ();
					TLorentzVector jet3 = TLorentzVector ();
					TLorentzVector jet4 = TLorentzVector ();
					jet1.SetPtEtaPhiE(jet_pt[goodjet_index[i]], jet_eta[goodjet_index[i]], jet_phi[goodjet_index[i]], jet_E[goodjet_index[i]]);
					jet2.SetPtEtaPhiE(jet_pt[goodjet_index[j]], jet_eta[goodjet_index[j]], jet_phi[goodjet_index[j]], jet_E[goodjet_index[j]]);
					jet3.SetPtEtaPhiE(jet_pt[goodjet_index[k]], jet_eta[goodjet_index[k]], jet_phi[goodjet_index[k]], jet_E[goodjet_index[k]]);
					jet4.SetPtEtaPhiE(jet_pt[goodjet_index[h]], jet_eta[goodjet_index[h]], jet_phi[goodjet_index[h]], jet_E[goodjet_index[h]]);
					
					float pt_sum = (jet1 + jet2 + jet3).Pt()/1000. ;
					
					if (pt_sum > hardest_pt_sum){
					  hardest_pt_sum=pt_sum;
					  hardest3mass = (jet1 + jet2 + jet3).M()/1000.;
					  fijet=i,sejet=j,thjet=k,fojet=h;

					  float pt2_sum12 = (jet1 + jet2).Pt()/1000.;
					  float pt2_sum13 = (jet1 + jet3).Pt()/1000.;
					  float pt2_sum23 = (jet2 + jet3).Pt()/1000.;
					  float maximum = std::max(pt2_sum12,std::max(pt2_sum13,pt2_sum23));
					  
					  if(maximum==pt2_sum13){
					    fijet=i,sejet=k,thjet=j;
					  }
					  if(maximum==pt2_sum23){
					    fijet=k,sejet=j,thjet=i;
					  }
					  
					}
				      }
				    }
				  }
				}


				TLorentzVector j1 = TLorentzVector ();
				TLorentzVector j2 = TLorentzVector ();
				j1.SetPtEtaPhiE(jet_pt[goodjet_index[fijet]], jet_eta[goodjet_index[fijet]], jet_phi[goodjet_index[fijet]], jet_E[goodjet_index[fijet]]);
				j2.SetPtEtaPhiE(jet_pt[goodjet_index[sejet]], jet_eta[goodjet_index[sejet]], jet_phi[goodjet_index[sejet]], jet_E[goodjet_index[sejet]]);
				
				float hardest2mass=(j1+j2).M()/1000.;
				
				if(hardest2mass > 50 && hardest2mass < 120){
				  hist_Wmass->Fill(hardest2mass, weight);}
				
				if(hardest3mass > 100 && hardest3mass < 250){
				  hist_Topmass->Fill(hardest3mass, weight);}
				
				
				TLorentzVector j3 = TLorentzVector ();
				TLorentzVector j4 = TLorentzVector ();
				j1.SetPtEtaPhiE(jet_pt[goodjet_index[thjet]], jet_eta[goodjet_index[thjet]], jet_phi[goodjet_index[thjet]], jet_E[goodjet_index[thjet]]);
				j2.SetPtEtaPhiE(jet_pt[goodjet_index[fojet]], jet_eta[goodjet_index[fojet]], jet_phi[goodjet_index[fojet]], jet_E[goodjet_index[fojet]]);

				float mazz = (j1 + j2 + j3 + j4 + Lepton_1 + MeT).M()/1000.;
				hist_mazz->Fill(mazz, weight);
                    
				double names_of_global_variable[]={InvMass1_inGeV, missingEt, vxp_z, (double)pvxp_n, Lepton1_MeT_MT/1000.,hardest3mass,hardest2mass,mazz};
				double names_of_leadlep_variable[]={Lepton_1.Pt()/1000., Lepton_1.Eta(), Lepton_1.E()/1000., Lepton_1.Phi(), lep_charge[goodlep_index], (double)lep_type[goodlep_index], lep_ptcone30[goodlep_index], lep_etcone20[goodlep_index], lep_z0[goodlep_index], lep_trackd0pvunbiased[goodlep_index]};
				double names_of_jet_variable[]={(double)jet_n, jet_pt[0]/1000., jet_eta[0], jet_m[0]/1000., jet_jvf[0], jet_MV1[0]};
				
				TString histonames_of_global_variable[]={"hist_vismass","hist_etmiss","hist_vxp_z","hist_pvxp_n", "hist_mt","hist_Topmass","hist_Wmass","hist_mazz"};
				TString histonames_of_leadlep_variable[]={"hist_leadleptpt", "hist_leadlepteta","hist_leadleptE","hist_leadleptphi","hist_leadleptch","hist_leadleptID","hist_leadlept_ptc","hist_leadleptetc","hist_leadlepz0","hist_leadlepd0"};
				TString histonames_of_jet_variable[]={"hist_n_jets","hist_leadjet_pt","hist_leadjet_eta","hist_leadjet_m", "hist_leadjet_jvf", "hist_leadjet_MV1"};
				
				int length_global = sizeof(names_of_global_variable)/sizeof(names_of_global_variable[0]);
				int length_leadlep = sizeof(names_of_leadlep_variable)/sizeof(names_of_leadlep_variable[0]);
				int length_leadjet = sizeof(names_of_jet_variable)/sizeof(names_of_jet_variable[0]);
				
				for (int i=0; i<length_global; i++)
				  {
				    FillHistogramsGlobal( names_of_global_variable[i], weight, histonames_of_global_variable[i]);
				  }
				for (int i=0; i<length_leadlep; i++)
				  {
				    FillHistogramsLeadlept( names_of_leadlep_variable[i], weight, histonames_of_leadlep_variable[i]);
				  }
				for (int i=0; i<length_leadjet; i++)
				  {
				    FillHistogramsLeadJet( names_of_jet_variable[i], weight, histonames_of_jet_variable[i]);
				  }
			      }
			  }
		      }
		  }
	      }
	    }
	  }
	}
      }
    }
  }
  //  std::cout<<cut1_mc<<std::endl;
  return kTRUE;
}

void ZPrimeResolvedAnalysis::SlaveTerminate()
{
    
}

void ZPrimeResolvedAnalysis::Terminate()
{
  TString filename_option = GetOption();
  printf("Writting with name option: %s \n", filename_option.Data());
  TString output_name="Output_ZPrimeResolvedAnalysis/"+filename_option+".root";
  const char* filename = output_name;

  TFile physicsoutput_ZPrimeResolved(filename,"recreate");
  WriteHistograms();
  physicsoutput_ZPrimeResolved.Close();

}
