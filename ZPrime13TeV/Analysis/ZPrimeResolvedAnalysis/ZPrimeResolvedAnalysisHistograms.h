////////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include <iostream>
// Methods
////////////////////////////////////////////////////////////////////////////////
void ZPrimeResolvedAnalysis::define_histograms()
{
  // Histograms
  // Global histograms
  hist_vismass      = new TH1F("hist_vismass",     "Visible Mass; M^{vis}_{ll};Events", 10, 0, 200);
  hist_etmiss       = new TH1F("hist_etmiss",      "Missing Transverse Momentum;E_{T}^{Miss} [GeV];Events", 20, 0, 200);
  hist_vxp_z        = new TH1F("hist_vxp_z",       "Primary Vertex Position; z_{Vertex}; Events", 40, -200, 200);
  hist_pvxp_n       = new TH1F("hist_pvxp_n",      "Number of Vertices; N_{vertex}; Events", 30, -0.5, 29.5);
  hist_mt           = new TH1F("hist_mt",          "Transverse Mass; M^{W}_{T};Events", 40, 0, 200);
  hist_Topmass      = new TH1F("hist_Topmass","Top Invariant Mass; M_{Top};Events", 40, 0, 300);
  hist_Wmass        = new TH1F("hist_Wmass","W Invariant Mass; M_{W};Events", 40, 0, 200);
  hist_mazz         = new TH1F("hist_mazz","Z' Invariant Mass; M_{Z'};Events", 40, 0, 3000);

  // Leading Lepton histograms
  hist_leadleptpt   = new TH1F("hist_leadleptpt",  "Leading Lepton Transverse Momentum;p_{T}^{leadlep} [GeV];Leptons", 20, 0, 200);
  hist_leadlepteta  = new TH1F("hist_leadlepteta", "Leading Lepton Pseudorapidity; #eta^{leadlep}; Leptons", 30, -3, 3);
  hist_leadleptE    = new TH1F("hist_leadleptE",   "Leading Lepton Energy; E^{leadlep} [GeV]; Leptons", 30, 0, 300);
  hist_leadleptphi  = new TH1F("hist_leadleptphi", "Leading Lepton Azimuthal Angle ; #phi^{leadlep}; Leptons", 32, -3.2, 3.2);
  hist_leadleptch   = new TH1F("hist_leadleptch",  "Leading Lepton Charge; Q^{leadlep}; Leptons", 7, -1.75, 1.75);
  hist_leadleptID   = new TH1F("hist_leadleptID",  "Leading Lepton Absolute PDG ID; |PDG ID|^{leadlep}; Leptons",  31, -0.5, 30.5);
  hist_leadlept_ptc = new TH1F("hist_leadlept_ptc", "Leading Lepton Relative Transverse Momentum Isolation; ptconerel30^{leadlep}; Leptons", 40, -0.1, 0.4);
  hist_leadleptetc  = new TH1F("hist_leadleptetc", "Leading Lepton Relative Transverse Energy Isolation; etconerel20^{leadlep}; Leptons", 40, -0.1, 0.4);
  hist_leadlepz0    = new TH1F("hist_leadlepz0",   "Leading Lepton z0 impact parameter; z_{0}^{leadlep} [mm]; Leptons", 40, -1, 1);
  hist_leadlepd0    = new TH1F("hist_leadlepd0",   "Leading Lepton d0 impact parameter; d_{0}^{leadlep} [mm]; Leptons", 40, -1, 1);
    
  // Jet variables histograms
  hist_n_jets           = new TH1F("hist_n_jets",      "Number of Jets;N_{jets};Events", 10, -0.5, 9.5);
  hist_leadjet_pt       = new TH1F("hist_jet_pt",      "Jet Transverse Momentum;p_{T}^{jet} [GeV];Jets", 20, 0, 200);
  hist_leadjet_m        = new TH1F("hist_jet_m",       "Jet Mass; m^{jet} [MeV]; Jets", 20, 0, 20);
  hist_leadjet_jvf      = new TH1F("hist_jet_jvf",     "Jet Vertex Fraction; JVF ; Jets", 20, 0, 1);
  hist_leadjet_eta      = new TH1F("hist_jet_eta",     "Jet Pseudorapidity; #eta^{jet}; Jets", 30, -3, 3);
  hist_leadjet_MV1      = new TH1F("hist_jet_MV1",     "Jet MV1; MV1 weight ; Jets", 20, 0, 1);    
    
}

////////////////////////////////////////////////////////////////////////////////
void ZPrimeResolvedAnalysis::FillOutputList()
{
  // histograms
   
  // Add Global histograms
  GetOutputList()->Add(hist_etmiss);
  GetOutputList()->Add(hist_vxp_z);
  GetOutputList()->Add(hist_pvxp_n);
  GetOutputList()->Add(hist_vismass);
  GetOutputList()->Add(hist_mt);    
  GetOutputList()->Add(hist_Wmass); 
  GetOutputList()->Add(hist_Topmass);
  GetOutputList()->Add(hist_mazz);

  // Add Leading Lepton histograms
  GetOutputList()->Add(hist_leadleptpt);
  GetOutputList()->Add(hist_leadlepteta);
  GetOutputList()->Add(hist_leadleptE);
  GetOutputList()->Add(hist_leadleptphi);
  GetOutputList()->Add(hist_leadleptch);
  GetOutputList()->Add(hist_leadleptID);
  GetOutputList()->Add(hist_leadlept_ptc);
  GetOutputList()->Add(hist_leadleptetc);
  GetOutputList()->Add(hist_leadlepz0);
  GetOutputList()->Add(hist_leadlepd0);
    
  // Add Jet variables histograms
  GetOutputList()->Add(hist_n_jets);
  GetOutputList()->Add(hist_leadjet_pt);
  GetOutputList()->Add(hist_leadjet_m);
  GetOutputList()->Add(hist_leadjet_jvf);
  GetOutputList()->Add(hist_leadjet_eta);
  GetOutputList()->Add(hist_leadjet_MV1);
    
}

////////////////////////////////////////////////////////////////////////////////
void ZPrimeResolvedAnalysis::WriteHistograms()
{
  // histograms

  // Write Global histograms
  hist_etmiss->Write();
  hist_vxp_z->Write();
  hist_pvxp_n->Write();
  hist_vismass->Write();
  hist_mt->Write();
  hist_Wmass->Write();
  hist_Topmass->Write();
  hist_mazz->Write();
    
  // Write Leading Lepton histograms
  hist_leadleptpt->Write();
  hist_leadlepteta->Write();
  hist_leadleptE->Write();
  hist_leadleptphi->Write();
  hist_leadleptch->Write();
  hist_leadleptID->Write();
  hist_leadlept_ptc->Write();
  hist_leadleptetc->Write();
  hist_leadlepz0->Write();
  hist_leadlepd0->Write();
    
  // Write Jet variables histograms
  hist_n_jets->Write();
  hist_leadjet_pt->Write();
  hist_leadjet_m->Write();
  hist_leadjet_jvf->Write();
  hist_leadjet_eta->Write();
  hist_leadjet_MV1->Write();    

}

void ZPrimeResolvedAnalysis::FillHistogramsGlobal( double m, float w , TString s)
{

  //Fill Global histograms
  if (s.Contains("hist_vismass")) hist_vismass->Fill(m,w);
  if (s.Contains("hist_etmiss")) hist_etmiss->Fill(m,w);
  if (s.Contains("hist_vxp_z")) hist_vxp_z->Fill(m,w);
  if (s.Contains("hist_pvxp_n")) hist_pvxp_n->Fill(m,w);
  if (s.Contains("hist_mt")) hist_mt->Fill(m,w);
  if (s.Contains("hist_Topmass")) hist_Topmass->Fill(m,w);
  if (s.Contains("hist_Wmass")) hist_Wmass->Fill(m,w);
  if (s.Contains("hist_mazz")) hist_mazz->Fill(m,w);

}

void ZPrimeResolvedAnalysis::FillHistogramsLeadlept( double m, float w , TString s)
{
  //Leading lepton histograms
  if (s.Contains("hist_leadleptpt")) hist_leadleptpt->Fill(m,w);
  if (s.Contains("hist_leadlepteta")) hist_leadlepteta->Fill(m,w);
  if (s.Contains("hist_leadleptE")) hist_leadleptE->Fill(m,w);
  if (s.Contains("hist_leadleptphi")) hist_leadleptphi->Fill(m,w);
  if (s.Contains("hist_leadleptch")) hist_leadleptch->Fill(m,w);
  if (s.Contains("hist_leadleptID")) hist_leadleptID->Fill(m,w);
  if (s.Contains("hist_leadlept_ptc")) hist_leadlept_ptc->Fill(m,w);
  if (s.Contains("hist_leadleptetc")) hist_leadleptetc->Fill(m,w);
  if (s.Contains("hist_leadlepz0")) hist_leadlepz0->Fill(m,w);
  if (s.Contains("hist_leadlepd0")) hist_leadlepd0->Fill(m,w);    
}


void ZPrimeResolvedAnalysis::FillHistogramsLeadJet( double m, float w , TString s)
{
  //Leading Jet histograms
  if (s.Contains("hist_n_jets")) hist_n_jets->Fill(m,w);
  if (s.Contains("hist_leadjet_pt")) hist_leadjet_pt->Fill(m,w);
  if (s.Contains("hist_leadjet_m")) hist_leadjet_m->Fill(m,w);
  if (s.Contains("hist_leadjet_jvf")) hist_leadjet_jvf->Fill(m,w);
  if (s.Contains("hist_leadjet_eta")) hist_leadjet_eta->Fill(m,w);
  if (s.Contains("hist_leadjet_MV1")) hist_leadjet_MV1->Fill(m,w);
}


////////////////////////////////////////////////////////////////////////////////
