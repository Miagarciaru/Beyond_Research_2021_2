//////////////////////////////////////////////////////////
#ifndef ZPrimeResolvedAnalysis_h
#define ZPrimeResolvedAnalysis_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TH1.h"
#include "vector"

class ZPrimeResolvedAnalysis : public TSelector {
  public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain

  //////////////////////////////////////////////////////////
  // histograms

  // Global variables histograms
  TH1F *hist_etmiss       = 0;
  TH1F *hist_vxp_z        = 0;
  TH1F *hist_pvxp_n       = 0;
  TH1F *hist_vismass      = 0;
  TH1F *hist_mt           = 0;    
  TH1F *hist_Wmass        = 0;
  TH1F *hist_Topmass      = 0;
  TH1F *hist_mazz         = 0;
    
  // Leading Lepton histograms
  TH1F *hist_leadleptpt   = 0;
  TH1F *hist_leadlepteta  = 0;
  TH1F *hist_leadleptE    = 0;
  TH1F *hist_leadleptphi  = 0;
  TH1F *hist_leadleptch   = 0;
  TH1F *hist_leadleptID   = 0;
  TH1F *hist_leadlept_ptc = 0;
  TH1F *hist_leadleptetc  = 0;
  TH1F *hist_leadlepz0    = 0;
  TH1F *hist_leadlepd0    = 0;

  // Jet variables histograms
  TH1F *hist_n_jets       = 0;
  TH1F *hist_leadjet_pt   = 0;
  TH1F *hist_leadjet_m    = 0;
  TH1F *hist_leadjet_jvf  = 0;
  TH1F *hist_leadjet_eta  = 0;
  TH1F *hist_leadjet_MV1  = 0;

  //////////////////////////////////////////////////////////
  // Declaration of leaf types
  Int_t           runNumber;
  Int_t           eventNumber;
  Int_t           channelNumber;
  Float_t         mcWeight;
  Int_t           pvxp_n;
  Float_t         vxp_z;
  Float_t         scaleFactor_PILEUP;
  Float_t         scaleFactor_ELE;
  Float_t         scaleFactor_MUON;
  Float_t         scaleFactor_BTAG;
  Float_t         scaleFactor_TRIGGER;
  Float_t         scaleFactor_JVFSF;
  Float_t         scaleFactor_ZVERTEX;
  Bool_t          trigE;
  Bool_t          trigM;
  Bool_t          passGRL;
  Bool_t          hasGoodVertex;
  UInt_t          lep_n;
  vector<bool>    *lep_truthMatched;
  vector<bool>    *lep_trigMatched;  
  vector<float>   *lep_pt;
  vector<float>   *lep_eta;  
  vector<float>   *lep_phi;  
  vector<float>   *lep_E;  
  vector<float>   *lep_z0;  
  vector<float>   *lep_charge;  //Int?
  vector<unsigned int>   *lep_type;  
  vector<unsigned int>   *lep_flag;  
  vector<float>   *lep_ptcone30; 
  vector<float>   *lep_etcone20;   
  vector<float>   *lep_trackd0pvunbiased;   
  vector<float>   *lep_tracksigd0pvunbiased; 
    
  Float_t         met_et;
  Float_t         met_phi;
  UInt_t          jet_n;
  UInt_t          alljet_n;
  vector<float>   *jet_pt;  
  vector<float>   *jet_eta;  
  vector<float>   *jet_phi;   
  vector<float>   *jet_E; 
  vector<float>   *jet_m; 
  vector<float>   *jet_jvf;  
  vector<int>     *jet_trueflav; 
  vector<bool>    *jet_truthMatched;   //Int
  vector<float>   *jet_SV0;   
  vector<float>   *jet_MV1; 
    
  /////////////////////////////////////////////////////////
    // List of branches
  TBranch        *b_runNumber;   //!
  TBranch        *b_eventNumber;   //!
  TBranch        *b_channelNumber;   //!
  TBranch        *b_mcWeight;   //!
  TBranch        *b_pvxp_n;   //!
  TBranch        *b_vxp_z;   //!  
  TBranch        *b_scaleFactor_PILEUP;   //!
  TBranch        *b_scaleFactor_ELE;   //! 
  TBranch        *b_scaleFactor_MUON;   //!
  TBranch        *b_scaleFactor_BTAG;   //!
  TBranch        *b_scaleFactor_TRIGGER;   //!
  TBranch        *b_scaleFactor_JVFSF;   //!
  TBranch        *b_scaleFactor_ZVERTEX;   //!
  TBranch        *b_trigE;   //!
  TBranch        *b_trigM;   //!  
  TBranch        *b_passGRL;   //!
  TBranch        *b_hasGoodVertex;   //!
  TBranch        *b_lep_n;   //!
  TBranch        *b_lep_truthMatched;   //!
  TBranch        *b_lep_trigMatched;   //!
  TBranch        *b_lep_pt;   //!
  TBranch        *b_lep_eta;   //!
  TBranch        *b_lep_phi;   //!
  TBranch        *b_lep_E;   //!
  TBranch        *b_lep_z0;   //!
  TBranch        *b_lep_charge;   //!    
  TBranch        *b_lep_type;   //!
  TBranch        *b_lep_flag;   //!
  TBranch        *b_lep_ptcone30;   //!
  TBranch        *b_lep_etcone20;   //!
  TBranch        *b_lep_trackd0pvunbiased;   //!
  TBranch        *b_lep_tracksigd0pvunbiased;   //!
  TBranch        *b_met_et;   //!
  TBranch        *b_met_phi;   //!
  TBranch        *b_jet_n;   //!
  TBranch        *b_alljet_n;   //!  
  TBranch        *b_jet_pt;   //!
  TBranch        *b_jet_eta;   //!
  TBranch        *b_jet_phi;   //!
  TBranch        *b_jet_E;   //!
  TBranch        *b_jet_m;   //!
  TBranch        *b_jet_jvf;   //!
  TBranch        *b_jet_trueflav;   //!
  TBranch        *b_jet_truthMatched;   //!
  TBranch        *b_jet_SV0;   //!
  TBranch        *b_jet_MV1;   //!

  ZPrimeResolvedAnalysis(TTree * =0) : fChain(0) { }
  virtual ~ZPrimeResolvedAnalysis() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();

  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual void    FillHistogramsGlobal( double m, float w , TString s);
  virtual void    FillHistogramsLeadlept( double m, float w , TString s);
  virtual void    FillHistogramsLeadJet( double m, float w , TString s);
    
    
  // Get Output List to save our histograms in the output file
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    define_histograms();
  virtual void    FillOutputList();
  virtual void    WriteHistograms();
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  int nEvents;

  ClassDef(ZPrimeResolvedAnalysis,0);
};

#endif

#ifdef ZPrimeResolvedAnalysis_cxx
void ZPrimeResolvedAnalysis::Init(TTree *tree)
{
    /*
   lep_truthMatched = 0;
   lep_trigMatched = 0;
   lep_pt = 0;
   lep_eta = 0;
   lep_phi = 0;
   lep_E = 0;
   lep_z0 = 0;
   lep_charge = 0;
   lep_type = 0;
   lep_isTightID = 0;
   lep_ptcone30 = 0;
   lep_etcone20 = 0;
   lep_trackd0pvunbiased = 0;
   lep_tracksigd0pvunbiased = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_E = 0;
   jet_jvt = 0;
   jet_trueflav = 0;
   jet_truthMatched = 0;
   jet_MV2c10 = 0;
   largeRjet_pt = 0;
   largeRjet_eta = 0;
   largeRjet_phi = 0;
   largeRjet_E = 0;
   largeRjet_m = 0;
   largeRjet_truthMatched = 0;
   largeRjet_D2 = 0;
   largeRjet_tau32 = 0;
   lep_pt_syst = 0;
   jet_pt_syst = 0;
   largeRjet_pt_syst = 0;
    */
    
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
  fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
  fChain->SetBranchAddress("channelNumber", &channelNumber, &b_channelNumber);
  fChain->SetBranchAddress("mcWeight", &mcWeight, &b_mcWeight);
  fChain->SetBranchAddress("pvxp_n", &pvxp_n, &b_pvxp_n);
  fChain->SetBranchAddress("vxp_z", &vxp_z, &b_vxp_z);
  fChain->SetBranchAddress("scaleFactor_PILEUP", &scaleFactor_PILEUP, &b_scaleFactor_PILEUP);
  fChain->SetBranchAddress("scaleFactor_ELE", &scaleFactor_ELE, &b_scaleFactor_ELE);
  fChain->SetBranchAddress("scaleFactor_MUON", &scaleFactor_MUON, &b_scaleFactor_MUON);
  fChain->SetBranchAddress("scaleFactor_BTAG", &scaleFactor_BTAG, &b_scaleFactor_BTAG);
  fChain->SetBranchAddress("scaleFactor_TRIGGER", &scaleFactor_TRIGGER, &b_scaleFactor_TRIGGER);
  fChain->SetBranchAddress("scaleFactor_JVFSF", &scaleFactor_JVFSF, &b_scaleFactor_JVFSF);
  fChain->SetBranchAddress("scaleFactor_ZVERTEX", &scaleFactor_ZVERTEX, &b_scaleFactor_ZVERTEX);
  fChain->SetBranchAddress("trigE", &trigE, &b_trigE);
  fChain->SetBranchAddress("trigM", &trigM, &b_trigM);
  fChain->SetBranchAddress("passGRL", &passGRL, &b_passGRL);
  fChain->SetBranchAddress("hasGoodVertex", &hasGoodVertex, &b_hasGoodVertex);
  fChain->SetBranchAddress("lep_n", &lep_n, &b_lep_n);
  fChain->SetBranchAddress("lep_truthMatched", lep_truthMatched, &b_lep_truthMatched);
  fChain->SetBranchAddress("lep_trigMatched", lep_trigMatched, &b_lep_trigMatched);
  fChain->SetBranchAddress("lep_pt", lep_pt, &b_lep_pt);
  fChain->SetBranchAddress("lep_eta", lep_eta, &b_lep_eta);
  fChain->SetBranchAddress("lep_phi", lep_phi, &b_lep_phi);
  fChain->SetBranchAddress("lep_E", lep_E, &b_lep_E);
  fChain->SetBranchAddress("lep_z0", lep_z0, &b_lep_z0);
  fChain->SetBranchAddress("lep_charge", lep_charge, &b_lep_charge);
  fChain->SetBranchAddress("lep_type", lep_type, &b_lep_type);
  fChain->SetBranchAddress("lep_flag", lep_flag, &b_lep_flag);
  fChain->SetBranchAddress("lep_ptcone30", lep_ptcone30, &b_lep_ptcone30);
  fChain->SetBranchAddress("lep_etcone20", lep_etcone20, &b_lep_etcone20);
  fChain->SetBranchAddress("lep_trackd0pvunbiased", lep_trackd0pvunbiased, &b_lep_trackd0pvunbiased);
  fChain->SetBranchAddress("lep_tracksigd0pvunbiased", lep_tracksigd0pvunbiased, &b_lep_tracksigd0pvunbiased);
  fChain->SetBranchAddress("met_et", &met_et, &b_met_et);
  fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
  fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
  fChain->SetBranchAddress("alljet_n", &alljet_n, &b_alljet_n);
  fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
  fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
  fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
  fChain->SetBranchAddress("jet_E", jet_E, &b_jet_E);
  fChain->SetBranchAddress("jet_m", jet_m, &b_jet_m);
  fChain->SetBranchAddress("jet_jvf", jet_jvf, &b_jet_jvf);
  fChain->SetBranchAddress("jet_trueflav", jet_trueflav, &b_jet_trueflav);
  fChain->SetBranchAddress("jet_truthMatched", jet_truthMatched, &b_jet_truthMatched);
  fChain->SetBranchAddress("jet_SV0", jet_SV0, &b_jet_SV0);
  fChain->SetBranchAddress("jet_MV1", jet_MV1, &b_jet_MV1);

}

Bool_t ZPrimeResolvedAnalysis::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef ZPrimeResolvedAnalysis_cxx
