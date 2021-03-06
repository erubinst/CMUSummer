#include "TChain.h"

#include "nTupleAnalysis/baseClasses/interface/jetData.h"

using namespace nTupleAnalysis;

//jet object
jet::jet(){}
jet::jet(UInt_t i, jetData* data){

  cleanmask = data->cleanmask[i];

  pt  = data->pt [i];
  eta = data->eta[i];
  phi = data->phi[i];
  m   = data->m  [i];
  p = TLorentzVector();
  p.SetPtEtaPhiM(pt, eta, phi, m);
  e = p.E();

  bRegCorr = data->bRegCorr[i];

  deepB     = data->deepB[i];
  CSVv2     = data->CSVv2[i];
  deepFlavB = data->deepFlavB[i];

}

jet::jet(TLorentzVector& vec, float tag){

  p = TLorentzVector(vec);
  pt  = p.Pt();
  eta = p.Eta();
  phi = p.Phi();
  m   = p.M();
  e   = p.E();

  bRegCorr = pt;

  deepB = tag;
  CSVv2 = tag;
  deepFlavB = tag;

}

void jet::bRegression(){
  p  *= bRegCorr;
  pt  = p.Pt();
  eta = p.Eta();
  phi = p.Phi();
  m   = p.M();
  e   = p.E();
}

jet::~jet(){}


//access tree
jetData::jetData(std::string name, TChain* tree){

  initBranch(tree, ("n"+name).c_str(), n );

  initBranch(tree, (name+"_cleanmask").c_str(), cleanmask );

  initBranch(tree, (name+"_pt"  ).c_str(), pt  );  
  initBranch(tree, (name+"_eta" ).c_str(), eta );  
  initBranch(tree, (name+"_phi" ).c_str(), phi );  
  initBranch(tree, (name+"_mass").c_str(), m   );  

  initBranch(tree, (name+"_bRegCorr").c_str(), bRegCorr );  

  initBranch(tree, (name+"_btagDeepB"    ).c_str(), deepB     );
  initBranch(tree, (name+"_btagCSVV2"    ).c_str(), CSVv2     );
  initBranch(tree, (name+"_btagDeepFlavB").c_str(), deepFlavB );
  //initBranch(tree, (name+"_").c_str(),  );

}

std::vector< std::shared_ptr<jet> > jetData::getJets(float ptMin, float ptMax, float etaMax, bool clean, float tagMin, std::string tagger, bool antiTag){
  
  std::vector< std::shared_ptr<jet> > outputJets;
  float *tag = CSVv2;
  if(tagger == "deepB")     tag = deepB;
  if(tagger == "deepFlavB") tag = deepFlavB;

  for(UInt_t i = 0; i < n; ++i){
    if(clean && cleanmask[i] == 0) continue;
    if(          pt[i]  <  ptMin ) continue;
    if(          pt[i]  >= ptMax ) continue;
    if(    fabs(eta[i]) > etaMax ) continue;
    if(antiTag^(tag[i]  < tagMin)) continue; // antiTag XOR (jet fails tagMin). This way antiTag inverts the tag criteria to select untagged jets
    outputJets.push_back(std::make_shared<jet>(jet(i, this)));
  }

  return outputJets;
}

std::vector< std::shared_ptr<jet> > jetData::getJets(std::vector< std::shared_ptr<jet> > inputJets, float ptMin, float ptMax, float etaMax, bool clean, float tagMin, std::string tagger, bool antiTag){
  
  std::vector< std::shared_ptr<jet> > outputJets;

  for(auto &jet: inputJets){
    if(clean && jet->cleanmask == 0) continue;
    if(         jet->pt   <  ptMin ) continue;
    if(         jet->pt   >= ptMax ) continue;
    if(    fabs(jet->eta) > etaMax ) continue;

    if(     tagger == "deepFlavB" && antiTag^(jet->deepFlavB < tagMin)) continue;
    else if(tagger == "deepB"     && antiTag^(jet->deepB     < tagMin)) continue;
    else if(tagger == "CSVv2"     && antiTag^(jet->CSVv2     < tagMin)) continue;
    outputJets.push_back(jet);
  }

  return outputJets;
}


jetData::~jetData(){}

