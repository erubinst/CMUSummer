#include "AnaExample/nTupleAnalysis/interface/eventData.h"

using namespace nTupleAnalysis;


eventData::eventData(TChain* t, bool mc, std::string y, bool d){
  std::cout << "eventData::eventData()" << std::endl;
  tree  = t;
  isMC  = mc;
  year  = y;
  debug = d;

  //std::cout << "eventData::eventData() tree->Lookup(true)" << std::endl;
  //tree->Lookup(true);
  std::cout << "eventData::eventData() tree->LoadTree(0)" << std::endl;
  tree->LoadTree(0);
  initBranch(tree, "run",             run);
  initBranch(tree, "luminosityBlock", lumiBlock);
  initBranch(tree, "event",           event);
  initBranch(tree, "PV_npvs",         nPVs);
  initBranch(tree, "PV_npvsGood",     nPVsGood);
  if(isMC){
    initBranch(tree, "genWeight", genWeight);
  }

  ////triggers
  ////trigObjs = new trigData("TrigObj", tree);
  //if(year=="2016"){
  //  initBranch(tree, "HLT_QuadJet45_TripleBTagCSV_p087",            HLT_4j45_3b087);
  //  initBranch(tree, "HLT_DoubleJet90_Double30_TripleBTagCSV_p087", HLT_2j90_2j30_3b087);
  //}
  //if(year=="2018"){
  //  initBranch(tree, "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5", HLT_HT330_4j_75_60_45_40_3b);
  //  initBranch(tree, "HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1",    HLT_4j_103_88_75_15_2b_VBF1);
  //  initBranch(tree, "HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2",              HLT_4j_103_88_75_15_1b_VBF2);
  //  initBranch(tree, "HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71",       HLT_2j116_dEta1p6_2b);
  //  initBranch(tree, "HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02",            HLT_J330_m30_2b);
  //  initBranch(tree, "HLT_PFJet500",            HLT_j500);
  //  initBranch(tree, "HLT_DiPFJetAve300_HFJEC", HLT_2j300ave);
  //  //                            HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v
  //  //                            HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v
  //  //                            HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_v
  //  //                            HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_v
  //  //                            HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_v
  //  // HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v
  //  // HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v
  //  // HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_v
  //}

  std::cout << "eventData::eventData() Initialize jets and muons" << std::endl;
  treeJets  = new jetData( "Jet",  tree);
  treeMuons = new muonData("Muon", tree);
} 


void eventData::update(int e){
  //if(e>2546040) debug = true;
  if(debug){
    std::cout<<"Get Entry "<<e<<std::endl;
    std::cout<<tree->GetCurrentFile()->GetName()<<std::endl;
    tree->Show(e);
  }
  Long64_t loadStatus = tree->LoadTree(e);
  if(loadStatus<0){
    std::cout << "Error "<<loadStatus<<" getting event "<<e<<std::endl; 
    return;
  }
  tree->GetEntry(e);
  if(debug) std::cout<<"Got Entry "<<e<<std::endl;

  if(debug) std::cout<<"Reset eventData"<<std::endl;

  ////Trigger
  //if(year=="2016"){
  //  passHLT = HLT_4j45_3b087 || HLT_2j90_2j30_3b087;
  //}
  //if(year=="2018"){
  //  passHLT = HLT_HT330_4j_75_60_45_40_3b || HLT_4j_103_88_75_15_2b_VBF1 || HLT_4j_103_88_75_15_1b_VBF2 || HLT_2j90_2j30_3b087 || HLT_J330_m30_2b || HLT_j500 || HLT_2j300ave;
  //}

  passHLT = true;

  //Objects
  if(debug) std::cout << "Get Jets\n";
  allJets = treeJets->getJets(20);
  selJets = treeJets->getJets(allJets, jetPtMin, 1e6, jetEtaMax, doJetCleaning);
  
  if(debug) std::cout << "Get Muons\n";
  allMuons = treeMuons->getMuons();
  isoMuons = treeMuons->getMuons(40, 2.4, 2, true);

  if(debug) std::cout<<"eventData updated\n";
  return;
}




void eventData::dump(){

  std::cout << "   Run: " << run    << std::endl;
  std::cout << " Event: " << event  << std::endl;  
  std::cout << "Weight: " << weight << std::endl;
  std::cout << " allJets: " << allJets .size() << " |  selJets: " << selJets .size() << std::endl;
  std::cout << "allMuons: " << allMuons.size() << " | isoMuons: " << isoMuons.size() << std::endl;

  return;
}

eventData::~eventData(){} 

