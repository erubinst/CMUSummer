// -*- C++ -*-
#if !defined(fourVectorHists_H)
#define fourVectorHists_H

#include <iostream>
#include <TH1F.h>
#include <TLorentzVector.h>
#include "PhysicsTools/FWLite/interface/TFileService.h"

namespace nTupleAnalysis {

  class fourVectorHists {
  public:
    
    TH1F* pt_s;
    TH1F* pt_m;
    TH1F* pt_l;

    TH1F* eta;
    TH1F* phi;

    TH1F* m;
    TH1F* m_l;
    TH1F* e;

    fourVectorHists(std::string,       TFileDirectory&, std::string title = "");
    fourVectorHists(std::string, fwlite::TFileService&, std::string title = "");

    void makeHists(std::string, TFileDirectory&, std::string title);
    void Fill(TLorentzVector&, float);

    ~fourVectorHists(); 

  };

}
#endif // fourVectorHists_H
