/*
 * PFOtoMCviaClusterLink: PFOs + Reco<->MC links, with 'selected' histograms
 * and runtime toggles for selection + debug counters.
 *
 * Selection (CLD-style, arXiv:1911.12230):
 *  - SelectByType:     reco/MC same coarse type (default true)
 *  - SelectByAngles:   |dphi| < MaxDeltaPhi && |dtheta| < MaxDeltaTheta  (defaults: 2 mrad / 1 mrad)
 *  - SelectPtForCharged: |pT(reco)/pT(true)-1| < MaxRelDpT (default 5%) for charged MC
 *
 * Debug:
 *  - DebugCounters: print rolling and final summaries
 *
 * Histograms (under HistPath, default "/PFO/PFOCluster"):
 *  Inclusive:  mcE, pfoE, mcCos, pfoCos, respE, effE_den/num, effCos_den/num, pidConf, respR
 *  Selected:   mcE_sel, pfoE_sel, mcCos_sel, pfoCos_sel, respE_sel,
 *              effE_sel_den/num, effCos_sel_den/num, respR_sel,
 *              dE_sel [GeV], dpT_sel [GeV] (charged MC only), dTheta_sel [rad], dPhi_sel [rad, signed]
 */

#include "k4FWCore/Consumer.h"

#include "Gaudi/Property.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/SmartIF.h"

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/RecoMCParticleLinkCollection.h"
#include "podio/ObjectID.h"

#include "TFile.h"
#include "TTree.h"

#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <cstdint>

struct PFONTupleWriter final
  : k4FWCore::Consumer<void(const edm4hep::ReconstructedParticleCollection&,
                            const edm4hep::MCParticleCollection&,
                            const edm4hep::RecoMCParticleLinkCollection&,
                            const edm4hep::RecoMCParticleLinkCollection&)> {

  PFONTupleWriter(const std::string& name, ISvcLocator* svcLoc)
    : Consumer(name, svcLoc, {
        KeyValues("InputPFOs",        {"PandoraPFOs"}),
        KeyValues("InputMCParticles", {"MCParticles"}),
        KeyValues("InputRecoMCLink",  {"RecoMCTruthLink"}),
        KeyValues("InputMCRecoLink",  {"MCTruthRecoLink"})
      }) {}


  Gaudi::Property<std::string> m_filePath{this, "FilePath", "test.root", "long blabla"};
  std::unique_ptr<TFile> m_file = nullptr;
  std::unique_ptr<TTree> m_tree = nullptr;
  std::vector<TBranch*> m_branches;

  StatusCode initialize() override {
    m_file = std::make_unique<TFile>(m_filePath.value().c_str(), "RECREATE");
    m_file->cd();
    m_tree = std::make_unique<TTree>("PFONtuple", "PFONtuple");
    
    m_branches.push_back(m_tree->Branch("test_branch", (double*)nullptr));

    return Consumer::initialize();
  }


  void operator()(const edm4hep::ReconstructedParticleCollection& pfos,
                  const edm4hep::MCParticleCollection& mcps,
                  const edm4hep::RecoMCParticleLinkCollection& RecoMCTruthLinks,
                  const edm4hep::RecoMCParticleLinkCollection& MCTruthRecoLinks) const override {
                    double foo = 42.;
                    m_branches[0]->SetAddress(&foo);
                    m_tree->Fill();
  }

  StatusCode finalize() override {
    m_tree->Write();
    return Consumer::finalize();
  }
};

DECLARE_COMPONENT(PFONTupleWriter)
