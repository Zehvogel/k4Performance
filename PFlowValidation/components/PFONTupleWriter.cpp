#include "k4FWCore/Consumer.h"

#include "Gaudi/Property.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/SmartIF.h"
#include "Gaudi/details/BranchWrapper.h"

#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/RecoMCParticleLinkCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "podio/ObjectID.h"
#include "podio/LinkNavigator.h"

#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

struct PFONTupleWriter final
    : k4FWCore::Consumer<void(const edm4hep::ReconstructedParticleCollection&, const edm4hep::MCParticleCollection&,
                              const edm4hep::RecoMCParticleLinkCollection&,
                              const edm4hep::RecoMCParticleLinkCollection&)> {

  PFONTupleWriter(const std::string& name, ISvcLocator* svcLoc)
      : Consumer(name, svcLoc,
                 {KeyValues("InputPFOs", {"PandoraPFOs"}), KeyValues("InputMCParticles", {"MCParticles"}),
                  KeyValues("InputRecoMCLink", {"RecoMCTruthLink"}),
                  KeyValues("InputMCRecoLink", {"MCTruthRecoLink"})}) {}

  Gaudi::Property<std::string> m_filePath{this, "FilePath", "test.root", "long blabla"};
  std::unique_ptr<TFile> m_file = nullptr;
  std::unique_ptr<TTree> m_tree = nullptr;
  mutable std::map<std::string, TBranch*> m_branches;

  StatusCode initialize() override {
    m_file = std::make_unique<TFile>(m_filePath.value().c_str(), "RECREATE");
    m_file->cd();
    m_tree = std::make_unique<TTree>("PFONtuple", "PFONtuple");

    m_branches["foo"] = m_tree->Branch("foo", (double*)nullptr);
    // TODO: initialize branches
    m_branches["mc_pdg"] = m_tree->Branch("mc_pdg", (int*)nullptr);
    m_branches["mc_charge"] = m_tree->Branch("mc_charge", (double*)nullptr);
    m_branches["mc_mass"] = m_tree->Branch("mc_mass", (double*)nullptr);
    m_branches["mc_E"] = m_tree->Branch("mc_E", (double*)nullptr);
    m_branches["mc_Px"] = m_tree->Branch("mc_Px", (double*)nullptr);
    m_branches["mc_Py"] = m_tree->Branch("mc_Py", (double*)nullptr);
    m_branches["mc_Pz"] = m_tree->Branch("mc_Pz", (double*)nullptr);
    m_branches["nPFOs"] = m_tree->Branch("nPFOs", (int*)nullptr);
    m_branches["max_c_weight"] = m_tree->Branch("max_c_weight", (float*)nullptr);
    m_branches["max_t_weight"] = m_tree->Branch("max_t_weight", (float*)nullptr);
    m_branches["has_track"] = m_tree->Branch("has_track", (bool*)nullptr);
    m_branches["has_cluster"] = m_tree->Branch("has_cluster", (bool*)nullptr);
    m_branches["track_equals_cluster"] = m_tree->Branch("track_equals_cluster", (bool*)nullptr);
    m_branches["has_reco_pfo"] = m_tree->Branch("has_reco_pfo", (bool*)nullptr);
    m_branches["pfo_pdg"] = m_tree->Branch("pfo_pdg", (int*)nullptr);
    m_branches["pfo_E"] = m_tree->Branch("pfo_E", (double*)nullptr);
    m_branches["pfo_Px"] = m_tree->Branch("pfo_Px", (double*)nullptr);
    m_branches["pfo_Py"] = m_tree->Branch("pfo_Py", (double*)nullptr);
    m_branches["pfo_Pz"] = m_tree->Branch("pfo_Pz", (double*)nullptr);

    return Consumer::initialize();
  }

  void fillTree(const std::map<std::string, void*>& fields) const {
    // TODO: grab lock here
    for (const auto& [name, ptr] : fields) {
      m_branches[name]->SetAddress(ptr);
    }
    m_tree->Fill();
  }

  /** Return track weight contribution encoded as trackwgt = (int(wgt)%10000)/1000. by
   * the encodeTrackAndClusterWeights function for the
   * ReconstructedParticle-MCParticle (or vise-versa) type of relations.
   * 
   */
  // STOLEN FROM https://github.com/iLCSoft/MarlinUtil/blob/master/source/include/MarlinUtil.h
  static inline float getTrackWeight(float encodedWeight) {
      return float( int(encodedWeight) % 10000 ) / 1000.f;
  }


  /** Return cluster weight contribution encoded as clusterwgt = (int(wgt)/10000)/1000. by
   * the encodeTrackAndClusterWeights function for the
   * ReconstructedParticle-MCParticle (or vise-versa) type of relations.
   */
  // STOLEN FROM https://github.com/iLCSoft/MarlinUtil/blob/master/source/include/MarlinUtil.h
  static inline float getClusterWeight(float encodedWeight) {
      return float( int(encodedWeight) / 10000 ) / 1000.f;
  }

  void operator()(const edm4hep::ReconstructedParticleCollection& PFOs, const edm4hep::MCParticleCollection& MCParticles,
                  const edm4hep::RecoMCParticleLinkCollection& RecoMCTruthLinks,
                  const edm4hep::RecoMCParticleLinkCollection& MCTruthRecoLinks) const override {

    // const auto linkNavigator = podio::LinkNavigator(recoMcLinks);
    // // For podio::LinkCollections with disparate types just use getLinked
    // const auto linkedRecs = linkNavigator.getLinked(mcParticle);

    std::map<std::string, void*> fields;

    // make navigators for both directions
    podio::LinkNavigator RecoMCTruthLinkNavigator(RecoMCTruthLinks);
    podio::LinkNavigator MCTruthRecoLinkNavigator(MCTruthRecoLinks);
    // Iterate over all mcparticles
    for (const auto& mcp : MCParticles) {
      // skip the ones not simulated by ddsim
      // i.e. only keep gen status one
      if (mcp.getGeneratorStatus() != 1) {
        continue;
      }
      double foo = 42.;
      fields["foo"] = &foo;
      //  TODO: check a bunch of the statuses and also store them in ntuple
      // get PDG and store it in ntuple
      int mc_pdg = mcp.getPDG();
      fields["mc_pdg"] = &mc_pdg;
      // need to figure out if there are any surprise particles
      // store charge, mass, 4 momentum in ntuple
      double mc_charge = mcp.getCharge();
      fields["mc_charge"] = &mc_charge;
      double mc_mass = mcp.getMass();
      fields["mc_mass"] = &mc_mass;
      double mc_E = mcp.getEnergy();
      fields["mc_E"] = &mc_E;
      double mc_Px = mcp.getMomentum().x;
      fields["mc_Px"] = &mc_Px;
      double mc_Py = mcp.getMomentum().y;
      fields["mc_Py"] = &mc_Py;
      double mc_Pz = mcp.getMomentum().z;
      fields["mc_Pz"] = &mc_Pz;
      // TODO: figure out how many tracker hits, ecal hits, hcal hits, muon hits this particle caused and store in ntuple
      // somehow also figure out if all of those are really from this particle or if there was a non-stored secondary

      // get related PFOs, how many are related by track, how many by cluster, are they the same, is the strongest the same
      auto linkedPFOs = MCTruthRecoLinkNavigator.getLinked(mcp);
      int nPFOs = linkedPFOs.size();
      fields["nPFOs"] = &nPFOs;
      float max_c_weight = 0.;
      fields["max_c_weight"] = &max_c_weight;
      float max_t_weight = 0.;
      fields["max_t_weight"] = &max_t_weight;
      const edm4hep::ReconstructedParticle* max_c_pfo = nullptr;
      const edm4hep::ReconstructedParticle* max_t_pfo = nullptr;
      for (const auto& [pfo, weight] : linkedPFOs) {
        // need to parse weight
        float cw = getClusterWeight(weight);
        float tw = getTrackWeight(weight);
        if (cw > max_c_weight) {
          max_c_weight = cw;
          max_c_pfo = &pfo;
        }
        if (tw > max_t_weight) {
          max_t_weight = tw;
          max_t_pfo = &pfo;
        }
      }
      // if the strongest is not the same which one is stronger? (closest in energy?)
      // continue with the strongest PFO
      // what to do depends a bit on the initial mc particle
      // TODO: how pure is the pfo? how many foreign tracks are in the track, how many in the cluster, what are the fractions
      const edm4hep::ReconstructedParticle* reco_pfo = nullptr;
      bool has_track = max_t_weight > 0.5;
      fields["has_track"] = &has_track;
      bool has_cluster = max_c_weight > 0.5;
      fields["has_cluster"] = &has_cluster;
      bool track_equals_cluster = max_t_pfo == max_c_pfo;
      fields["track_equals_cluster"] = &track_equals_cluster;
      if (has_cluster) {
        reco_pfo = max_c_pfo;
      }
      if (mc_charge != 0. && has_track) {
        // if it exists use the track instead of the cluster for charged particles
        reco_pfo = max_t_pfo;
      }

      edm4hep::ReconstructedParticle dummy_pfo;
      bool has_reco_pfo = !!reco_pfo;
      fields["has_reco_pfo"] = &has_reco_pfo;
      if (!has_reco_pfo) {
        reco_pfo = &dummy_pfo;
      }

      // what is the PDG assigned by pandora
      int pfo_pdg = reco_pfo->getPDG();
      fields["pfo_pdg"] = &pfo_pdg;
      // store 4 momentum in ntuple
      double pfo_E = reco_pfo->getEnergy();
      fields["pfo_E"] = &pfo_E;
      double pfo_Px = reco_pfo->getMomentum().x;
      fields["pfo_Px"] = &pfo_Px;
      double pfo_Py = reco_pfo->getMomentum().y;
      fields["pfo_Py"] = &pfo_Py;
      double pfo_Pz = reco_pfo->getMomentum().z;
      fields["pfo_Pz"] = &pfo_Pz;

      fillTree(fields);

    }

  }

  StatusCode finalize() override {
    m_tree->Write();
    return Consumer::finalize();
  }
};

DECLARE_COMPONENT(PFONTupleWriter)
