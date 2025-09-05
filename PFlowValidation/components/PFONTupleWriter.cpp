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
  std::vector<Gaudi::details::BranchWrapper> m_branchWrappers;
  mutable std::map<std::string, std::pair<std::string, std::any>> m_fields;
  // std::map<const char*, std::any> m_fields;
  mutable bool m_tree_initialised = false;

  StatusCode initialize() override {
    m_file = std::make_unique<TFile>(m_filePath.value().c_str(), "RECREATE");
    m_file->cd();
    m_tree = std::make_unique<TTree>("PFONtuple", "PFONtuple");

    // m_branches.push_back(m_tree->Branch("test_branch", (double*)nullptr));

    return Consumer::initialize();
  }

  void createBranches() {
    m_branchWrappers.reserve(m_fields.size());
    for (const auto& [name, value]: m_fields) {
      const auto& typeName = value.first;
      m_branchWrappers.emplace_back(m_tree, typeName, name, "", this->name());
    }
    m_tree_initialised = true;
  }

  void fillTree() {
    // TODO: grab mutex
    if (!m_tree_initialised) {
      createBranches();
    }
    size_t idx = 0;
    for (const auto& [name, value] : m_fields) {
      // Access the corresponding branch wrapper
      if (idx < m_branchWrappers.size()) {
        auto& branchWrapper = m_branchWrappers[idx];
        // Set the address or value for the branch here, e.g.:
        if (value.second.type() == typeid(double)) {
          double val = std::any_cast<double>(value.second);

        } else if (value.second.type() == typeid(int)) {

        } else {
          // never happens
        }
        // branchWrapper.setBranchData(std::any_cast<value.second.type()>(value.second));
      }
      ++idx;
    }
    m_tree->Fill();
  }

  void operator()(const edm4hep::ReconstructedParticleCollection& PFOs, const edm4hep::MCParticleCollection& MCParticles,
                  const edm4hep::RecoMCParticleLinkCollection& RecoMCTruthLinks,
                  const edm4hep::RecoMCParticleLinkCollection& MCTruthRecoLinks) const override {
    // double foo = 42.;
    // m_branches[0]->SetAddress(&foo);
    // m_tree->Fill();

    // const auto linkNavigator = podio::LinkNavigator(recoMcLinks);
    // // For podio::LinkCollections with disparate types just use getLinked
    // const auto linkedRecs = linkNavigator.getLinked(mcParticle);


    // make navigators for both directions
    podio::LinkNavigator RecoMCTruthLinkNavigator(RecoMCTruthLinks);
    podio::LinkNavigator MCTruthRecoLinkNavigator(MCTruthRecoLinks);
    // Iterate over all mcparticles
    for (const auto& mcp : MCParticles) {
      // skip the ones not simulated by ddsim
      // i.e. only keep gen status one or created in simulation
      if (mcp.getGeneratorStatus() != 1 && !mcp.isCreatedInSimulation()) {
        continue;
      }
      m_fields["foo"] = std::make_pair("double", 42.);
      // check a bunch of the statuses and also store them in ntuple
      // get PDG and store it in ntuple
      // need to figure out if there are any surprise particles
      // store charge in ntuple
      // store mass, 4 momentum in ntuple
      // figure out how many tracker hits, ecal hits, hcal hits, muon hits this particle caused and store in ntuple
      // somehow also figure out if all of those are really from this particle or if there was a non-stored secondary
      // get related PFOs, how many are related by track, how many by cluster, are they the same, is the strongest the same
      // if the strongest is not the same which one is stronger? (closest in energy?)
      // continue with the strongest PFO
      // how pure is the pfo? how many foreign tracks are in the track, how many in the cluster, what are the fractions
      // what is the PDG assigned by pandora
      // store 4 momentum in ntuple

    }

  }

  StatusCode finalize() override {
    m_tree->Write();
    return Consumer::finalize();
  }
};

DECLARE_COMPONENT(PFONTupleWriter)
