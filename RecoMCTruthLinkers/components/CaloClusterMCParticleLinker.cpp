#include "Gaudi/Property.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/ClusterMCParticleLinkCollection.h"


#include "k4FWCore/Transformer.h"

#include <string>

/**
 * Gaudi functional doing a dummy association between calo clusters and MCParticles based on an angular matching
 */

struct CaloClusterMCParticleLinker final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::ClusterMCParticleLinkCollection>(
          const edm4hep::MCParticleCollection&, const edm4hep::ClusterCollection&)> {
  CaloClusterMCParticleLinker(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {KeyValues("InputMCParticles", {"MCParticles"}),
                          KeyValues("InputCaloClusters", {"CaloClusters"})},
                         {KeyValues("OutputCaloClusterMCParticleLink", {"CaloClusterMCParticleLink"})}) {}

  std::tuple<edm4hep::ClusterMCParticleLinkCollection> operator()(const edm4hep::MCParticleCollection&, const edm4hep::ClusterCollection&) const override {
    info() << "TEST" << endmsg;
  }

  private:

    /// Configurable property to decide whether to produce validation plots
    Gaudi::Property<bool> m_extrapolateToECal{
      this, "ProduceValidationPlots", true, "Decide whether to produce validation plots or not"
    };

    /// Criteria on the angular distance (dR = sqrt(d_theta * d_theta + d_phi * d_phi)) to consider the matching valid
    Gaudi::Property<float> m_cutoff_dR{
      this, "cutoff_dR", 0.03, "Cut off on the angular distance between cluster and MC particles to be associated"
    };


};

DECLARE_COMPONENT(CaloClusterMCParticleLinker)
