#ifndef MICROPHYSICS_HPP
#define MICROPHYSICS_cPP

// Athena headers
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include "../athena.hpp"

class MeshBlock;
class ParameterInput;
template<typename T> class AthenaArray;

class Microphysics {
public:
  // data
  int ncycle;

  // functions
  Microphysics(MeshBlock *pmb, ParameterInput *pin);
  ~Microphysics();
  Real GetRd() const {return Rd_;}
  Real GetCvRatio(int n) const {return rcv_[n];}
  Real GetLatent(int n) const {return latent_[n];}
  Real GetMassRatio(int n) const {return eps_[n];}
  void SaturationAdjust(AthenaArray<Real> &w, AthenaArray<Real> &u);

  Real Temperature(Real prim[NHYDRO]) const {
    Real qa = 1.;
    for (int n = ICD; n < ICD + NVAPOR; ++n)
      qa -= prim[n];
    Real qb = 1.;
    for (int n = 1; n < ITR; ++n)
      qb += prim[n]*(eps_[n] - 1.);
    return (prim[IPR]*qb)/(prim[IDN]*Rd_*qa);
  }

  // Volume heat capacity
  Real HeatCapacity(Real prim[NHYDRO]) const {
    Real qa = 1., qb = 1.;
    Real gamma = pmy_block_->peos->GetGamma();
    for (int n = 1; n < ITR; ++n) {
      qa += prim[n]*(rcv_[n] - 1.);
      qb += prim[n]*(eps_[n] - 1.);
    }
    return Rd_*qa/qb/(gamma - 1.);
  }

private:
  MeshBlock *pmy_block_;
  Real Rd_;
  Real eps_[1+2*NVAPOR];
  Real rcv_[1+2*NVAPOR];
  Real latent_[1+2*NVAPOR];
};

#endif