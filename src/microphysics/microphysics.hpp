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
  Real GetTerminalVelocity() const {return termv_;}

  Real Temperature(Real prim[NHYDRO]) const {
    Real qa = 1.;
    for (int n = ICD; n < ICD + NVAPOR; ++n)
      qa -= prim[n];
    Real qb = 1.;
    for (int n = 1; n < ITR; ++n)
      qb += prim[n]*(eps_[n] - 1.);
    return (prim[IPR]*qb)/(prim[IDN]*Rd_*qa);
  }

  // heat capacity
  Real Cp(Real prim[NHYDRO]) const {
    Real qa = 1., qb = 1.;
    Real gamma = pmy_block_->peos->GetGamma();
    for (int n = 1; n < ITR; ++n) {
      qa += prim[n]*(rcv_[n] - 1.);
      qb += prim[n]*(eps_[n] - 1.);
    }
    return gamma/(gamma - 1.)*Rd_*qa/qb;
  }

  void SaturationAdjustment(AthenaArray<Real> &w,
    int is, int ie, int js, int je, int ks, int ke);
  void Precipitation(AthenaArray<Real> &w, AthenaArray<Real> const &u, Real dt,
    int is, int ie, int js, int je, int ks, int ke);
  void Evaporation(AthenaArray<Real> &w, AthenaArray<Real> const &u, Real dt);

private:
  MeshBlock *pmy_block_;
  Real Rd_;                   // gas constant of dry air
  Real eps_[1+2*NVAPOR];
  Real rcv_[1+2*NVAPOR];
  Real latent_[1+2*NVAPOR];

  Real termv_;                // terminal velocity
  Real autoc_;                // auto-conversion time
  Real evapr_;                // evaporation rate
  Real tiny_number_;          // very small number
};

#endif
