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

  // heat capacity
  /*Real Cp(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const {
    Real gamma = pmy_block_->peos->GetGamma();
    Real fsig = 1.;
    for (int n = 1; n < ITR; ++n)
      fsig += w(n,k,j,i)*(rcp_[n] - 1.);
    return gamma/(gamma - 1.)*Rd_*fsig;
  }*/
  Real Cp(Real prim[]) const {
    Real gamma = pmy_block_->peos->GetGamma();
    Real fsig = 1.;
    for (int n = 1; n < ITR; ++n)
      fsig += prim[n]*(rcp_[n] - 1.);
    return gamma/(gamma - 1.)*Rd_*fsig;
  }

  void CalculateTemperature(AthenaArray<Real> const& u);
  void Evaporation(AthenaArray<Real> &u, Real dt);
  void SaturationAdjustment(AthenaArray<Real> &u);
  void Precipitation(AthenaArray<Real> &u, Real dt);

  // temperature
  AthenaArray<Real> T;

private:
  MeshBlock *pmy_block_;
  AthenaArray<bool> recondense_;

  Real Rd_;                   // gas constant of dry air
  Real eps_[1+2*NVAPOR];
  Real rcp_[1+2*NVAPOR];
  Real beta_[1+2*NVAPOR];
  Real t3_[1+NVAPOR];
  Real p3_[1+NVAPOR];

  Real termv_;                // terminal velocity
  Real autoc_;                // auto-conversion time
  Real evapr_;                // evaporation rate
  Real tiny_number_;          // very small number

  Real latent_[1+2*NVAPOR];
  Real rcv_[1+2*NVAPOR];
};

#endif
