#ifndef MICROPHYSICS_HPP
#define MICROPHYSICS_HPP
#include <cfloat>

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

  // access functions
  Real GetRd() const {return Rd_;}
  Real GetCvRatio(int n) const {return rcv_[n];}
  Real GetCpRatio(int n) const {return rcp_[n];}
  Real GetLatent(int n, Real temp = 0.) const {
    return latent_[n] - beta_[n]*Rd_/eps_[n]*temp;
  }
  Real GetMassRatio(int n) const {return eps_[n];}
  Real GetTerminalVelocity() const {return termv_;}

  // microphysics functions
  void CalculateTP(AthenaArray<Real> const& u);
  void Evaporation(AthenaArray<Real> &u, Real dt);
  void SaturationAdjustment(AthenaArray<Real> &u);
  void Precipitation(AthenaArray<Real> &u, Real dt);

  // thermodynamics functions
  Real Chi(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Cp(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Cv(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real TempV(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Temp(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Theta(Real p0, AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real ThetaV(Real p0, AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real ThetaE(Real p0, AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real MSE(Real grav, AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  void SetPrimitive(Real const prim[], AthenaArray<Real>& w, int i, int j = 0, int k = 0) const;
  void DryAdiabat(AthenaArray<Real>& w, Real T0, Real P0, Real grav,
    int k, int j, int i0, int is, int ie, Real ptop = 0., Real pbot = FLT_MAX) const;
  void MoistAdiabat(AthenaArray<Real>& w, Real T0, Real P0, Real grav,
    int k, int j, int i0, int is, int ie, Real ptop = 0., Real pbot = FLT_MAX) const;

  // state variables
  AthenaArray<Real> T, P;

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

  Real latent_[1+2*NVAPOR];
  Real rcv_[1+2*NVAPOR];
};

#endif
