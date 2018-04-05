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
  Real GetDp() const {return Dp_;}

  // microphysics functions
  void CalculateTP(AthenaArray<Real> const& u);
  void Evaporation(AthenaArray<Real> &u, Real dt);
  void SaturationAdjustment(AthenaArray<Real> &u);
  void Precipitation(AthenaArray<Real> &u, Real dt);
  Real TerminalVelocity(Real dm, Real rho, Real grav) {
    return -1./18.*dm*dm*rhol_*grav/(rho*Kv_);
  }

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

  // conversion function
  void Prim2Hydro(Real const prim[], AthenaArray<Real>& w, int i, int j = 0, int k = 0) const;
  void Hydro2Prim(Real prim[], AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;

  // initialization function
  void DryAdiabat(AthenaArray<Real>& w, Real T0, Real P0, Real grav,
    int is, int ie, int i0, int j = 0, int k = 0, Real ptop = 0., Real pbot = FLT_MAX) const;
  void MoistAdiabat(AthenaArray<Real>& w, Real T0, Real P0, Real grav,
    int is, int ie, int i0, int j = 0, int k = 0, Real ptop = 0., Real pbot = FLT_MAX) const;


private:
  MeshBlock *pmy_block_;
  AthenaArray<bool> recondense_;
  AthenaArray<Real> temp_, pres_;

  Real Rd_;                   // gas constant of dry air
  Real eps_[1+2*NVAPOR];
  Real rcp_[1+2*NVAPOR];
  Real beta_[1+2*NVAPOR];
  Real t3_[1+NVAPOR];
  Real p3_[1+NVAPOR];

  Real autoc_;  // auto-conversion time
  Real Kw_;     // mass diffusion coefficient
  Real Kt_;     // heat diffusion coefficient
  Real Kv_;     // kinematic viscosity
  Real Dp_;     // diameter of precipitation particles
  Real rhol_;   // density of precipitation particles

  Real latent_[1+2*NVAPOR];
  Real rcv_[1+2*NVAPOR];
};

#endif
