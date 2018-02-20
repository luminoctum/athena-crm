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

  // access functions
  Real GetRd() const {return Rd_;}
  Real GetCvRatio(int n) const {return rcv_[n];}
  Real GetCpRatio(int n) const {return rcp_[n];}
  Real GetLatent(int n) const {return latent_[n];}
  Real GetMassRatio(int n) const {return eps_[n];}
  Real GetTerminalVelocity() const {return termv_;}

  // microphysics functions
  void CalculateTemperature(AthenaArray<Real> const& u);
  void Evaporation(AthenaArray<Real> &u, Real dt);
  void SaturationAdjustment(AthenaArray<Real> &u);
  void Precipitation(AthenaArray<Real> &u, Real dt);

  // thermodynamics functions
  Real Chi(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Cp(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Cv(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Tempv(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Temp(AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Theta(Real p0, AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real Thetav(Real p0, AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  Real MSE(Real grav, AthenaArray<Real> const& w, int i, int j = 0, int k = 0) const;
  void DryAdiabat(AthenaArray<Real>& w, Real T0, Real P0, Real grav,
    int k, int j, int i0, int is, int ie, bool isothermal = false) const;

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

// frequently used inline functions
inline Real SatVaporPresIdeal(Real t, Real p3, Real beta, Real delta) {
  return p3*exp((1. - 1./t)*beta - delta*log(t));
}

// alpha = L/cv evaluated at current temperature
inline Real GasCloudIdeal(Real const prim[], int iv, int ic,
  Real p3, Real t3, Real alpha, Real beta, Real delta)
{
  Real xv = prim[iv];
  Real xc = prim[ic];
  Real t = prim[IDN]/t3;
  Real s = SatVaporPresIdeal(t,p3,beta,delta)/prim[IPR];

  // if saturation vapor pressure is larger than the total pressure
  // evaporate all condensates
  if (s > 1.) return -xc;

  Real g = 1.;
  for (int n = ICD; n < ICD + NVAPOR; ++n) g -= prim[n];
  g -= xv;

  Real s1 = s/(1. - s);
  Real rate = (xv - s1*g)/(1. + alpha*g*(beta/t - delta)*s1/(1. - s));

  // condensate at most xv vapor
  if (rate > 0.) rate = std::min(rate, xv);

  // evaporate at most xc cloud
  if (rate < 0.) rate = std::max(rate, -xc);
  return rate;
}

#endif
