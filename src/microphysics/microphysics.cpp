#include "microphysics.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../mesh/mesh.hpp"

Microphysics::Microphysics(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  Rd_ = pin->GetReal("microphysics", "Rd");
  ncycle = pin->GetOrAddInteger("microphysics", "ncycle", 1);

  std::vector<std::string> str;

  // read eps, size should be NVAPOR
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "eps", ""), str);
  for (int n = 0; n < ITR; ++n) eps_[n] = 1.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n) {
    eps_[1+n] = atof(str[n].c_str());
    eps_[1+NVAPOR+n] = atof(str[n].c_str());
  }

  // read rcv, size should be NVAPOR
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "rcv", ""), str);
  for (int n = 0; n < ITR; ++n) rcv_[n] = 1.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n) {
    rcv_[1+n] = atof(str[n].c_str());
    rcv_[1+NVAPOR+n] = atof(str[n].c_str());
  }

  // read latent, size should be NVAPOR, unit is [J/kg]
  if (NVAPOR > 0)
    SplitString(pin->GetOrAddString("microphysics", "latent", ""), str);
  for (int n = 0; n < ITR; ++n) latent_[n] = 0.;
  for (int n = 0; n < std::min(NVAPOR, (int)str.size()); ++n) {
    latent_[1+n] = atof(str[n].c_str());
    latent_[1+NVAPOR+n] = atof(str[n].c_str());
  }
}

Microphysics::~Microphysics()
{}
