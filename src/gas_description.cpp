#include "gas_description.h"
#include "models_math.h"
#include "models_exceptions.h"

#include <iostream>
#include <boost/array.hpp>

//================================
// constgasparameters ctor
//================================

real_gas_models::constgasparameters::constgasparameters(float vk, float pk, float tk, float mol,
                                                float adi, float cv, float b, float R, float af)
  :V_K(vk),P_K(pk),T_K(tk),molecularmass(mol),adiabatic_index(adi),
                         cv(cv),beta_kr(b),R(R),acentricfactor(af) {}

//================================
// constgasparameters getInstance
//================================

std::shared_ptr<real_gas_models::constgasparameters> real_gas_models::constgasparameters::getInstance(float vk,
                                                 float pk, float tk, float mol, float adi, float cv,float af) {
  bool hasError = ((adi==1.0)||(adi<=0.0)||(vk<=0.0)||(pk<=0.0)||
                   (tk<=0.0)||(cv<=0.0)||(mol<=0.0)||(af<=0.0));

  hasError = hasError || (((!std::isfinite(adi))||(!std::isfinite(vk))||(!std::isfinite(pk))||(!std::isfinite(tk))
                           ||(!std::isfinite(cv))||(!std::isfinite(mol))||(!std::isfinite(af))));
  if(hasError)
    throw modelExceptions("Negative or zero parameters in struct constgasparameters");
  const float tempb = std::pow(2.0/(adi+1.0),adi/(adi-1.0)),
              tempR = 8314.4599/mol;
  return std::shared_ptr<constgasparameters>(new constgasparameters(vk,pk,tk,mol,adi,cv,tempb,tempR,af));
}

//================================
// constgasparameters getInstance
//================================

std::shared_ptr<real_gas_models::constgasparameters> real_gas_models::constgasparameters::getInstance(std::vector<float> vec) {
  return getInstance(vec[0],vec[1],vec[2],vec[3],vec[4],vec[5],vec[6]);
}

//================================
// parameters ctors
//================================

real_gas_models::parameters::parameters(float v, float p, float t)
  :volume(v),pressure(p),temperature(t){}

real_gas_models::parameters::parameters()
  :parameters(0,0,0){}

