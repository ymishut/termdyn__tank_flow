#ifndef DYNAMIC_MODELING_CALCULATE
#define DYNAMIC_MODELING_CALCULATE

#include "dynamic_modeling_models.h"
#include "model_general.h"

namespace real_gas_models {


  //================================
  // setStartVariables templates
  //================================


  template <class DerivateFunctorType>
  void setStartVariables (float &pm, float &pl, float ac);

  // for inflow outbl>pin
  template <>
  inline void setStartVariables<DerivateFunctorInflow> (float &pm, float &pl,float ac) {
    (void) pl;
    (void) pm;
    pm*=(1.0-ac);                                      // now pm < 1
  }

  template <>
  inline void setStartVariables<DerivateFunctorOutflow> (float &pm, float &pl, float ac) {
    std::swap(pl,pm);
    pl*=(1.0+ac);              // now pl > 1
  }

}
#endif // DYNAMIC_MODELING_CALCULATE

