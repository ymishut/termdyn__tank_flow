#ifndef DYNAMIC_MODELING_DIFFEQUATIONS
#define DYNAMIC_MODELING_DIFFEQUATIONS

#include "model_general.h"
#include "gas_description_static.h"

namespace real_gas_models {

  //================================
  // balloon struct
  //================================

  struct balloon {
    const float capacity,                //m^3
                nozzleSq;                //m^2

    balloon (float V, float F);
  };

  //================================
  // dXdt struct
  //================================

  struct dXdt {
  protected:
    const balloon &bl_;
    gasparameters &outbl_;                                                                      // parameters of gas outside the balloon

    void calculateGW(float &G, float &w, float p_less, float p_more, float vol_const);          // G - [kg/sec]; w - [m/sec]

    dXdt(const balloon &bl, gasparameters &outbl);
    virtual void operator ()(const difresult_t &x, difresult_t &dxdt, float=0) = 0;             // x = { 1/volume , temperature }
  public:

    virtual ~dXdt(){}
  };
}

#endif // DYNAMIC_MODELING_DIFFEQUATIONS

