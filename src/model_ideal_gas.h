#ifndef MODEL_IDEAL_GAS
#define MODEL_IDEAL_GAS

#include "model_general.h"

namespace real_gas_models {
  class Ideal_gas_equation;

  //================================
  // idealGas
  //================================

  class idealGas : public modelGeneral {
    friend class Ideal_gas_equation;

  private:
    idealGas(modelName mn, std::shared_ptr<constgasparameters> &cgp);

  public:
    bool isValid() const;
    void dynamicflowAccept(DerivateFunctor &df);
    void setTemperature(float v,float p);

    void setVolume(float p,float t);
    void setPressure(float v,float t);
    float getVolume(float p,float t)    const;
    float getPressure(float v,float t)  const;

  };
}

#endif // MODEL_IDEAL_GAS

