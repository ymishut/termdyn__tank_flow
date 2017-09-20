#ifndef MODEL_PENG_ROBINSON
#define MODEL_PENG_ROBINSON

#include "model_general.h"

namespace real_gas_models {
  class Peng_Robinson_equation;

  //================================
  // Peng_Robinson
  //================================

  class Peng_Robinson :public modelGeneral {
    friend class Peng_Robinson_equation;

  private:
    float modelCoefA_,
          modelCoefB_,
          modelCoefK_;

    Peng_Robinson(modelName mn, std::shared_ptr<constgasparameters> &cgp);

  public:
    bool isValid() const;
    void dynamicflowAccept(class DerivateFunctor &df);
    void setVolume(float p,float t);
    void setPressure(float v,float t);
    float getVolume(float p,float t)    const;
    float getPressure(float v,float t)  const;

    float getCoefA()          const;
    float getCoefB()          const;
    float getCoefK()          const;

  };
}
#endif // MODEL_PENG_ROBINSON

