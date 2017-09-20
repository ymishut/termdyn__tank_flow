#ifndef MODEL_REDLICH_KWONG_H
#define MODEL_REDLICH_KWONG_H

#include "model_general.h"

namespace real_gas_models {
  class Redlich_Kwong_equation;

  //================================
  // Redlich_Kwong2
  //================================

  class Redlich_Kwong2 :public modelGeneral {
    friend class Redlich_Kwong_equation;

  private:
    float modelCoefA_,
          modelCoefB_;

    Redlich_Kwong2(modelName mn, std::shared_ptr<constgasparameters> &cgp);

  public:
    bool isValid() const;
    void dynamicflowAccept(class DerivateFunctor &df);
    void setVolume(float p,float t);
    void setPressure(float v,float t);
    float getVolume(float p,float t)    const;
    float getPressure(float v,float t)  const;

    float getCoefA()   const;
    float getCoefB()   const;

  };
}

#endif // MODEL_REDLICH_KWONG_H

