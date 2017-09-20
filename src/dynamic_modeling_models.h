#ifndef DYNAMIC_MODELING_MODELS
#define DYNAMIC_MODELING_MODELS

#include <memory>
#include <functional>
#include <boost/array.hpp>

#include "dynamic_modeling_diffequations.h"

namespace real_gas_models {

  struct balloon;
  class gasparameters;
  class modelGeneral;
  class idealGas;
  class Redlich_Kwong2;
  class Peng_Robinson;

  //================================
  // dXdtobserver struct
  //================================

  struct dXdtobserver {
    std::shared_ptr <modelGeneral> pmg_;
    void operator ()(const difresult_t &x, float t);
    explicit dXdtobserver(std::shared_ptr <modelGeneral> &pmg);
  };

  //================================
  // DerivateFunctorInflow  !IN
  //================================

  class DerivateFunctorInflow: public DerivateFunctor {
    const balloon &bl_;
    gasparameters &outbl_;

  public:
    static constexpr int pressCoef = -1;                              // process will stoped then pressure in balloon will be
    difEquat_t derivfunc_;                                                //    equal outbl_.pressure*0.98

  public:
    DerivateFunctorInflow(const balloon &bl, gasparameters &outbl, modelGeneral* pmg);

    void getFunctor( idealGas &mg);
    void getFunctor( Redlich_Kwong2 &mg);
    void getFunctor( Peng_Robinson &mg);

  };

  //================================
  // DerivateFunctorOutflow  !OUT
  //================================

  class DerivateFunctorOutflow: public DerivateFunctor {
    const balloon &bl_;
    gasparameters &outbl_;

  public:
    static constexpr int pressCoef = 1;
    difEquat_t derivfunc_;

  public:
    DerivateFunctorOutflow(const balloon &bl, gasparameters &outbl, modelGeneral* pmg);

    void getFunctor( idealGas &mg);
    void getFunctor( Redlich_Kwong2 &mg);
    void getFunctor( Peng_Robinson &mg);
  };
}

#endif // DYNAMIC_MODELING_MODELS

