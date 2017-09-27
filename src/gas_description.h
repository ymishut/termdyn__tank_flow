#ifndef SRC_GAS_DESCRIPTION_H_
#define SRC_GAS_DESCRIPTION_H_

#include <memory>
#include <vector>
#include <string>

#include <boost/noncopyable.hpp>

namespace real_gas_models {
  //================================
  // state_phase enum || stateToString
  //================================
  // SCF: t>T_K, p>P_K;    GAS: p_binodal < p < p_K, t>t_binodal;
  
  enum class state_phase {SCF = 0, LIQUID, LIQ_STEAM, GAS};    

  // LIQUID: p<P_K; v<vleft;     
  // in perspective:  LIQ_STEAM: p<p_binodal, vleft < v < vrigth;
  // There are not metastable states (between binodal and spinodal)
  
  static const std::vector<std::string> stateToString {
    "SCF", "LIQUID", "LIQ_STEAM", "GAS"};           

  //================================
  // parameters struct
  //================================

  struct parameters {
    float   volume,               // m^3 / kg
            pressure,             // Pa
            temperature;          // K

    parameters(float v, float p, float t);
    parameters();
  };

  //================================
  // constgasparameters struct
  //================================

  // gas paramaters depending on the physics characteristics
  struct constgasparameters : private boost::noncopyable {           
    const float  V_K,              // K point parameters (critical point)
                 P_K,
                 T_K,
                 molecularmass,
                 adiabatic_index,                    // not const IRL
                 cv,               // heat capacity      not const IRL
                 beta_kr,          // beta_kr=beta_kr(adiabatic_index)
                                   // = [0.0, ... ,1.0]
                                     // if (pressure_1/pressure_2 < beta_kr):
                                     // flow velocity = const, (maximum)
                                     // else : flow velocity ~ p,adiabatic_index and etc
                                     // P.S. look dynamic_modeling*.*
                 R,                // gas constant
                 acentricfactor;

  public:
    static std::shared_ptr<constgasparameters> getInstance(float vk, float pk,
                           float tk, float mol, float adi, float cv, float af);
    static std::shared_ptr<constgasparameters> getInstance(
                                     std::vector<float> vec);

  private:
    constgasparameters(float vk, float pk, float tk, float mol, float adi,
                                      float cv, float b, float R, float af);
  };

  //================================
  // Igasparameters interface
  //================================

    class Igasparameters {
    public:
      virtual float cgetV_K()            const = 0;
      virtual float cgetP_K()            const = 0;
      virtual float cgetT_K()            const = 0;
      virtual float cgetMolecularMass()  const = 0;
      virtual float cgetAdiabatic()      const = 0;
      virtual float cgetCV()             const = 0;
      virtual float cgetBeta()           const = 0;
      virtual float cgetR()              const = 0;
      virtual float cgetAcentricFactor() const = 0;
      virtual float cgetVolume()         const = 0;
      virtual float cgetPressure()       const = 0;
      virtual float cgetTemperature()    const = 0;
      virtual parameters cgetParameters()const = 0;
      virtual state_phase cgetState()    const = 0;
      virtual std::shared_ptr<constgasparameters> cgetConstparameters() 
                                                               const = 0;

      virtual void csetParameters(float, float, float, state_phase) = 0;
      virtual ~Igasparameters() {}
    };
}  // namespace real_gas_models

#endif  // SRC_GAS_DESCRIPTION_H_

