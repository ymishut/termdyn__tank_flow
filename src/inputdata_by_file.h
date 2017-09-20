#ifndef INPUTDATA_BY_FILE
#define INPUTDATA_BY_FILE

#include <fstream>
#include <vector>
#include <string>
#include <boost/noncopyable.hpp>

namespace real_gas_models {

  //================================
  // InputData
  //================================

  class InputData : private boost::noncopyable {

  public:
    InputData(const std::string &filename);

    std::vector <float> getConstgasparameters();
    std::pair <float,float> getBalloonVF();
    std::vector <float> getBalloonParameters();
    std::string getFunction();
    std::string getEquation();
    std::string getFlowType();
    static void makeExampleInput();

  private:
    std::ifstream infile_;
    std::vector <std::string> inputParameters_;

  };
}

#endif // INPUTDATA_BY_FILE

