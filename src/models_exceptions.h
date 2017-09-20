#ifndef MODELS_EXCEPTIONS
#define MODELS_EXCEPTIONS

#include <exception>
#include <string>

namespace real_gas_models {

class modelExceptions: public std::exception {
protected:
  std::string message_;
public:

  //================================
  // modelExceptions
  //================================

  template <class T>
  explicit modelExceptions(T&& m)
    :message_(std::forward<T>(m)) {}

  virtual const char* what() const noexcept{
    return message_.c_str();
  }
  virtual ~modelExceptions() noexcept{}
};
}

#endif // MODELS_EXCEPTIONS

