#ifndef LOAD_CONFIG
#define LOAD_CONFIG

namespace real_gas_models {

  //================================
  // SetConfigure
  //================================

  class SetConfigure {
    static bool loaded;

    static void loadConf();
    SetConfigure() {}
  public:
    static void setConf();
  };
}

#endif // LOAD_CONFIG

