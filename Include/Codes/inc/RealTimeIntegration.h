#ifndef REALTIMEINTEGRATION_H
#define REALTIMEINTEGRATION_H

#include <curl/curl.h>
#include <json.hpp>
#include <fstream>

#include <TString.h>
#include <iostream>

using json = nlohmann::json;

namespace API
{
  class RTI
  {
  private:
    CURL *
        _curlHandle;
    CURLcode
        _result;
    TString
        _URLschema,
        _interface;

    std::vector<std::string>
        _MultiURL,
        _pdgID;

    std::vector<std::string>
                      _resptable;

    std::string
            _response,
            _parPath;

    Int_t
        _pdgidMax;

    json
      _responseJSON,
      _parametersJSON;

  public:
    RTI(const char *URLschema, std::string parPath);
    RTI(TString URLschema, std::string parPath);
    RTI();

    void setURL(const char *URLschema);
    void setURL(TString URLschema);
    void setOpt();

    void setParsPath(std::string parsPath);

    void setMultiURL(TString Paths::propName);

    void getAPICall();
    void getMultiAPICall();
    Double_t getValue(Int_t i);

    void createJSONPars();

    static size_t write_to_string(void *ptr, size_t size, size_t count, void *stream);
  };

}

#endif // !REALTIMEINTEGRATION_H