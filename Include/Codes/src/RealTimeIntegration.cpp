#include <RealTimeIntegration.h>

API::RTI::RTI(const char *URLschema, std::string parPath) : _URLschema(URLschema), _parPath(parPath)
{
  _curlHandle = curl_easy_init();
}

API::RTI::RTI(TString URLschema, std::string parPath) : _URLschema(URLschema), _parPath(parPath)
{
  _curlHandle = curl_easy_init();
}

API::RTI::RTI()
{
  _curlHandle = curl_easy_init();
}

void API::RTI::setURL(const char *URLschema)
{
  _URLschema = URLschema;
}

void API::RTI::setURL(TString URLschema)
{
  _URLschema = URLschema;
}

void API::RTI::setOpt()
{
  // curl_easy_setopt(_curlHandle, CURLOPT_URL, _URLschema.Data());
  curl_easy_setopt(_curlHandle, CURLOPT_WRITEFUNCTION, &API::RTI::write_to_string);
  curl_easy_setopt(_curlHandle, CURLOPT_WRITEDATA, &_response);
}

void API::RTI::setParsPath(std::string parsPath)
{
  _parPath = parsPath;
}

void API::RTI::setMultiURL(TString propName)
{
  std::ifstream fprop(propName);
  json MultiURL = json::parse(fprop);

  _URLschema = (std::string)MultiURL["PDGURL"];
  _interface = (std::string)MultiURL["interfaces"]["summary"];
  _pdgidMax = (Int_t)MultiURL["pdgids"]["values"].size();

  TString
      combined = "";

  for (Int_t i = 0; i < _pdgidMax; i++)
  {
    _pdgID.push_back((std::string)MultiURL["pdgids"]["values"][i]);
    combined = _URLschema + _interface + _pdgID[i];
    _MultiURL.push_back((std::string)combined);
  }
}

void API::RTI::getAPICall()
{
  _result = curl_easy_perform(_curlHandle);
}

void API::RTI::getMultiAPICall()
{
  for (Int_t i = 0; i < _pdgidMax; i++)
  {
    curl_easy_setopt(_curlHandle, CURLOPT_URL, _MultiURL[i].c_str());
    _result = curl_easy_perform(_curlHandle);
    _resptable.push_back(_response);
    _response = "";
  }
}

Double_t API::RTI::getValue(Int_t i)
{
  _responseJSON = json::parse(_resptable[i]);

  Int_t 
      length = _responseJSON["pdg_values"].size(),
      nonCPTindex = 0;

  if(length > 1)
  {
    for (Int_t i = 0; i < length; i++)
      if(_responseJSON["pdg_values"][i]["comment"] == "Not assuming CPT")
        nonCPTindex = i;
  };

  return _responseJSON["pdg_values"][nonCPTindex]["value"];
}

void API::RTI::createJSONPars()
{
  std::ofstream finalfile(_parPath);
  
  for(Int_t i = 0; i < _pdgidMax; i++)
    _parametersJSON["values"][_pdgID[i]] = getValue(i);
  
  finalfile << _parametersJSON.dump(4);
  finalfile.close();
}

size_t API::RTI::write_to_string(void *ptr, size_t size, size_t count, void *stream)
{
  ((std::string *)stream)->append((char *)ptr, 0, size * count);
  return size * count;
}
