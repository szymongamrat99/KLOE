#include <PdgApiClient.h>
#include <stdexcept>

PdgApiClient::PdgApiClient() {
    // curl_easy_init można robić tutaj lub w metodzie
}

PdgApiClient::~PdgApiClient() {}

size_t PdgApiClient::WriteCallback(void* contents, size_t size, size_t nmemb, std::string* s) {
    size_t newLength = size * nmemb;
    s->append((char*)contents, newLength);
    return newLength;
}

std::string PdgApiClient::fetchParticleJson(const std::string& pdgid, int year) {
    CURL* curl = curl_easy_init();
    std::string response;

    if (curl) {
        // Dynamiczne budowanie URL na podstawie Twoich potrzeb
        std::string url = "https://pdgapi.lbl.gov/summaries/" + pdgid + "/" + std::to_string(year);
        
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
        curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L); // Timeout 10s

        CURLcode res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);

        if (res != CURLE_OK) {
            throw std::runtime_error("CURL failed: " + std::string(curl_easy_strerror(res)));
        }
    }
    else{
      throw std::runtime_error("Failed to initialize CURL");
    }
    return response;
}