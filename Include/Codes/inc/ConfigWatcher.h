#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
#include <mutex>
#include <atomic>
#include <boost/filesystem.hpp>
#include <nlohmann/json.hpp>

namespace fs = boost::filesystem;
using json = nlohmann::json;

class ConfigWatcher
{
public:
  ConfigWatcher(const std::string &filename)
      : filename_(filename), stopFlag_(false)
  {
    lastWriteTime_ = fs::last_write_time(filename_);
    loadConfig();
  }

  void start();
  void stop();

  json getConfig();

private:
  void loadConfig()
  {
    std::ifstream file(filename_);
    if (!file)
    {
      std::cerr << "[Watcher] Nie można otworzyć pliku: " << filename_ << std::endl;
      return;
    }

    try
    {
      json newConfig;
      file >> newConfig;
      std::lock_guard<std::mutex> lock(mutex_);
      config_ = std::move(newConfig);
    }
    catch (const std::exception &e)
    {
      std::cerr << "[Watcher] Błąd parsowania JSON: " << e.what() << std::endl;
      // zachowujemy starą konfigurację
    }
  }

  std::string filename_;
  json config_;
  std::mutex mutex_;
  time_t lastWriteTime_;
  std::atomic<bool> stopFlag_;
  std::thread watcherThread_;
};
