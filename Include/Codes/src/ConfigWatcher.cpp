#include <ConfigWatcher.h>

void ConfigWatcher::start()
{
  watcherThread_ = std::thread([this]()
      {
            while (!stopFlag_) {
                auto currentWriteTime = fs::last_write_time(filename_);
                if (currentWriteTime != lastWriteTime_) {
                    lastWriteTime_ = currentWriteTime;
                    loadConfig();
                    std::cout << "[Watcher] Konfiguracja przeładowana." << std::endl;
                }
                std::this_thread::sleep_for(std::chrono::seconds(1));
            } });
}

void ConfigWatcher::stop()
{
  stopFlag_ = true;
  if (watcherThread_.joinable())
    watcherThread_.join();
}

json ConfigWatcher::getConfig()
{
  std::lock_guard<std::mutex> lock(mutex_);
  return config_;
}