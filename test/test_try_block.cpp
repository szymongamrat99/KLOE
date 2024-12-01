#include <iostream>
#include <json.hpp>

using namespace std;
using json = nlohmann::json;

namespace ns
{
  struct person
  {
    std::string name;
    std::string address;
    int age;
  };

  void to_json(json &j, const person &p)
  {
    j = json{{"name", p.name}, {"address", p.address}, {"age", p.age}};
  }

  void from_json(const json &j, person &p)
  {
    j.at("name").get_to(p.name);
    j.at("address").get_to(p.address);
    j.at("age").get_to(p.age);
  }
};

int main()
{
  ns::person p = {"Andrzej Wajda", "Zwyciestwa 22/4", 82};

  json j;

  j["Name"] = "Andrzej Wajda";
  j["Addresses"]["Zwycięstwa 22/4"] = "Zwycięstwa 22/4";
  j["Addresses"]["Dupy 22/4"] = "Dupy 22/4";
  j["Cars"] = "Old";

  std::cout << j;

  return 0;
}