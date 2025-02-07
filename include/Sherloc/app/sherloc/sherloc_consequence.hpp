#pragma once

#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <map>
#include <algorithm>


namespace Sherloc::app::sherloc {

class SherlocConsequence {
  // TODO: 
public:

  std::map< std::string, std::vector< std::string > > miss_hgvs;

  std::map< std::string, std::vector< int > > miss_pos;

  std::map< std::string, std::vector< int > > null_map;

  std::map< std::string, std::vector< int > > nucli_map;

  SherlocConsequence() {}
  SherlocConsequence(SherlocConsequence&&) = default;
  SherlocConsequence& operator = (SherlocConsequence&&) = default;

  SherlocConsequence(const SherlocConsequence&) = default;
  SherlocConsequence& operator = (const SherlocConsequence&) = default;

  SherlocConsequence(const std::string& filename) {
    std::ifstream is;
    std::vector< std::string > file_vec;
    std::string str;

    is.open(filename);

    while (std::getline(is, str))
      file_vec.emplace_back(str);

    is.close();

    for (auto& file : file_vec) {
      is.open(file);

      while (std::getline(is, str)) {
        std::vector< std::string > vec;
        boost::split(vec, str, boost::is_any_of("\t"));
        if (vec[4] != "pathogenic" && vec[4] != "likely_pathogenic") continue;

        auto it_nucli = nucli_map.find(vec[6]);
        if (it_nucli != nucli_map.end()) {
          it_nucli->second.emplace_back(std::stoi(vec[9]));
          std::sort(it_nucli->second.begin(), it_nucli->second.end());
        } else {
          std::vector< int > v;
          v.emplace_back(std::stoi(vec[9]));
          nucli_map[vec[6]] = v;
        }

        if (vec[1] != ".") {
          std::vector< std::string > hgvsp;
          boost::split(hgvsp, vec[1], boost::is_any_of(":"));
          auto it = miss_hgvs.find(hgvsp[0]);
          if (it != miss_hgvs.end()) it->second.emplace_back(hgvsp[1]);
          else {
            std::vector< std::string > v;
            v.emplace_back(hgvsp[1]);
            miss_hgvs[hgvsp[0]] = v;
          }
        }

        std::vector< std::string > type;
        boost::split(type, vec[5], boost::is_any_of(";"));
        bool miss = false;
        for (auto& i : type)
          if (i == "missense_variant") miss = true;

        if (miss == true) {
          auto it = miss_pos.find(vec[7]);
          if (it == miss_pos.end())
            miss_pos[vec[7]] = std::vector<int>(1, std::stoi(vec[9]));
          else {
            it->second.emplace_back(std::stoi(vec[9]));
            std::sort(it->second.begin(), it->second.end());
          }
        }


        bool null = false;

        for (auto& i : type)
          if (i == "stop_gained" || i == "frameshift_variant") null = true;

        if (null == true) {
          auto it = null_map.find(vec[7]);
          std::vector< std::string > hgvsg;
          boost::split(hgvsg, vec[2], boost::is_any_of("."));
          int pos = 0;
          for (auto& i : hgvsg[1]) {
            if (i < '0' || i > '9')  break;
            pos = pos * 10 + (i - '0');
          }
          if (it != null_map.end()) {
            it->second.emplace_back(pos);
            std::sort(it->second.begin(), it->second.end());
          } else {
            std::vector< int > v;
            v.emplace_back(pos);
            null_map[vec[7]] = v;
          }

        }
      }

      is.close();
    }
  }

  bool get_miss(const std::string& hgvsp) {
    std::vector< std::string > vec;
    boost::split(vec, hgvsp, boost::is_any_of(":"));
    auto it = miss_hgvs.find(vec[0]);
    if (it == miss_hgvs.end())  return false;
    for (auto& i : it->second)
      if (i == vec[1]) return true;

    return false;
  }

  bool get_miss(const std::string& trans, const int& pos) {
    auto it = miss_pos.find(trans);
    if (it == miss_pos.end()) return false;
    auto it2 = std::find(it->second.begin(), it->second.end(), pos);
    if (it2 == it->second.end()) return false;
    return true;
  }

  bool get_nucli(const std::string& gene, const int& pos) {
    auto it = nucli_map.find(gene);
    if (it == nucli_map.end()) return false;
    auto it2 = std::find(it->second.begin(), it->second.end(), pos);
    if (it2 == it->second.end()) return false;
    return true;
  }

  bool get_null(const std::string& trans, const int& pos) {
    auto it = null_map.find(trans);
    if (it == null_map.end())  return false;
    if (it->second[it->second.size() - 1] >= pos)  return true;
    return false;
  }
};

}
