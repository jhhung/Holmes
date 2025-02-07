#pragma once

#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/DB/vcf.hpp>
#include <nlohmann/json.hpp>

namespace Sherloc::Patient {

class OtherPatient {
  // TODO: 
public:

  std::vector<DB::DataBaseVcf> op_vec;

  OtherPatient() = default;

  OtherPatient(const std::vector< std::string >& vec) {
    for (auto& i : vec){
      op_vec.emplace_back();
      op_vec.back().from(i);
    }
  }

  OtherPatient(OtherPatient&&) = default;
  OtherPatient& operator = (OtherPatient&&) = default;

  OtherPatient(const OtherPatient&) = default;
  OtherPatient& operator = (const OtherPatient&) = default;

  OtherPatient(const nlohmann::json& js) {
    for (size_t i {}; i < js.size(); ++i) {
      op_vec.emplace_back();
      op_vec.back().from(js[i]);
    }
  }

  int get_observation(const SherlocMember& sher_mem) {
    int num = 0;
    for (auto& i : op_vec)  if (i.find(sher_mem).empty == false) ++num;
    return num;
  }
};

}
