#pragma once

#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/disease_database.hpp>
#include <Sherloc/DB/vcf.hpp>
#include <Sherloc/Patient/other_patient.hpp>
#include <Sherloc/Attr/allele.hpp>
#include <nlohmann/json.hpp>

namespace Sherloc::Patient {

class Patient {
public:
  std::string vcf_file;
  std::string name;
  std::vector<SherlocMember> sher_mems;
  bool sex;
  bool is_sick;
  bool is_denovo;
  DB::DataBaseVcf dad_vcf, mom_vcf;
  bool is_dad_sick, is_mom_sick;
  OtherPatient healthy_family, sick_family;
  Patient() : is_sick(false), is_dad_sick(false), is_mom_sick(false) {}

  Patient(const Patient& rhs) = default;
  Patient& operator =(const Patient& rhs) = default;

  Patient(Patient&& rhs) = default;
  Patient& operator =(Patient&& rhs) = default;

  Patient(const nlohmann::json& js) {
    is_dad_sick = js["is_dad_sick"];
    is_mom_sick = js["is_mom_sick"];
    sex = js["sex"];
    is_sick = js["is_sick"];
    is_denovo = js["is_denovo"];

    std::vector< std::string > v;
    vcf_file = js["patient_vcf"];
    boost::split(v, vcf_file, boost::is_any_of("/"));
    name = v.back();

    dad_vcf.from(js["dad_vcf"]);
    mom_vcf.from(js["mom_vcf"]);

    std::vector< std::string > vec;
    for (size_t i{}; i < js["healthy_family"].size(); ++i) {
      vec.emplace_back(js["healthy_family"][i]);
    }
    healthy_family = OtherPatient(vec);
    vec.clear();

    for (size_t i{}; i < js["sick_family"].size(); ++i) {
      vec.emplace_back(js["sick_family"][i]);
    }
    sick_family = OtherPatient(vec);
  }

  auto check_allele_denovo(const SherlocMember& sher_mem, bool sex = false){
    using namespace Attr;
    if(sex){
      // male variant on chrY must from dad or denovo
      if(sher_mem.chr == "Y" and !dad_vcf.empty()){ 
        auto dad_allele = dad_vcf.find(sher_mem);
        return dad_allele.empty ? Allele::IsDeNovo : Allele::NotDeNovo;
      }

      // male variant on chrX must from mom or denovo
      if(sher_mem.chr == "X" and !mom_vcf.empty()){
        auto mom_allele = mom_vcf.find(sher_mem);
        return mom_allele.empty ? Allele::IsDeNovo : Allele::NotDeNovo;
      }
    }

    // TODO: the probability of both chromosomes have denovo variant on the same allele is low
    // Should we consider that?

    // Unknown if parents' vcfs not available
    // Just to distinguish the situation that is "Not found" or "VCF itself is not provided"
    if(dad_vcf.empty() and mom_vcf.empty())
      return Allele::Unknown;

    // other situation: heterozygous allele and both parents vcfs are available.
    auto dad_allele = dad_vcf.find(sher_mem);
    auto mom_allele = mom_vcf.find(sher_mem);
    return (!dad_allele.empty or !mom_allele.empty) ? Allele::NotDeNovo : Allele::IsDeNovo;
  }
  
  auto check_allele_denovo(size_t idx){
    return check_allele_denovo(sher_mems[idx], this->sex);
  }
};

}
