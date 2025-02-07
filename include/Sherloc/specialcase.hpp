#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>
#include <Sherloc/variant.hpp>
#include <spdlog/spdlog.h>

namespace Sherloc {

class SpecialCase {
public:
  // chromosome
  std::string chr;

  // position
  size_t pos;

  // reference 
  std::string ref;

  // alteration
  std::string alt;

  std::string attribute;

  SpecialCase(const SpecialCase&) = default;
  SpecialCase& operator =(const SpecialCase&) = default;

  SpecialCase(SpecialCase&&) = default;
  SpecialCase& operator =(SpecialCase&&) = default;

  // constructor with subject input vcf
  SpecialCase(const std::string& chr0, const int& pos0, const std::string& ref0, const std::string& alt0, const std::string& attribute0) {
    chr = chr0;
    pos = pos0;
    ref = ref0;
    alt = alt0;

    if(ref.size() != alt.size()){ // indel vcf style to vep style
      pos++;
      ref = ref.substr(1);
      alt = alt.substr(1);
      if(ref.size() == 0) ref = "-";
      if(alt.size() == 0) alt = "-";
    }
    attribute = attribute0;
  }

  bool operator== (SpecialCase& a) {
    if (chr == a.chr && pos == a.pos) {
      return true;
    } else {
      return false;
    }
  }

  bool operator< (SpecialCase& a) {
    if (chr < a.chr) {
      return true;
    } else if (chr == a.chr && pos < a.pos) {
      return true;
    } else {
      return false;
    }
  }

  bool operator> (SpecialCase& a) {
    if (chr > a.chr) {
      return true;
    } else if (chr == a.chr && pos > a.pos) {
      return true;
    } else {
      return false;
    }
  }

  void print() {
    std::cout << chr << '\t' << pos << '\t' << ref << '\t' << alt;
    std::cout << "\n";
  }

  static auto load_special_cases(const std::filesystem::path& specialcase_file) {
    std::vector<SpecialCase> special_case_list;
    if(specialcase_file.empty()){
      SPDLOG_INFO("No specialcase_file loaded");
      return special_case_list;
    }
    std::ifstream fin(specialcase_file);
    if(!fin.is_open()){
      SPDLOG_ERROR("specialcase_file: '{}'can't be opened!",
        specialcase_file.c_str());
      exit(1);
    }
    std::string s;
    while (getline(fin, s)) {
      if (s.starts_with('#'))
        continue;
      std::vector<std::string> col;
      boost::split(col, s, Attr::delimiter('\t'));
      std::string chr_no;
      if (col[0].starts_with("chr")) {
        chr_no = std::move(col[0].substr(3));
      } else {
        chr_no = col[0];
      }
      special_case_list.emplace_back(chr_no, std::stoi(col[1]), col[2], col[3], col[4]);
    }
    return special_case_list;
  }
};

}
