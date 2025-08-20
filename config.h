#pragma once

#include <iostream>
#include <type_traits>
#include "enum.h"

template <typename T>
concept Enum = std::is_enum_v<T>;

template<Enum E>
std::ifstream& operator>>(std::ifstream& is, E& eval) {
    typename std::underlying_type<E>::type val;
    if (is >> val) eval = static_cast<E>(val);
    return is;
}

struct Config {
    Config() = default;
    Config(const std::string& fileName) {
        std::ifstream is(fileName);
        is >> mode >> dist >> cdist >> prec
            >> nsrcs >> L >> EPS >> maxNodeParts  >> evalDirect;
    }

    Mode mode;
    Dist dist;
    ChargeDist cdist;
    Precision prec;
    int nsrcs;
    double L;
    double EPS;
    int maxNodeParts;
    bool evalDirect;
};

const std::string makeFname(const Config& config) {
    std::string distStr =
        [&]() -> std::string {
            switch (config.dist) {
                case Dist::UNIFORM:  return "uniform";
                case Dist::GAUSSIAN: return "gauss";
                case Dist::GRID:     return "grid";
            }
        }();
    std::string cdistStr = 
        [&]() -> std::string { 
            switch (config.cdist) {
                case ChargeDist::PLUS:  return "plus";
                case ChargeDist::MINUS: return "minus"; 
                case ChargeDist::DIP:   return "dip";
                case ChargeDist::QUAD:  return "quad";
                case ChargeDist::OCT:   return "oct";
                case ChargeDist::RAND:  return "rand";
            }
        }();

    return "config/part3D/" + distStr + "_" + cdistStr + ".txt";
}