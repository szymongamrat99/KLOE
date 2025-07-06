// ConstantsKeys.h
#pragma once

#include <string>

/**
 * @brief Provides string constants for accessing physical constants in JSON.
 * Using these constants prevents typos and enables easier refactoring.
 * All keys are defined using JSON Pointer syntax for clarity in nesting.
 */
namespace ConstantsKeys
{
  static const std::string Re = "/S013EPS";
  static const std::string Im = "/S013EPI";
  static const std::string K0mass = "/S011M";
  static const std::string TauS = "/S012T";
  static const std::string TauL = "/S013T";
  static const std::string deltaM = "/S013D";
  static const std::string modEps = "/S013EP";
  static const std::string phiPM = "/S013F+-";
  static const std::string phi00 = "/S013FOO";
}