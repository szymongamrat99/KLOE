// tests/kloe_class.cpp
//
// This file contains unit tests for the kloe_class functions using Google Test.

// Include the Google Test framework
#define BOOST_TEST_MODULE LibRecTest
#include <boost/test/included/unit_test.hpp>

// Include the header of the code we want to test
#include <kloe_class.h>

struct pm00Fixture
{
  // Konstruktor armatury (Setup): Wywoływany przed każdym testem w tym suite.
  pm00Fixture()
  {
    BOOST_TEST_MESSAGE("Setting up KLOE::pm00 fixture...");
    analyzer = new KLOE::pm00(); // Utwórz nową instancję KLOE::pm00
                                 // Możesz tutaj dodać domyślne dane, jeśli większość testów ich potrzebuje
  }

  // Destruktor armatury (TearDown): Wywoływany po każdym teście w tym suite.
  ~pm00Fixture()
  {
    BOOST_TEST_MESSAGE("Tearing down KLOE::pm00 fixture...");
    delete analyzer; // Posprzątaj instancję KLOE::pm00
    analyzer = nullptr;
  }

  // Obiekt(y), które mają być testowane i dostępne dla przypadków testowych.
  KLOE::pm00 *analyzer;
};

// BOOST_AUTO_TEST_CASE defines an individual test case.
// It will be automatically discovered by Boost.Test runner.
BOOST_FIXTURE_TEST_CASE(test_add_positive_numbers, pm00Fixture)
{
  // BOOST_TEST is a universal assertion.
  // It is equivalent to EXPECT_EQ in GTest.
  BOOST_TEST( true );
};
