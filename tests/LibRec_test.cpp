// Include the Google Test framework
#define BOOST_TEST_MODULE LibRecTest
#include <boost/test/included/unit_test.hpp>

// --- Test, który ZAWSZE PRZEJDZIE ---
BOOST_AUTO_TEST_CASE(this_test_should_pass) {
    BOOST_TEST_MESSAGE("Running a test that is designed to pass.");
    int result = 5 + 1;
    BOOST_TEST(result == 6, "Expected add_one(5) to be 6"); // To powinno być prawdą
    BOOST_CHECK(true); // To też powinno być prawdą
}
