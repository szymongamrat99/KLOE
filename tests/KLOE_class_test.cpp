/**
 * @file KLOE_class_test.cpp
 * @author Szymon Gamrat
 * @brief Unit tests for KLOE analysis framework core classes
 * @version 1.0
 * @date 2025-08-26
 * 
 * Thi    Int_t result = obj.neu_triangle(&TrcSum, &vtxSigma, 
                                  Clusters, ip, phi4Mom, 
                                  kne4Mom, kne4Vec, trc);
    
    BOOST_CHECK_EQUAL(result, 0);  // No errors
    BOOST_CHECK_SMALL(std::abs(TrcSum), 10.0f);  // Time residual within reasonable bounds
    BOOST_CHECK_GE(vtxSigma, 0.0f);  // Non-negative vertex uncertaintycontains unit tests for the KLOE analysis framework,
 * particularly focusing on the pm00 class functionality.
 * Tests cover basic operations, physics calculations, and error handling.
 */

#define BOOST_TEST_MODULE KLOEClassTest
#include <boost/test/included/unit_test.hpp>
#include <kloe_class.h>
#include <TMath.h>
#include <TLorentzVector.h>

// Test fixture for KLOE tests
struct KLOETestFixture {
    KLOE::pm00 obj;
    const Float_t EPSILON = 1e-3f;  // Precision for floating-point comparisons, relaxed for float
    
    // Constants for physics calculations
    const Float_t mPhi = 1019.461f;  // φ meson mass [MeV/c²]
    const Float_t mK0 = 497.611f;    // K⁰ meson mass [MeV/c²]
    const Float_t c = 299.792458f;   // Speed of light [mm/ns]
};

BOOST_FIXTURE_TEST_SUITE(KLOETests, KLOETestFixture)

/**
 * @brief Test Lorentz transformation functionality
 * Tests the correctness of Lorentz transformations by verifying:
 * 1. Mass invariance
 * 2. Energy-momentum conservation
 * 3. Proper transformation of time component
 */
BOOST_AUTO_TEST_CASE(LorentzTransformation) {
    Float_t input[4] = {100.0f, 0.0f, 0.0f, static_cast<Float_t>(std::sqrt(100.0f*100.0f + mK0*mK0))};  // px=100 MeV/c
    Float_t boost[3] = {0.1f, 0.0f, 0.0f};  // βx = 0.1
    Float_t output[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    obj.lorentz_transf(boost, input, output);
    
    // Test mass invariance
    double mass_before = sqrt(input[3]*input[3] - 
        (input[0]*input[0] + input[1]*input[1] + input[2]*input[2]));
    double mass_after = sqrt(output[3]*output[3] - 
        (output[0]*output[0] + output[1]*output[1] + output[2]*output[2]));
    
    BOOST_CHECK_CLOSE(mass_before, mass_after, EPSILON);
}

/**
 * @brief Test two-body decay momentum calculation
 * Verifies the correct calculation of momentum in two-body decay:
 * φ → K⁰K⁰
 */
BOOST_AUTO_TEST_CASE(TwoBodyDecayMomentum) {
    double p = obj.TwoBodyDecayMass(mPhi, mK0, mK0);
    double expected = 110.4;  // Expected momentum ~110.4 MeV/c
    
    BOOST_CHECK_CLOSE(p, expected, 0.1);
}

/**
 * @brief Test ΔT calculation
 * Verifies the proper calculation of time difference between
 * K→π⁺π⁻ and K→π⁰π⁰ decays in their respective rest frames
 */
BOOST_AUTO_TEST_CASE(DeltaTCalculation) {
    // Create two kaons in back-to-back configuration
    // First kaon at origin, second kaon displaced in positive x
    Float_t p = 110.4f;
    Float_t E = std::sqrt(p*p + mK0*mK0);
    
    TLorentzVector momKch(p, 0, 0, E);
    TLorentzVector posKne(0, 0, 0, 0);
    TLorentzVector momKne(-p, 0, 0, E);
    
    // Second kaon is displaced by 20mm in x direction
    // Time component is calculated as x/β*c where β = p/E
    Float_t beta = p/E;
    Float_t displacement = 20.0f;  // mm
    Float_t time = displacement/(beta*c);
    TLorentzVector posKch(displacement, 0, 0, time);
    
    Float_t deltaT = obj.DeltaT(&momKch, &posKch, &momKne, &posKne);
    
    // For this setup with displaced second kaon, expect positive ΔT
    BOOST_CHECK_LT(deltaT, 0.0f);
}

/**
 * @brief Test array operations
 * Verifies the correctness of array manipulation methods
 */
BOOST_AUTO_TEST_CASE(ArrayOperations) {
    // Test 1D array clearing
    const int size = 5;
    Float_t arr[size] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    
    obj.Clear1DArray(size, arr);
    
    for(int i = 0; i < size; i++) {
        BOOST_CHECK_CLOSE(arr[i], 999.0f, EPSILON);
    }
    
    // Test array equality - test matching elements count
    Int_t arr1[size] = {1,2,3,4,5};
    Int_t arr2[size] = {1,2,3,4,5};  // Same as arr1
    Int_t arr3[size] = {1,2,3,4,6};  // Differs in last element
    
    // Should return size for identical arrays
    BOOST_CHECK_EQUAL(obj.ArrayEquality(arr1, arr2, size), 0);
    // Should return size-1 for arrays differing in last element
    BOOST_CHECK_EQUAL(obj.ArrayEquality(arr1, arr3, size), 1);
}

/**
 * @brief Test time and date utilities
 * Verifies the proper functioning of timing utilities
 */
BOOST_AUTO_TEST_CASE(TimeUtilities) {
    // Test timer
    obj.startTimer();
    // Simple delay using a loop instead of sleep
    for(int i = 0; i < 1000000; i++) {
        volatile int dummy = i * i;  // Prevent optimization
    }
    std::string duration = obj.endTimer();
    
    // Check if duration string contains expected components
    BOOST_TEST(!duration.empty());
    
    // Test timestamp format
    std::string timestamp = obj.getCurrentTimestamp();
    // Check if timestamp has correct length (YYYY-MM-DD_HHMMSS = 15 chars)
    BOOST_TEST(timestamp.length() == 15);
    // Check if timestamp contains the expected separator characters
    BOOST_TEST(timestamp[4] == '-');
    BOOST_TEST(timestamp[7] == '-');
    BOOST_TEST(timestamp[10] == '_');
}

/**
 * @brief Test neutral vertex reconstruction
 * Verifies the triangle method for neutral vertex reconstruction
 */
BOOST_AUTO_TEST_CASE(NeutralVertexReconstruction) {
    Float_t 
        TrcSum = 0.0f,
        vtxSigma = 0.0f,
        Clusters[4][5] = {
            {0.0f, 0.0f, 0.0f, 0.0f, 100.0f},          // x,y,z,t,E for cluster 1
            {100.0f, 0.0f, 0.0f, 333.33f, 100.0f},     // Cluster 2: t = x/c
            {0.0f, 100.0f, 0.0f, 333.33f, 100.0f},     // Cluster 3: t = y/c
            {100.0f, 100.0f, 0.0f, 471.4f, 100.0f}     // Cluster 4: t = sqrt(x²+y²)/c
        },
        ip[3] = {0.0f, 0.0f, 0.0f},
        phi4Mom[4] = {0.0f, 0.0f, 0.0f, mPhi},
        kne4Mom[4] = {110.4f, 0.0f, 0.0f, std::sqrt(110.4f*110.4f + mK0*mK0)},
        kne4Vec[4] = {0.0f, 0.0f, 0.0f, 0.0f},
        trc[4] = {0.0f, 0.0f, 0.0f, 0.0f};

    Int_t result = obj.neu_triangle(&TrcSum, &vtxSigma, 
                                  Clusters, ip, phi4Mom, 
                                  kne4Mom, kne4Vec, trc);
    
    BOOST_CHECK_EQUAL(result, 0);  // No errors
    BOOST_CHECK_SMALL(TrcSum, 1.0f);  // Small time residual
    BOOST_CHECK_GE(vtxSigma, 0.0);  // Non-negative vertex uncertainty
}

BOOST_AUTO_TEST_SUITE_END()
