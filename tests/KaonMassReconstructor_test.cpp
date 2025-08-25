/**
 * @file KaonMassReconstructor_test.cpp
 * @brief Unit tests for KaonMassReconstructor class
 */

#define BOOST_TEST_MODULE KaonMassReconstructorTest
#include <boost/test/included/unit_test.hpp>
#include "KaonMassReconstructor.h"
#include "klspm00.hpp"
#include <cmath>

struct KaonMassReconstructorFixture {
    KLOE::pm00 obj;
    const Double_t mPiCh = 139.57061;  // Charged pion mass [MeV/c²]
    const Double_t mK0 = 497.611;      // K⁰ meson mass [MeV/c²]
    const Double_t EPSILON = 1e-4;     // Precision for floating-point comparisons

    // Helper function to create test kaon data
    std::vector<Double_t> createKaonBoost(Double_t px, Double_t py, Double_t pz, 
                                        Double_t E, Double_t vx, Double_t vy, Double_t vz) {
        return {px, py, pz, E, 0.0, 0.0, vx, vy, vz};  // 4-momentum and vertex
    }

    // Helper function to create test tracks
    std::vector<std::vector<Double_t>> createTracks(
        Double_t px1, Double_t py1, Double_t pz1, Double_t E1,
        Double_t px2, Double_t py2, Double_t pz2, Double_t E2) {
        return {
            {px1, py1, pz1, E1},
            {px2, py2, pz2, E2}
        };
    }

    // Helper function to check if two vectors are approximately equal
    void checkVectorsEqual(const std::vector<Double_t>& expected, 
                          const std::vector<Double_t>& actual, 
                          Double_t tolerance = 1e-4) {
        BOOST_REQUIRE_EQUAL(expected.size(), actual.size());
        for(size_t i = 0; i < expected.size(); ++i) {
            BOOST_CHECK_CLOSE(expected[i], actual[i], tolerance);
        }
    }
};

BOOST_FIXTURE_TEST_SUITE(KaonMassReconstructorTests, KaonMassReconstructorFixture)

BOOST_AUTO_TEST_CASE(BackToBackPionsTest) {
    // Test case with back-to-back pions in kaon rest frame
    Double_t p = 110.4;  // Momentum magnitude for K⁰→π⁺π⁻
    Double_t E = std::sqrt(p*p + mK0*mK0);

    // Kaon moving along x-axis
    auto kaonBoost = createKaonBoost(p, 0, 0, E, 20.0, 0, 0);
    std::vector<Double_t> ip = {0, 0, 0};  // Interaction point at origin

    // Two pions with opposite momenta
    auto tracks = createTracks(
        p/2, 0, 0, std::sqrt(pow(p/2, 2) + mPiCh*mPiCh),
        -p/2, 0, 0, std::sqrt(pow(p/2, 2) + mPiCh*mPiCh)
    );

    auto result = KLOE::KaonMassReconstructor::reconstructKaonMass(
        kaonBoost, ip, tracks, mPiCh, obj
    );

    // Check kaon mass
    BOOST_CHECK_CLOSE(result.KaonTwoBody[5], mK0, EPSILON);
    
    // Check momentum conservation
    Double_t totalPx = result.track1TwoBody[0] + result.track2TwoBody[0];
    Double_t totalPy = result.track1TwoBody[1] + result.track2TwoBody[1];
    Double_t totalPz = result.track1TwoBody[2] + result.track2TwoBody[2];
    
    BOOST_CHECK_SMALL(totalPx - kaonBoost[0], EPSILON);
    BOOST_CHECK_SMALL(totalPy, EPSILON);
    BOOST_CHECK_SMALL(totalPz, EPSILON);
}

BOOST_AUTO_TEST_CASE(PerpendicularPionsTest) {
    // Test case with pions at 90 degrees in kaon rest frame
    Double_t p = 110.4;
    Double_t E = std::sqrt(p*p + mK0*mK0);
    Double_t pPi = p/std::sqrt(2.0);  // Pion momentum

    auto kaonBoost = createKaonBoost(0, p, 0, E, 0, 20.0, 0);
    std::vector<Double_t> ip = {0, 0, 0};

    // Pions at 90 degrees to each other
    auto tracks = createTracks(
        pPi, 0, 0, std::sqrt(pPi*pPi + mPiCh*mPiCh),
        0, pPi, 0, std::sqrt(pPi*pPi + mPiCh*mPiCh)
    );

    auto result = KLOE::KaonMassReconstructor::reconstructKaonMass(
        kaonBoost, ip, tracks, mPiCh, obj
    );

    // Check kaon mass
    BOOST_CHECK_CLOSE(result.KaonTwoBody[5], mK0, EPSILON);

    // Check that pions are perpendicular in CM frame
    Double_t dotProduct = 
        result.track1TwoBody[0] * result.track2TwoBody[0] +
        result.track1TwoBody[1] * result.track2TwoBody[1] +
        result.track1TwoBody[2] * result.track2TwoBody[2];
    
    BOOST_CHECK_SMALL(dotProduct, EPSILON);
}

BOOST_AUTO_TEST_CASE(BoostInvariantMassTest) {
    // Test that reconstructed mass is boost invariant
    std::vector<Double_t> boosts = {0.1, 0.3, 0.5, 0.7, 0.9};  // Different boost velocities
    Double_t p = 110.4;
    Double_t baseMass = 0;

    for(Double_t beta : boosts) {
        Double_t gamma = 1.0/std::sqrt(1.0 - beta*beta);
        Double_t E = gamma * std::sqrt(p*p + mK0*mK0);
        Double_t pz = gamma * beta * std::sqrt(p*p + mK0*mK0);

        auto kaonBoost = createKaonBoost(0, 0, pz, E, 0, 0, 20.0);
        std::vector<Double_t> ip = {0, 0, 0};

        // Back-to-back pions in rest frame
        auto tracks = createTracks(
            p/2, 0, 0, std::sqrt(pow(p/2, 2) + mPiCh*mPiCh),
            -p/2, 0, 0, std::sqrt(pow(p/2, 2) + mPiCh*mPiCh)
        );

        auto result = KLOE::KaonMassReconstructor::reconstructKaonMass(
            kaonBoost, ip, tracks, mPiCh, obj
        );

        if(baseMass == 0) {
            baseMass = result.KaonTwoBody[5];
        } else {
            BOOST_CHECK_CLOSE(result.KaonTwoBody[5], baseMass, EPSILON);
        }
    }
}

BOOST_AUTO_TEST_CASE(EnergyConservationTest) {
    Double_t p = 110.4;
    Double_t E = std::sqrt(p*p + mK0*mK0);

    auto kaonBoost = createKaonBoost(p, p, 0, std::sqrt(2*p*p + mK0*mK0), 20.0, 20.0, 0);
    std::vector<Double_t> ip = {0, 0, 0};

    auto tracks = createTracks(
        p/2, p/2, 0, std::sqrt(p*p/2 + mPiCh*mPiCh),
        -p/2, -p/2, 0, std::sqrt(p*p/2 + mPiCh*mPiCh)
    );

    auto result = KLOE::KaonMassReconstructor::reconstructKaonMass(
        kaonBoost, ip, tracks, mPiCh, obj
    );

    Double_t totalE = result.track1TwoBody[3] + result.track2TwoBody[3];
    BOOST_CHECK_CLOSE(totalE, kaonBoost[3], EPSILON);
}

BOOST_AUTO_TEST_SUITE_END()
