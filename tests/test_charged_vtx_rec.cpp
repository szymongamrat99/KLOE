#define BOOST_TEST_MODULE ChargedVtxRecTests
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <vector>
#include <array>

#include "../Include/Codes/inc/charged_mom.h"
#include "../Include/Codes/inc/closest_approach.h"
#include "../Include/Codes/inc/const.h"

// ============================================================
//  Helper: create a temporary logger for unit tests
// ============================================================
static ErrorHandling::ErrorLogs &testLogger()
{
  static ErrorHandling::ErrorLogs logger("/tmp/kloe_test_logs/", 0, "");
  return logger;
}

// ============================================================
//  SUITE 1: charged_mom – obliczanie pędu z krzywizny
// ============================================================
BOOST_AUTO_TEST_SUITE(ChargedMomSuite)

BOOST_AUTO_TEST_CASE(mode1_pion_phi0_cot1)
{
  // curv=2, phi=0, cot=1 → px=500, py≈0, pz=500
  auto &logger = testLogger();
  double mom[4] = {};
  KLOE::ChargedVtxRec<double, int>::charged_mom(2.0, 0.0, 1.0, mom, 1, logger);

  BOOST_CHECK_CLOSE(mom[0], 500.0, 1e-4);
  BOOST_CHECK_SMALL(mom[1], 1e-9);
  BOOST_CHECK_CLOSE(mom[2], 500.0, 1e-4);

  const double p2 = mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2];
  const double expectedE = std::sqrt(p2 + PhysicsConstants::mPiCh * PhysicsConstants::mPiCh);
  BOOST_CHECK_CLOSE(mom[3], expectedE, 1e-4);
}

BOOST_AUTO_TEST_CASE(mode2_electron)
{
  auto &logger = testLogger();
  double mom[4] = {};
  KLOE::ChargedVtxRec<double, int>::charged_mom(4.0, M_PI / 2.0, 0.0, mom, 2, logger);

  // px = cos(π/2)*250 ≈ 0, py = sin(π/2)*250 = 250, pz = 0
  BOOST_CHECK_SMALL(mom[0], 0.01);
  BOOST_CHECK_CLOSE(mom[1], 250.0, 1e-4);
  BOOST_CHECK_SMALL(mom[2], 1e-9);

  const double p2 = mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2];
  const double expectedE = std::sqrt(p2 + PhysicsConstants::mElec * PhysicsConstants::mElec);
  BOOST_CHECK_CLOSE(mom[3], expectedE, 1e-4);
}

BOOST_AUTO_TEST_CASE(mode3_muon)
{
  auto &logger = testLogger();
  double mom[4] = {};
  KLOE::ChargedVtxRec<double, int>::charged_mom(5.0, M_PI, 2.0, mom, 3, logger);

  // px = cos(π)*200 = -200, py ≈ 0, pz = 2*200 = 400
  BOOST_CHECK_CLOSE(mom[0], -200.0, 1e-3);
  BOOST_CHECK_SMALL(mom[1], 0.01);
  BOOST_CHECK_CLOSE(mom[2], 400.0, 1e-4);

  const double p2 = mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2];
  const double expectedE = std::sqrt(p2 + PhysicsConstants::mMuon * PhysicsConstants::mMuon);
  BOOST_CHECK_CLOSE(mom[3], expectedE, 1e-4);
}

BOOST_AUTO_TEST_CASE(curv_zero_guard)
{
  // curv=0 → guard clause powinien wrócić bez modyfikacji mom
  auto &logger = testLogger();
  double mom[4] = {-1.0, -1.0, -1.0, -1.0};
  KLOE::ChargedVtxRec<double, int>::charged_mom(0.0, 0.5, 1.0, mom, 1, logger);

  // mom nie powinien być nadpisany
  BOOST_CHECK_EQUAL(mom[0], -1.0);
  BOOST_CHECK_EQUAL(mom[3], -1.0);
}

BOOST_AUTO_TEST_CASE(negative_curvature_uses_abs)
{
  // Ujemna krzywizna → abs(curv) powinien dać te same wielkości pędów
  auto &logger = testLogger();
  double momPos[4] = {}, momNeg[4] = {};
  KLOE::ChargedVtxRec<double, int>::charged_mom(3.0, 0.7, 0.5, momPos, 1, logger);
  KLOE::ChargedVtxRec<double, int>::charged_mom(-3.0, 0.7, 0.5, momNeg, 1, logger);

  for (int i = 0; i < 4; i++)
    BOOST_CHECK_CLOSE(momPos[i], momNeg[i], 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  SUITE 2: IPBoostCorr – przecięcie prostej z płaszczyzną
// ============================================================
BOOST_AUTO_TEST_SUITE(IPBoostCorrSuite)

BOOST_AUTO_TEST_CASE(happy_path)
{
  auto &logger = testLogger();
  KLOE::ChargedVtxRec<double, int> rec(logger);

  // Linia: (0,1,2) + t*(1,0,0), płaszczyzna: normalna (1,0,0) przez (5,0,0)
  // t = ((X_line - X_plane) · n) / (vec_line · n) = (0-5)/1 = -5
  // punkt = (0,1,2) + (-5)*(1,0,0) = (-5,1,2)
  double X_line[3] = {0.0, 1.0, 2.0};
  double vec_line[3] = {1.0, 0.0, 0.0};
  double X_plane[3] = {5.0, 0.0, 0.0};
  double vec_plane[3] = {1.0, 0.0, 0.0};
  double out[3] = {};

  const int rc = rec.IPBoostCorr(X_line, vec_line, X_plane, vec_plane, out);
  BOOST_CHECK_EQUAL(rc, 0);
  BOOST_CHECK_CLOSE(out[0], -5.0, 1e-6);
  BOOST_CHECK_CLOSE(out[1], 1.0, 1e-6);
  BOOST_CHECK_CLOSE(out[2], 2.0, 1e-6);
}

BOOST_AUTO_TEST_CASE(denom_zero)
{
  auto &logger = testLogger();
  KLOE::ChargedVtxRec<double, int> rec(logger);

  // vec_line ⊥ vec_plane → dot=0
  double X_line[3] = {0.0, 0.0, 0.0};
  double vec_line[3] = {1.0, 0.0, 0.0};
  double X_plane[3] = {0.0, 0.0, 0.0};
  double vec_plane[3] = {0.0, 1.0, 0.0};
  double out[3] = {};

  const int rc = rec.IPBoostCorr(X_line, vec_line, X_plane, vec_plane, out);
  BOOST_CHECK_EQUAL(rc, static_cast<int>(ErrorHandling::ErrorCodes::DENOM_EQ_ZERO));
}

BOOST_AUTO_TEST_CASE(vector_overload)
{
  auto &logger = testLogger();
  KLOE::ChargedVtxRec<double, int> rec(logger);

  // Linia: (0,0,0) + t*(0,0,1), płaszczyzna: normalna (0,0,1) przez (0,0,10)
  // t = ((0,0,0)-(0,0,10))·(0,0,1) / ((0,0,1)·(0,0,1)) = -10/1 = -10
  // punkt = (0,0,0) + (-10)*(0,0,1) = (0,0,-10)
  double X_line[3] = {0.0, 0.0, 0.0};
  double vec_line[3] = {0.0, 0.0, 1.0};
  double X_plane[3] = {0.0, 0.0, 10.0};
  double vec_plane[3] = {0.0, 0.0, 1.0};
  std::vector<double> out(3, 0.0);

  const int rc = rec.IPBoostCorr(X_line, vec_line, X_plane, vec_plane, out);
  BOOST_CHECK_EQUAL(rc, 0);
  BOOST_CHECK_CLOSE(out[2], -10.0, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  SUITE 3: closest_approach – punkt najkrótszej odległości
// ============================================================
BOOST_AUTO_TEST_SUITE(ClosestApproachSuite)

BOOST_AUTO_TEST_CASE(parallel_offset_lines)
{
  // Dwie proste wzdłuż osi Z, odległe o 2 w X
  Float_t p1[3] = {0, 0, 0}, v1[3] = {0, 0, 1};
  Float_t p2[3] = {2, 0, 0}, v2[3] = {0, 0, 1};
  Float_t ip[3] = {};

  // Uwaga: linie równoległe → iloczyn wektorowy = 0, wynik niezdefiniowany
  // Test sprawdza, że funkcja nie crashuje
  closest_approach(p1, v1, p2, v2, ip);
  // Nie sprawdzamy wartości, bo n_len_2 = 0 → dzielenie przez zero
  BOOST_CHECK(true); // no crash
}

BOOST_AUTO_TEST_CASE(crossing_lines_at_origin)
{
  // Linia 1: (0,0,0) + t*(1,0,0), linia 2: (0,0,0) + s*(0,1,0)
  Float_t p1[3] = {0, 0, 0}, v1[3] = {1, 0, 0};
  Float_t p2[3] = {0, 0, 0}, v2[3] = {0, 1, 0};
  Float_t ip[3] = {};

  closest_approach(p1, v1, p2, v2, ip);

  BOOST_CHECK_SMALL(ip[0], 1e-5f);
  BOOST_CHECK_SMALL(ip[1], 1e-5f);
  BOOST_CHECK_SMALL(ip[2], 1e-5f);
}

BOOST_AUTO_TEST_CASE(skew_lines_known_result)
{
  // Linia 1: (0,0,0) + t*(1,0,0) → oś X
  // Linia 2: (0,1,0) + s*(0,0,1) → równoległa do Z, przesunięta w Y o 1
  // Najbliższy punkt na L1 to (0,0,0)
  Float_t p1[3] = {0, 0, 0}, v1[3] = {1, 0, 0};
  Float_t p2[3] = {0, 1, 0}, v2[3] = {0, 0, 1};
  Float_t ip[3] = {};

  closest_approach(p1, v1, p2, v2, ip);

  BOOST_CHECK_SMALL(ip[0], 1e-5f);
  BOOST_CHECK_SMALL(ip[1], 1e-5f);
  BOOST_CHECK_SMALL(ip[2], 1e-5f);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  SUITE 4: findKClosestRec – szukanie wierzchołka z 2 trackami
//           Używamy wariantu std::vector + PxTv/PyTv/PzTv,
//           żeby uniknąć zależności od Utils::properties
// ============================================================
BOOST_AUTO_TEST_SUITE(FindKClosestRecSuite)

BOOST_AUTO_TEST_CASE(no_vertices_returns_error)
{
  auto &logger = testLogger();

  Int_t nv = 0, ntv = 0, mode = 1;
  Int_t iv[1] = {0};
  double IP[3] = {0, 0, 0};
  double CurV[1] = {1.0};
  double PxTv[1] = {}, PyTv[1] = {}, PzTv[1] = {};
  double xv[1] = {}, yv[1] = {}, zv[1] = {};

  KLOE::ChargedVtxRec<double, int> rec(nv, ntv, iv, IP, CurV, PxTv, PyTv, PzTv, xv, yv, zv, mode, logger);

  std::vector<double> KchRec(9, 0.0), trk1(4, 0.0), trk2(4, 0.0);
  std::vector<Int_t> vtaken;

  auto err = rec.findKClosestRec(KchRec, trk1, trk2, vtaken, logger);
  BOOST_CHECK_EQUAL(static_cast<int>(err),
                    static_cast<int>(ErrorHandling::ErrorCodes::NO_VTX_WITH_OPPOSITE_TRACKS));
}

BOOST_AUTO_TEST_CASE(single_vtx_two_opposite_tracks)
{
  auto &logger = testLogger();

  // 1 wierzchołek, 2 tracki z przeciwnymi krzywiznami
  // iv = {1, 1} → oba tracki przypisane do wierzchołka 0 (iv jest 1-based)
  Int_t nv = 1, ntv = 2, mode = 1;
  Int_t iv[2] = {1, 1};
  double IP[3] = {0.0, 0.0, 0.0};
  double CurV[2] = {2.0, -2.0}; // przeciwne znaki → ok

  // pędy: trk1 = (100, 200, 50), trk2 = (150, -100, 80)
  double PxTv[2] = {100.0, 150.0};
  double PyTv[2] = {200.0, -100.0};
  double PzTv[2] = {50.0, 80.0};

  double xv[1] = {1.0}, yv[1] = {2.0}, zv[1] = {3.0};

  KLOE::ChargedVtxRec<double, int> rec(nv, ntv, iv, IP, CurV, PxTv, PyTv, PzTv, xv, yv, zv, mode, logger);

  std::vector<double> KchRec(9, 0.0), trk1(4, 0.0), trk2(4, 0.0);
  std::vector<Int_t> vtaken;

  auto err = rec.findKClosestRec(KchRec, trk1, trk2, vtaken, logger);
  BOOST_CHECK_EQUAL(static_cast<int>(err), static_cast<int>(ErrorHandling::ErrorCodes::NO_ERROR));

  // vtaken: wierzchołek 0, tracki 0 i 1
  BOOST_CHECK_EQUAL(vtaken[0], 0);
  BOOST_CHECK_EQUAL(vtaken[1], 0);
  BOOST_CHECK_EQUAL(vtaken[2], 1);

  // KchRec[0..2] = suma pędów tracków
  BOOST_CHECK_CLOSE(KchRec[0], 100.0 + 150.0, 1e-4);
  BOOST_CHECK_CLOSE(KchRec[1], 200.0 + (-100.0), 1e-4);
  BOOST_CHECK_CLOSE(KchRec[2], 50.0 + 80.0, 1e-4);

  // KchRec[3] = suma energii (E_pi = sqrt(p^2 + m_pi^2))
  double E1 = std::sqrt(100.*100. + 200.*200. + 50.*50. + PhysicsConstants::mPiCh * PhysicsConstants::mPiCh);
  double E2 = std::sqrt(150.*150. + 100.*100. + 80.*80. + PhysicsConstants::mPiCh * PhysicsConstants::mPiCh);
  BOOST_CHECK_CLOSE(KchRec[3], E1 + E2, 1e-4);

  // KchRec[4] = |p_total|
  double ptot = std::sqrt(KchRec[0]*KchRec[0] + KchRec[1]*KchRec[1] + KchRec[2]*KchRec[2]);
  BOOST_CHECK_CLOSE(KchRec[4], ptot, 1e-4);

  // KchRec[5] = masa inwariantna = sqrt(E^2 - p^2)
  double minv = std::sqrt(KchRec[3]*KchRec[3] - ptot*ptot);
  BOOST_CHECK_CLOSE(KchRec[5], minv, 1e-4);

  // KchRec[6..8] = pozycja wierzchołka
  BOOST_CHECK_CLOSE(KchRec[6], 1.0, 1e-6);
  BOOST_CHECK_CLOSE(KchRec[7], 2.0, 1e-6);
  BOOST_CHECK_CLOSE(KchRec[8], 3.0, 1e-6);

  // trk1, trk2 powinny mieć pędy i energię
  BOOST_CHECK_CLOSE(trk1[0], 100.0, 1e-4);
  BOOST_CHECK_CLOSE(trk2[0], 150.0, 1e-4);
  BOOST_CHECK_CLOSE(trk1[3], E1, 1e-4);
  BOOST_CHECK_CLOSE(trk2[3], E2, 1e-4);
}

BOOST_AUTO_TEST_CASE(same_sign_curvatures_no_vertex)
{
  // Dwa tracki ze znakiem krzywizny takim samym → nie tworzą wierzchołka
  auto &logger = testLogger();

  Int_t nv = 1, ntv = 2, mode = 1;
  Int_t iv[2] = {1, 1};
  double IP[3] = {0, 0, 0};
  double CurV[2] = {3.0, 3.0}; // ten sam znak!

  double PxTv[2] = {100, 200}, PyTv[2] = {100, 200}, PzTv[2] = {100, 200};
  double xv[1] = {0}, yv[1] = {0}, zv[1] = {0};

  KLOE::ChargedVtxRec<double, int> rec(nv, ntv, iv, IP, CurV, PxTv, PyTv, PzTv, xv, yv, zv, mode, logger);

  std::vector<double> KchRec(9, 0.0), trk1(4, 0.0), trk2(4, 0.0);
  std::vector<Int_t> vtaken;

  auto err = rec.findKClosestRec(KchRec, trk1, trk2, vtaken, logger);
  BOOST_CHECK_EQUAL(static_cast<int>(err),
                    static_cast<int>(ErrorHandling::ErrorCodes::NO_VTX_WITH_OPPOSITE_TRACKS));
}

BOOST_AUTO_TEST_CASE(single_track_per_vtx_no_match)
{
  // 2 wierzchołki, ale każdy ma tylko 1 track → mapTmp[iv] != 2
  auto &logger = testLogger();

  Int_t nv = 2, ntv = 2, mode = 1;
  Int_t iv[2] = {1, 2}; // track 0→vtx0, track 1→vtx1
  double IP[3] = {0, 0, 0};
  double CurV[2] = {2.0, -2.0};

  double PxTv[2] = {100, 200}, PyTv[2] = {100, 200}, PzTv[2] = {100, 200};
  double xv[2] = {0, 5}, yv[2] = {0, 5}, zv[2] = {0, 5};

  KLOE::ChargedVtxRec<double, int> rec(nv, ntv, iv, IP, CurV, PxTv, PyTv, PzTv, xv, yv, zv, mode, logger);

  std::vector<double> KchRec(9, 0.0), trk1(4, 0.0), trk2(4, 0.0);
  std::vector<Int_t> vtaken;

  auto err = rec.findKClosestRec(KchRec, trk1, trk2, vtaken, logger);
  BOOST_CHECK_EQUAL(static_cast<int>(err),
                    static_cast<int>(ErrorHandling::ErrorCodes::NO_VTX_WITH_OPPOSITE_TRACKS));
}

BOOST_AUTO_TEST_CASE(two_vtx_selects_closest_to_ip)
{
  // 2 wierzchołki, oba z 2 trackami o przeciwnych krzywiznach
  // findKClosestRec powinien wybrać wierzchołek bliższy IP
  auto &logger = testLogger();

  Int_t nv = 2, ntv = 4, mode = 1;
  // tracki 0,1 → vtx 0; tracki 2,3 → vtx 1
  Int_t iv[4] = {1, 1, 2, 2};
  double IP[3] = {0, 0, 0};
  double CurV[4] = {2.0, -2.0, 3.0, -3.0};

  double PxTv[4] = {100, 200, 50, 80};
  double PyTv[4] = {100, -100, 50, -80};
  double PzTv[4] = {50, 80, 30, 40};

  // vtx 0 daleko (10,10,10), vtx 1 blisko (1,1,1)
  double xv[2] = {10.0, 1.0};
  double yv[2] = {10.0, 1.0};
  double zv[2] = {10.0, 1.0};

  KLOE::ChargedVtxRec<double, int> rec(nv, ntv, iv, IP, CurV, PxTv, PyTv, PzTv, xv, yv, zv, mode, logger);

  std::vector<double> KchRec(9, 0.0), trk1(4, 0.0), trk2(4, 0.0);
  std::vector<Int_t> vtaken;

  auto err = rec.findKClosestRec(KchRec, trk1, trk2, vtaken, logger);
  BOOST_CHECK_EQUAL(static_cast<int>(err), static_cast<int>(ErrorHandling::ErrorCodes::NO_ERROR));

  // Powinien wybrać vtx 1 (bliższy IP)
  BOOST_CHECK_EQUAL(vtaken[0], 1);

  // Pozycja wierzchołka
  BOOST_CHECK_CLOSE(KchRec[6], 1.0, 1e-6);
  BOOST_CHECK_CLOSE(KchRec[7], 1.0, 1e-6);
  BOOST_CHECK_CLOSE(KchRec[8], 1.0, 1e-6);

  // Pędy tracków z vtx 1 (indeksy 2 i 3)
  BOOST_CHECK_CLOSE(trk1[0], 50.0, 1e-4);
  BOOST_CHECK_CLOSE(trk2[0], 80.0, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()

// ============================================================
//  SUITE 5: KaonMomFromBoost – pęd kaonu z boostu phi
// ============================================================
BOOST_AUTO_TEST_SUITE(KaonMomFromBoostSuite)

BOOST_AUTO_TEST_CASE(basic_positive_result)
{
  auto &logger = testLogger();
  KLOE::ChargedVtxRec<double, int> rec(logger);

  // pKaon: kierunek kaonu (znormalizujemy do osi X), pboost: boost phi
  double pKaon[9] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // phi boost: px=15, py=0, pz=0, E=1020 (blisko progu phi→KK)
  double pboost[4] = {15.0, 0.0, 0.0, 1020.0};
  double pKaonBoost[9] = {};

  int rc = rec.KaonMomFromBoost(pKaon, pboost, pKaonBoost);
  BOOST_CHECK_EQUAL(rc, 0);

  // Wynik: pKaonBoost powinien mieć sensowną masę = mK0
  BOOST_CHECK_CLOSE(pKaonBoost[5], PhysicsConstants::mK0, 1e-4);

  // Energia: E = sqrt(p^2 + mK0^2)
  double p_mod = pKaonBoost[4];
  double E_expected = std::sqrt(p_mod * p_mod + PhysicsConstants::mK0 * PhysicsConstants::mK0);
  BOOST_CHECK_CLOSE(pKaonBoost[3], E_expected, 1e-4);
}

BOOST_AUTO_TEST_CASE(vector_overload_same_result)
{
  auto &logger = testLogger();
  KLOE::ChargedVtxRec<double, int> rec(logger);

  double pKaonArr[9] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 6.0, 7.0};
  std::vector<double> pKaonVec(pKaonArr, pKaonArr + 9);
  double pboost[4] = {15.0, 0.0, 0.0, 1020.0};
  double pKaonBoostArr[9] = {};
  std::vector<double> pKaonBoostVec(9, 0.0);

  int rc1 = rec.KaonMomFromBoost(pKaonArr, pboost, pKaonBoostArr);
  int rc2 = rec.KaonMomFromBoost(pKaonVec, pboost, pKaonBoostVec);

  BOOST_CHECK_EQUAL(rc1, rc2);
  for (int i = 0; i < 9; i++)
    BOOST_CHECK_CLOSE(pKaonBoostArr[i], pKaonBoostVec[i], 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

