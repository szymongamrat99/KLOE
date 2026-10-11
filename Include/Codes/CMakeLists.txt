# CMakeLists.txt w katalogu MyProject/Analyses/AnalysisA/

# 1. Definicja nagłówków ze ścieżkami względnymi do katalogu źródłowego
set(LibRec_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/const.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/ErrorLogs.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/ConfigManager.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/AnalysisManager.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/kloe_class.h
    ${CMAKE_CURRENT_SOURCE_DIR}/inc/reconstructor.h
)

# 2. Wygenerowanie słownika ROOT z pełnymi ścieżkami
root_generate_dictionary(LibRecDict 
    ${LibRec_HEADERS}
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/inc/LinkDef.h
)

# Deklarujemy bibliotekę dla Analizy A (dodając wygenerowany słownik LibRecDict.cxx)
add_library(LibRec SHARED 
            src/charged_mom.cpp
            src/chi2_dist.cpp
            src/clear_variables.cpp
            src/closest_approach.cpp
            src/ConfigWatcher.cpp
            src/ConstraintsOmega.cpp
            src/ConstraintsTrilateration.cpp
            src/ConstraintsSignal.cpp
            src/ConstraintsTest.cpp
            src/ConstraintsPM.cpp
            src/ConstraintsNeutral.cpp
            src/CustomGraph.cpp
            src/cylinder_intersection.cpp
            src/DoubleRatio.cpp
            src/EfficiencyAnalysis.cpp
            src/fraction_fitter.cpp
            src/interf_function.cpp
            src/interference.cpp
            src/KinFitter.cpp
            src/kloe_class.cpp
            src/MeasQualityGraph.cpp
            src/MomentumSmearing.cpp
            src/neu_triangle.cpp
            src/neutral_mom.cpp
            src/pi0_photon_pair.cpp
            src/plane_intersection.cpp
            src/reconstructor.cpp
            src/SplitFileWriter.cpp
            src/triple_gaus.cpp
            src/uncertainties.cpp
            src/ConfigManager.cpp
            src/GeneratedVariables.cpp
            src/FileManager.cpp
            src/StatisticalCutter.cpp
            src/StatisticalCutter_Report.cpp
            src/HistManager.cpp
            src/TripleGaussFitter.cpp
            src/DoublGaussFitter.cpp
            src/double_gaus.cpp
            src/KaonMassReconstructor.cpp
            src/trilaterationKinFit.cpp
            src/signalKinFit.cpp
            src/omegaKinFit.cpp
            src/pmKinFit.cpp
            src/neutralKinFit.cpp
            src/NeutralReconstruction.cpp
            src/VariableConfig.cpp
            src/DataAccessWrapper.cpp
            src/const.cpp
            src/AnalysisManager.cpp
            src/PdgApiClient.cpp
            src/PdgManager.cpp
            LibRecDict.cxx # <--- TUTAJ DODANO SŁOWNIK ROOT
            )

# Poinformuj, gdzie są nagłówki dla tej biblioteki.
target_include_directories(LibRec PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
    ${CMAKE_CURRENT_SOURCE_DIR}/../FortranAnalysis/inc
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${Boost_INCLUDE_DIRS}
    ${CURL_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}/nlohmann   # <--- Zapewnia dostęp do pliku json.hpp bezpośrednio
    ${ROOT_INCLUDE_DIRS}          # <--- Zapewnia dostęp do nlohmann/json.hpp
)

# Linkowanie bibliotek
target_link_libraries(LibRec PRIVATE
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
    ${CURL_LIBRARIES}
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)