# CMakeLists.txt w katalogu MyProject/Analyses/AnalysisA/

# Deklarujemy bibliotekę dla Analizy A
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
            src/ConfigManager.cpp
            src/trilaterationKinFit.cpp
            src/signalKinFit.cpp
            src/omegaKinFit.cpp
            src/NeutralReconstruction.cpp
            src/VariableConfig.cpp
            src/DataAccessWrapper.cpp
            src/const.cpp
            src/AnalysisManager.cpp
            src/PdgApiClient.cpp
            src/PdgManager.cpp
            )

# Poinformuj, gdzie są nagłówki dla tej biblioteki.
# PUBLIC: oznacza, że katalog z nagłówkami (tu: bieżący katalog, '.') będzie
# również widoczny dla celów, które linkują do AnalysisALib.

target_include_directories(LibRec PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
    ${CMAKE_CURRENT_SOURCE_DIR}/../FortranAnalysis/inc
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${Boost_INCLUDE_DIRS}
    ${CURL_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIRS}/nlohmann
    # Jeśli AnalysisA.cpp używa nagłówków z src/ lub Subanalysis/,
    # trzeba je tu również dodać (np. parent_directory/src)
    # Przykład: target_include_directories(AnalysisALib PRIVATE ${CMAKE_SOURCE_DIR}/src)
)

# Jeśli AnalysisALib samo potrzebuje bibliotek ROOT, musisz je tu linkować.
# Zmienna ${ROOT_LIBRARIES} jest dziedziczona z głównego CMakeLists.txt.
target_link_libraries(LibRec PRIVATE
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
    # ✅ DODAJ OpenMP do LibRec (gdzie jest KinFitter)
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# # ✅ DODAJ definicję _OPENMP jeśli OpenMP jest dostępny
# if(OpenMP_CXX_FOUND)
#     target_compile_definitions(LibRec PRIVATE _OPENMP)
# endif()

# Możesz tu także zdefiniować lokalny program wykonawczy, np. do testów.
# add_executable(run_analysis_A AnalysisA.cpp)
# target_link_libraries(run_analysis_A PRIVATE AnalysisALib ${ROOT_LIBRARIES})