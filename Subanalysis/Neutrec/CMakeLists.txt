# CMakeLists.txt w katalogu MyProject/Analyses/AnalysisA/

# Deklarujemy bibliotekę dla Analizy A
add_library(NeutRec SHARED 
            app/neutrec_main.cpp
            src/triangle_neu_rec.cpp
            src/comp_of_methods.cpp
            src/trilateration_neu_rec_kin_fit_aux.cpp
            src/trilaterationKinFit.cpp)

# Poinformuj, gdzie są nagłówki dla tej biblioteki.
# PUBLIC: oznacza, że katalog z nagłówkami (tu: bieżący katalog, '.') będzie
# również widoczny dla celów, które linkują do AnalysisALib.
target_include_directories(NeutRec PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
    # Jeśli AnalysisA.cpp używa nagłówków z src/ lub Subanalysis/,
    # trzeba je tu również dodać (np. parent_directory/src)
    # Przykład: target_include_directories(AnalysisALib PRIVATE ${CMAKE_SOURCE_DIR}/src)
)

# Jeśli AnalysisALib samo potrzebuje bibliotek ROOT, musisz je tu linkować.
# Zmienna ${ROOT_LIBRARIES} jest dziedziczona z głównego CMakeLists.txt.
target_link_libraries(NeutRec PRIVATE
    ${ROOT_LIBRARIES}
    LibRec # Jeśli AnalysisALib zależy od SubanalysisLib
)

# Możesz tu także zdefiniować lokalny program wykonawczy, np. do testów.
# add_executable(run_analysis_A AnalysisA.cpp)
# target_link_libraries(run_analysis_A PRIVATE AnalysisALib ${ROOT_LIBRARIES})