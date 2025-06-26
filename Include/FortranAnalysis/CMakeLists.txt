# CMakeLists.txt w katalogu MyProject/Analyses/AnalysisA/

# Deklarujemy bibliotekę dla Analizy A
add_library(LibFort SHARED 
            src/analysismodule.f
            src/analysis_library.f)

target_include_directories(LibFort PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/inc
)

target_link_libraries(LibFort PUBLIC 
    mathlib
    kernlib
    packlib
    pawlib
    graflib
    kernlib)