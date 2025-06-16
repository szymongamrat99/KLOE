# MyProject/tests/CMakeLists.txt

add_executable(LibRecTest LibRec_test.cpp)

# Jeśli Catch2 został pobrany jako pojedynczy plik nagłówkowy (catch.hpp):
# Musisz dodać ścieżkę do katalogu, w którym znajduje się ten plik nagłówkowy.
# FetchContent_Populate ustawia zmienną ${Catch2_SOURCE_DIR} na ścieżkę do pobranego źródła.
target_include_directories(LibRecTest PRIVATE
    ${CMAKE_SOURCE_DIR}/../../Include/Codes/inc # Zapewnia dostęp do MyUtility.h
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${CURL_INCLUDE_DIRS}
)

# Linkowanie bibliotek
# W przypadku Catch2 jako pliku header-only, nie linkujesz biblioteki Catch2.
# Dyrektywa #define CATCH_CONFIG_MAIN w pliku testowym tworzy funkcję main()
# i kompiluje kod Catch2 bezpośrednio w Twój plik wykonywalny testu.
target_link_libraries(LibRecTest PRIVATE
    # Jeśli testujesz RootReaderLib:
    # RootReaderLib
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
    ${CURL_LIBRARIES}
    LibRec
    # ... inne biblioteki
)

# Zarejestruj program testowy w CTest
add_test(NAME LibRecTest COMMAND LibRecTest)