# CMakeLists.txt dla przykładów aplikacji desktopowych

cmake_minimum_required(VERSION 3.13)

# Demo aplikacji (zawsze kompiluje się)
add_executable(desktop_app_demo desktop_app_demo.cpp)
target_compile_features(desktop_app_demo PRIVATE cxx_std_14)

# Podstawowy przykład (zawsze kompiluje się)
add_executable(desktop_app_basic desktop_app_basic.cpp)
target_compile_features(desktop_app_basic PRIVATE cxx_std_14)

# Przykład Qt (tylko jeśli Qt jest dostępne)
find_package(Qt5 QUIET COMPONENTS Core Widgets)
if(Qt5_FOUND)
    message(STATUS "Qt5 znalezione - kompilowanie przykładu Qt")
    
    add_executable(desktop_app_qt desktop_app_qt.cpp)
    target_compile_definitions(desktop_app_qt PRIVATE QT_AVAILABLE)
    target_link_libraries(desktop_app_qt Qt5::Core Qt5::Widgets)
    
    # Automatyczne przetwarzanie MOC dla Qt
    set_target_properties(desktop_app_qt PROPERTIES AUTOMOC ON)
else()
    message(STATUS "Qt5 nie znalezione - przykład Qt będzie podstawowy")
    add_executable(desktop_app_qt desktop_app_qt.cpp)
endif()

# Przykład ROOT (tylko jeśli ROOT jest dostępne)
find_package(ROOT QUIET COMPONENTS Core Gpad Graf Gui)
if(ROOT_FOUND)
    message(STATUS "ROOT znalezione - kompilowanie przykładu ROOT GUI")
    
    add_executable(desktop_app_root desktop_app_root.cpp)
    target_compile_definitions(desktop_app_root PRIVATE ROOT_AVAILABLE)
    target_include_directories(desktop_app_root PRIVATE ${ROOT_INCLUDE_DIRS})
    target_link_libraries(desktop_app_root ${ROOT_LIBRARIES})
    
    # ROOT może wymagać specjalnych flag kompilacji
    target_compile_features(desktop_app_root PRIVATE cxx_std_14)
else()
    message(STATUS "ROOT nie znalezione - przykład ROOT będzie podstawowy")
    add_executable(desktop_app_root desktop_app_root.cpp)
endif()

# Instalacja przykładów
install(TARGETS desktop_app_demo desktop_app_basic desktop_app_qt desktop_app_root
        RUNTIME DESTINATION bin/examples)

# Informacje o kompilacji
message(STATUS "=== Kompilowanie przykładów aplikacji desktopowych ===")
message(STATUS "Qt5 dostępne: ${Qt5_FOUND}")
message(STATUS "ROOT dostępne: ${ROOT_FOUND}")
message(STATUS "Wszystkie przykłady będą skompilowane, ale funkcjonalność zależy od dostępnych bibliotek.")