cmake_minimum_required(VERSION 3.10)
project(MN13 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Biblioteka statyczna ---
add_library(mynum STATIC source.cpp)
target_include_directories(mynum PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# --- Program main.cpp ---
add_executable(main main.cpp)
target_link_libraries(main mynum)

# --- Catch2 testy ---
enable_testing()
add_executable(tests tests.cpp catch.hpp)
target_link_libraries(tests mynum)
add_test(NAME all COMMAND tests)

