cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
## Use -std=c++14 as default.
set(CMAKE_CXX_STANDARD 14)
## Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
## Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Test_AES_example
        LANGUAGES CXX)

find_package(helib REQUIRED)
# add_executable(TEST_AES simpleAES.cpp homAES.cpp Test_AES.cpp)
# target_link_libraries(TEST_AES PUBLIC helib)

add_executable(IDH_CLIENT1 simpleAES.cpp homAES.cpp IDH_client1.cpp)
target_link_libraries(IDH_CLIENT1 PUBLIC helib)

add_executable(IDH_CLIENT2 simpleAES.cpp homAES.cpp IDH_client2.cpp)
target_link_libraries(IDH_CLIENT2 PUBLIC helib)

add_executable(IDVP_CLIENT simpleAES.cpp homAES.cpp IDVP_client.cpp)
target_link_libraries(IDVP_CLIENT PUBLIC helib)
