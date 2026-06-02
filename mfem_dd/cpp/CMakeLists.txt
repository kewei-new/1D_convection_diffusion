cmake_minimum_required(VERSION 3.20)
project(mfem_dd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(DEFINED MFEM_ROOT AND NOT DEFINED MFEM_DIR)
  list(APPEND CMAKE_PREFIX_PATH
       "${MFEM_ROOT}"
       "${MFEM_ROOT}/build"
       "${MFEM_ROOT}/lib/cmake/mfem")
endif()

list(APPEND CMAKE_PREFIX_PATH
     "${CMAKE_CURRENT_LIST_DIR}/../../../mfem-4.9/build"
     "${CMAKE_CURRENT_LIST_DIR}/../../mfem-4.9/build")

find_package(MFEM CONFIG REQUIRED)

if(NOT TARGET mfem)
  message(FATAL_ERROR "MFEMConfig.cmake was found, but target 'mfem' is not available.")
endif()

function(add_mfem_dd_app name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE mfem)
  target_compile_features(${name} PRIVATE cxx_std_17)
  file(TO_CMAKE_PATH "${CMAKE_CURRENT_LIST_DIR}/../matlab" MFEM_DD_MATLAB_DIR)
  target_compile_definitions(${name} PRIVATE
    "MFEM_DD_MATLAB_DIR=\"${MFEM_DD_MATLAB_DIR}\"")
endfunction()

add_mfem_dd_app(dd1d_mms apps/dd1d_mms.cpp)
add_mfem_dd_app(dd2d_mms apps/dd2d_mms.cpp)
add_mfem_dd_app(dd_device apps/dd_device.cpp)

enable_testing()
add_test(NAME dd1d_mms_smoke COMMAND dd1d_mms -n 8 -o 2)
add_test(NAME dd1d_mms_legacy_baseline COMMAND dd1d_mms -n 20 -o 3 -b legacy_baseline)
add_test(NAME dd1d_mms_native_mfem COMMAND dd1d_mms -n 20 -o 3 -b native_mfem)
add_test(NAME dd2d_mms_smoke COMMAND dd2d_mms -n 4 -o 1)
add_test(NAME dd_device_smoke COMMAND dd_device -n 8 -o 1)
add_test(NAME dd_device_legacy_baseline COMMAND dd_device -n 16 -o 1 -b legacy_baseline)
add_test(NAME dd_device_native_mfem COMMAND dd_device -n 16 -o 1 -b native_mfem)
add_test(NAME dd_device_native_operator_snapshot COMMAND dd_device -n 160 -o 2 -b native_operator_snapshot)
