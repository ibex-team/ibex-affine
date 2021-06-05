# Looking for cppunit (the script FindCppUnit.cmake should have been installed
# by Ibex and the CMAKE_MODULE_PATH updated accordingly by ibex_init_common)
find_package (CppUnit)

# Only add the test targets if cppunit was found
if (CPPUNIT_FOUND)
  # Compile common stuff for the tests
  add_library (test_common utest.cpp utils.cpp utils.h)
  target_link_libraries (test_common PUBLIC ibex-affine)

  set (TESTS_LIST TestAffineArith)# TestAffineEval TestAffineBase TestAffineVector )

  foreach (test ${TESTS_LIST})
    # /!\ The test and the target building the executable have the same name
    add_executable (${test} ${test}.cpp ${test}.h)
    target_link_libraries (${test} test_common CppUnit::CppUnit)
    add_dependencies (check ${test})
    add_test (${test} ${test})
  endforeach ()
else ()
  message (STATUS "Will not run tests, required cppunit library was not found")
  set (MSG "No tests will be run as CMake failed to find the cppunit library \
during the configuration step. If you want to run the tests, you need to \
install cppunit. Then, you need to re-run the configuration step of CMake. You \
can use the CPPUNIT_DIR variable to indicate to CMake were to find it")
  add_custom_target (notest COMMAND ${CMAKE_COMMAND} -E echo ${MSG})
  add_dependencies (check notest)
endif ()
