# source files of libibex-affine
list (APPEND SRC ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine2_fAF2.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine2_fAF2.h
                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine2.h
                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine3_fAFFullI.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine3_fAFFullI.h
                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineDomain.h
                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine.h
                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineMatrix.h
                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineVector.h
                ${CMAKE_CURRENT_SOURCE_DIR}/function/ibex_AffineEval.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/function/ibex_AffineEval.h
                ${CMAKE_CURRENT_SOURCE_DIR}/numeric/ibex_LinearizerAffine2.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/numeric/ibex_LinearizerAffine2.h
                )

# Create the target for libibex-affine
add_library (ibex-affine ${SRC})
target_include_directories (ibex-affine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                          ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic
                                          ${CMAKE_CURRENT_SOURCE_DIR}/function
                                          ${CMAKE_CURRENT_SOURCE_DIR}/numeric)
target_link_libraries (ibex-affine PUBLIC Ibex::ibex)

# installation of libibex-affine files
list_filter_header (HDR ${SRC}) # Ibex should have installed this function
install (TARGETS ibex-affine DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (FILES ${HDR} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ibex-affine)