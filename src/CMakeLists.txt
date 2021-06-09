# source files of libibex-affine
list (APPEND SRC 	${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine2_fAF2.cpp
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine2_fAF2.h
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine3_fAFFullI.cpp
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine3_fAFFullI.h
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineMain.h
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineDomain.h
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineVar.cpp
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineVar.h
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineMatrix.h
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_AffineVector.h
	                ${CMAKE_CURRENT_SOURCE_DIR}/arithmetic/ibex_Affine.h
	                ${CMAKE_CURRENT_SOURCE_DIR}/function/ibex_AffineEval.h
	                ${CMAKE_CURRENT_SOURCE_DIR}/numeric/ibex_LinearizerAffine2.cpp
	                ${CMAKE_CURRENT_SOURCE_DIR}/numeric/ibex_LinearizerAffine2.h
                )

# Create the target for libibex-affine
add_library (ibex-affine ${SRC})
target_include_directories (ibex-affine PUBLIC
                "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/arithmetic>"
                "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/function>"
                "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/numeric>"
                "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ibex-affine>")
target_link_libraries (ibex-affine PUBLIC Ibex::ibex)

# installation of libibex-affine files
ibex_list_filter_header (HDR ${SRC}) # Ibex should have installed this function
install (TARGETS ibex-affine DESTINATION ${CMAKE_INSTALL_LIBDIR}
                              EXPORT AffineExport)
install (FILES ${HDR} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ibex-affine)
