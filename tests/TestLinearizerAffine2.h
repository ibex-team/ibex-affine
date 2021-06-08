//============================================================================
//                                  I B E X                                   
// File        : TestLinearizerAffine2.h
// Author      : Gilles Chabert and Jordan NININ
// Copyright   : ENSTA Bretagne (France)
// License     : See the LICENSE file
// Created     : June 6, 20120
//============================================================================

#ifndef __TEST_LINEARIZER_AFFINE2_H__
#define __TEST_LINEARIZER_AFFINE2_H__

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "utils.h"
#include "ibex_LPSolver.h"
#include "ibex_LinearizerAffine2.h"

namespace ibex {

class TestLinearizerAffine2 : public CppUnit::TestFixture {
public:

	CPPUNIT_TEST_SUITE(TestLinearizerAffine2);
	

#ifndef __IBEX_NO_LP_SOLVER__

		CPPUNIT_TEST(fixbug01);

#endif //__IBEX_NO_LP_SOLVER__

	CPPUNIT_TEST_SUITE_END();


	void fixbug01();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestLinearizerAffine2);


} // end namespace ibex
#endif // __TEST_CTC_POLYTOPE_HULL_H__
