/* ============================================================================
 * I B E X - Test of the Linearizer of Affine2 forms
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENCE.
 *
 * Author(s)   : Jordan Ninin
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */

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
