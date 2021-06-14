/* ============================================================================
 * I B E X - Test of the Forward-Backward with Affine2 forms
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENCE.
 *
 * Author(s)   : Jordan Ninin
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */

#ifndef __TEST_CTC_AFFINE_FWD_BWD_H__
#define __TEST_CTC_AFFINE_FWD_BWD_H__

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "utils.h"
#include "ibex.h"


#ifndef SRCDIR_TESTS
  #define SRCDIR_TESTS "../../tests"
#endif

namespace ibex {

class TestCtcAffineFwdBwd : public CppUnit::TestFixture {

public:

	CPPUNIT_TEST_SUITE(TestCtcAffineFwdBwd);
	CPPUNIT_TEST(sqrt_issue28);
	CPPUNIT_TEST(atan2_issue134);
	CPPUNIT_TEST(solve);
	CPPUNIT_TEST(func02);
	CPPUNIT_TEST(ponts30);
	CPPUNIT_TEST_SUITE_END();

	void sqrt_issue28();
	void atan2_issue134();
	void ponts30();
	void func02();
	void solve();

};

CPPUNIT_TEST_SUITE_REGISTRATION(TestCtcAffineFwdBwd);


} // namespace ibex

#endif // __TEST_CTC_AFFINE_FWD_BWD_H__
