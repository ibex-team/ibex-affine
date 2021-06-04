/* ============================================================================
 * I B E X - Affine evaluator Test
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENSE
 *
 * Author(s)   : Jordan NININ
 * Created     : Juin 20, 2021
 * ---------------------------------------------------------------------------- */


#ifndef __TEST_AFFINE2_H__
#define __TEST_AFFINE2_H__

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "ibex_AffineEval.h"
#include "ibex.h"
#include "utils.h"

using namespace ibex;

template<class T>
class TestAffineEval : public CppUnit::TestFixture {
public:
	CPPUNIT_TEST_SUITE(TestAffineEval<T>);
	CPPUNIT_TEST(test_root2);
	CPPUNIT_TEST(test_root4);
	CPPUNIT_TEST(test_root5);
	CPPUNIT_TEST(test_pow2);
	CPPUNIT_TEST(test_pow5);
	CPPUNIT_TEST(test_powINT1);
	CPPUNIT_TEST(test_powINT2);

	CPPUNIT_TEST(test_sqrt);
	CPPUNIT_TEST(test_exp);
	CPPUNIT_TEST(test_log);
	CPPUNIT_TEST(test_inv);
	CPPUNIT_TEST(test_cos);
	CPPUNIT_TEST(test_sin);
	CPPUNIT_TEST(test_tan);
	CPPUNIT_TEST(test_abs);
	CPPUNIT_TEST(test_acos);
	CPPUNIT_TEST(test_asin);
	CPPUNIT_TEST(test_atan);
	CPPUNIT_TEST(test_cosh);
	CPPUNIT_TEST(test_sinh);
	CPPUNIT_TEST(test_tanh);
	CPPUNIT_TEST_SUITE_END();


	bool check_af2 (Function& f, IntervalVector& I);
	bool check_af2 (Function& f, Interval& I);

	void test_root2();
	void test_root4();
	void test_root5();
	void test_pow2();
	void test_pow4();
	void test_pow5();
	void test_powINT1();
	void test_powINT2();


	void test_sqrt();
	void test_exp();
	void test_log();
	void test_inv();
	void test_cos();
	void test_sin();
	void test_tan();
	void test_abs();
	void test_acos();
	void test_asin();
	void test_atan();
	void test_cosh();
	void test_sinh();
	void test_tanh();


	void test01();
	void test02();

};

CPPUNIT_TEST_SUITE_REGISTRATION(TestAffineEval<AF_Default>);
CPPUNIT_TEST_SUITE_REGISTRATION(TestAffineEval<AF_Other>);



#ifdef _IBEX_WITH_AFFINE_EXTENDED_

CPPUNIT_TEST_SUITE_REGISTRATION(TestAffineEval<AF_fAF2>);
CPPUNIT_TEST_SUITE_REGISTRATION(TestAffineEval<AF_fAFFullI>);
CPPUNIT_TEST_SUITE_REGISTRATION(TestAffineEval<AF_iAF>);
CPPUNIT_TEST_SUITE_REGISTRATION(TestAffineEval<AF_fAF1>);
CPPUNIT_TEST_SUITE_REGISTRATION(TestAffineEval<AF_No>);
CPPUNIT_TEST_SUITE_REGISTRATION(TestAffineEval<AF_sAF>);

#endif

//#include "TestAffineEval.cpp"

#endif
