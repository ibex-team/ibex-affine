/* ============================================================================
 * I B E X - Affine Evaluator Test
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENSE
 *
 * Author(s)   : Jordan NININ
 * Created     : Juin 20, 2021
 * ---------------------------------------------------------------------------- */

#include "TestAffineEval.h"

#include "ibex_Function.h"
#include "ibex_CtcFwdBwd.h"
#include "ibex_CtcFixPoint.h"
#include "ibex_LargestFirst.h"
#include "ibex_CellStack.h"
#include "ibex_Solver.h"

#include "ibex_AffineEval.h"
//using namespace std;

template<class T>
void TestAffineEval<T>::test01() {
	const ExprSymbol& x = ExprSymbol::new_("x",Dim::col_vec(2));
	Function f(x,x[0]*pow(x[1],2)+exp(x[1]*x[0]));
	IntervalVector itv(2,Interval(1,2));
	CPPUNIT_ASSERT(check_af2(f,itv));

}

template<class T>
void TestAffineEval<T>::test02() {
	const ExprSymbol& x = ExprSymbol::new_();
	Function f(x,cosh(x)-x);
	IntervalVector itv(1,Interval(1,2));
	CPPUNIT_ASSERT(check_af2(f,itv));

}

template<class T>
void TestAffineEval<T>::test_pow2() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,pow(x,2));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}

template<class T>
void TestAffineEval<T>::test_pow4() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,pow(x,4));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}



template<class T>
void TestAffineEval<T>::test_pow5() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,pow(x,5));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}


template<class T>
void TestAffineEval<T>::test_root2() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,pow(x,1/2));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}

template<class T>
void TestAffineEval<T>::test_root4() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,pow(x,1/4));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}



template<class T>
void TestAffineEval<T>::test_root5() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,pow(x,1/5));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}



template<class T>
void TestAffineEval<T>::test_powINT1() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,pow(x,Interval(2,3)));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}

template<class T>
void TestAffineEval<T>::test_powINT2() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,pow(x,Interval(-2,3)));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}


template<class T>
void TestAffineEval<T>::test_sqrt() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,sqrt(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}

template<class T>
void TestAffineEval<T>::test_exp() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,exp(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-10,50000000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-50000000,1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.00001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1000000000.,1000000000.001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}

template<class T>
void TestAffineEval<T>::test_log() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,log(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.00001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1000000000.,1000000000.001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10,200000000);
	CPPUNIT_ASSERT(check_af2(f,itv));
}

template<class T>
void TestAffineEval<T>::test_inv() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,1.0/x);
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10,20);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.00001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10000000,10000000.0001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-20000000,-1);
	CPPUNIT_ASSERT(check_af2(f,itv));
}

template<class T>
void TestAffineEval<T>::test_cos() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,cos(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.00001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.0000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));


	itv =Interval::PI;
	CPPUNIT_ASSERT(check_af2(f,itv));
	for (int i =-20;i<34;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<15;i++) {
		itv =Interval(i*Interval::PI.lb()/6,(i+1)*Interval::PI.ub()/6);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+3)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+2)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+3)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<30;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+7)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<34;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}

}

template<class T>
void TestAffineEval<T>::test_sin() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,sin(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.00001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));


	itv =Interval::PI;
	CPPUNIT_ASSERT(check_af2(f,itv));
	for (int i =-20;i<34;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<15;i++) {
		itv =Interval(i*Interval::PI.lb()/6,(i+1)*Interval::PI.ub()/6);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+3)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+2)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+3)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<30;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+7)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<34;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}


}

template<class T>
void TestAffineEval<T>::test_tan() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,tan(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));

	itv =Interval::PI;
	CPPUNIT_ASSERT(check_af2(f,itv));
	for (int i =-20;i<34;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<15;i++) {
		itv =Interval(i*Interval::PI.lb()/6,(i+1)*Interval::PI.ub()/6);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+3)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+2)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+3)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<30;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+7)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<34;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}

}

template<class T>
void TestAffineEval<T>::test_abs() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,abs(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10000000,10000000.0001);
}

template<class T>
void TestAffineEval<T>::test_acos() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,acos(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-10,20);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-10,0.2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(0.5,0.7);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10000000,10000000.0001);

	itv =Interval::PI;
	CPPUNIT_ASSERT(check_af2(f,itv));
	for (int i =-20;i<34;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<15;i++) {
		itv =Interval(i*Interval::PI.lb()/6,(i+1)*Interval::PI.ub()/6);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+3)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+2)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+3)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<30;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+7)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<34;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}

}

template<class T>
void TestAffineEval<T>::test_asin() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,asin(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-10,20);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-10,0.2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(0.5,0.7);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10000000,10000000.0001);

	itv =Interval::PI;
	CPPUNIT_ASSERT(check_af2(f,itv));
	for (int i =-20;i<34;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<15;i++) {
		itv =Interval(i*Interval::PI.lb()/6,(i+1)*Interval::PI.ub()/6);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+3)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+2)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+3)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<30;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+7)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<34;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}

}

template<class T>
void TestAffineEval<T>::test_atan() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,atan(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10000000,10000000.0001);

	itv =Interval::PI;
	CPPUNIT_ASSERT(check_af2(f,itv));
	for (int i =-20;i<34;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<15;i++) {
		itv =Interval(i*Interval::PI.lb()/6,(i+1)*Interval::PI.ub()/6);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+3)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+2)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<20;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+3)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<40;i++) {
		itv =Interval(i*Interval::PI.lb()/16,(i+4)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<30;i++) {
		itv =Interval(i*Interval::PI.lb()/8,(i+7)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<34;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/16,(i+1)*Interval::PI.ub()/16);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<18;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/8,(i+1)*Interval::PI.ub()/8);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}
	for (int i =-20;i<10;i++) {
		itv =Interval(0.1+i*Interval::PI.lb()/4,(i+1)*Interval::PI.ub()/4);
		CPPUNIT_ASSERT(check_af2(f,itv));
	}

}

template<class T>
void TestAffineEval<T>::test_cosh() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,cosh(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10000000,10000000.0001);
}

template<class T>
void TestAffineEval<T>::test_sinh() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,sinh(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10000000,10000000.0001);
}



template<class T>
void TestAffineEval<T>::test_tanh() {
	const ExprSymbol& x = ExprSymbol::new_();
	Interval itv;
	Function f(x,tanh(x));
	itv =Interval(1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,POS_INFINITY);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,-2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(NEG_INFINITY,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::ALL_REALS;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval::EMPTY_SET;
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,2);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(-1,-0.5);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,200000);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1,1.000000000001);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(1);
	CPPUNIT_ASSERT(check_af2(f,itv));
	itv =Interval(10000000,10000000.0001);
}



template<class T>
bool TestAffineEval<T>::check_af2 (Function& f, IntervalVector& I){
	return check_af2(f,I[0]);
}

template<class T>
bool TestAffineEval<T>::check_af2 (Function& f, Interval& I){
	double n = 100; // number of try
	AffineMain<T> faa;
	Interval itv2;
	Interval itv;
	AffineMain<T> aitv ;

	AffineEval<T> eval_af(f);

	if (!((I.is_unbounded())||(I.is_degenerated())||(!I.is_bisectable())||(I.is_empty()))){
		for (double ii= I.lb(); ii<I.ub(); ii += I.diam()/n) {
			itv2 =f.eval(IntervalVector(1,Interval(ii)));

			itv = eval_af.eval(IntervalVector(1,Interval(ii))).i();
			faa = eval_af.af2.top->i();
			aitv = eval_af.eval(AffineVarMainVector<T>(1,Interval(ii))).i();


			if (!(itv2.is_subset(faa.itv())))
			{
				std::cout  << " DEP = "<< ii<< "  "  << f<< std::endl;
				std::cout  << " RES = "<< itv2 << " /// "<< itv << " ///// " << faa << std::endl;
				std::cout  << " RES = "<< itv2 << " ///// " << faa << std::endl;

				return false;
			}
			if (!(itv2.is_subset(aitv.itv())))
			{
				std::cout  << " DEP = "<< ii<< "  "  << f<< std::endl;
				std::cout  << " RES = "<< itv2 << " /// "<< itv << " ///// " << aitv << std::endl;
				std::cout  << " RES = "<< itv2 << " ///// " << aitv << std::endl;

				return false;
			}
		}
	}


	itv2 =f.eval(IntervalVector(1,I));

	itv = eval_af.eval(IntervalVector(1,I)).i();
	faa = eval_af.af2.top->i();

	aitv = eval_af.eval(AffineVarMainVector<T>(1,I)).i();

	if (!(itv2.is_subset(faa.itv())))
	{
		std::cout  << " DEP = "<< I<< "  "  << f<< std::endl;
		std::cout  << " RES = "<< itv2 << " /// "<< itv << " ///// " << faa << std::endl;
		std::cout  << " RES = "<< itv2 << " ///// " << faa << std::endl;

		return false;
	}
	if (!(itv2.is_subset(aitv.itv())))
	{
		std::cout  << " DEP = "<< I<< "  "  << f<< std::endl;
		std::cout  << " RES = "<< itv2 << " /// "<< itv << " ///// " << aitv << std::endl;
		std::cout  << " RES = "<< itv2 << " ///// " << aitv << std::endl;

		return false;
	}
/*	if (faa.size()>0) {
		const ExprSymbol& x = ExprSymbol::new_();
		Function lininf(x,(faa.val(0)+faa.val(1)*(2* (x)-(I.lb()+I.ub()))/(I.diam())));

		Function linsup(x, lininf(x)+faa.err().lb());

		Function c_inf(x,lininf(x)+faa.err().lb()  -f(x));
		Function c_sup(x, f(x) -linsup(x));

		CtcFwdBwd ct1(c_inf,GT,AFFINE2_MODE);
		CtcFixPoint ft1(ct1,0.1);

		CtcFwdBwd ct2(c_sup,GT,AFFINE2_MODE);
		CtcFixPoint ft2(ct2,0.1);

		LargestFirst bbb;
		CellStack ccc;
		Solver sol1(ft1,bbb, ccc, 0.1 );
		Solver sol2(ft2,bbb, ccc, 0.1 );


		std::cout  << " DEP = "<< I<< "  "  << f<< std::endl;
		std::cout  << " RES = "<< itv2 << " /// "<< itv << " ///// " << faa << std::endl;
		std::cout<< " PAVING = "<< I.diam()<< " /// " << sol1.solve(IntervalVector(1,I)).size() << " et "<<sol2.solve(IntervalVector(1,I)).size() << std::endl;
//		return ((sol1.solve(IntervalVector(1,I)).empty()) && (sol2.solve(IntervalVector(1,I)).empty()));
	}
*/
	return true;

}

template<class T>
void TestAffineEval<T>::issue242() {
	Function f("x[3]","-x");
	IntervalVector x(3,Interval::one());
	AffineEval<T> evalf(f);
	IntervalVector res = evalf.eval(x).v();
	CPPUNIT_ASSERT(almost_eq(res,-x,0));

	AffineMainVector<T> resaf = evalf.eval(AffineVarMainVector<T>(x)).v();
	CPPUNIT_ASSERT(almost_eq(resaf.itv(),-x,0));
}
template<class T>
void TestAffineEval<T>::eval_components01() {
	const ExprSymbol& x = ExprSymbol::new_("x");
	const ExprSymbol& y = ExprSymbol::new_("y");
	const ExprSymbol& z = ExprSymbol::new_("z");
	const ExprNode& e1=x+3*y;
	const ExprNode& e2=y-2*x;
	const ExprNode& e3=e1*e2;
	Function f(x,y,z,Return(e3+1,e2+3,e3-2,e3-4));

	Interval vx=Interval::one();
	Interval vy=2*Interval::one();
	IntervalVector box(3);
	box[0]=vx;
	box[1]=vy;

	BitSet components=BitSet::empty(4);
	components.add(0);
	components.add(2);


	AffineEval<T> evalf(f);
	IntervalVector res = evalf.eval(box,components).v();
	CPPUNIT_ASSERT(res.size()==2);
	CPPUNIT_ASSERT(res[0]==(vx+3*vy)*(vy-2*vx)+1);
	CPPUNIT_ASSERT(res[1]==(vx+3*vy)*(vy-2*vx)-2);

	AffineMainVector<T> resaf = evalf.eval(AffineVarMainVector<T>(box),components).v();

	CPPUNIT_ASSERT(resaf.size()==2);
	CPPUNIT_ASSERT(resaf[0]==(vx+3*vy)*(vy-2*vx)+1);
	CPPUNIT_ASSERT(resaf[1]==(vx+3*vy)*(vy-2*vx)-2);
}
template<class T>
void TestAffineEval<T>::eval_components02() {
	Dim d=Dim::matrix(3,3);
	const ExprSymbol& x = ExprSymbol::new_("x",Dim::col_vec(2));
	const ExprSymbol& y = ExprSymbol::new_("y",d);
	const ExprSymbol& z = ExprSymbol::new_("z",d);

	Function f(x,y,z,Return(x[1],transpose(y[DoubleIndex::one_row(d,1)]),z[DoubleIndex::one_col(d,2)]));

	IntervalVector box(20);
	for (int i=0; i<20; i++) box[i]=Interval(i,i);

	BitSet components=BitSet::empty(9);
	components.add(0);
	components.add(2);
	components.add(4);
	components.add(6);

	AffineEval<T> evalf(f);
	IntervalVector res = evalf.eval(box,components).v();

	CPPUNIT_ASSERT(res.size()==4);
	CPPUNIT_ASSERT(res[0]==1);
	CPPUNIT_ASSERT(res[1]==6);
	CPPUNIT_ASSERT(res[2]==13);
	CPPUNIT_ASSERT(res[3]==19);

	AffineMainVector<T> resaf = evalf.eval(AffineVarMainVector<T>(box),components).v();
	CPPUNIT_ASSERT(resaf.size()==4);
	CPPUNIT_ASSERT(resaf[0]==1);
	CPPUNIT_ASSERT(resaf[1]==6);
	CPPUNIT_ASSERT(resaf[2]==13);
	CPPUNIT_ASSERT(resaf[3]==19);
}
template<class T>
void TestAffineEval<T>::matrix_components() {
	const ExprSymbol& x = ExprSymbol::new_();
	Function f(x,Return(
			Return(x  ,x+1,x+2,ExprVector::ROW),
			Return(x+3,x+4,x+5,ExprVector::ROW),
			Return(x+6,x+7,x+8,ExprVector::ROW),
			ExprVector::COL));

	double _M1[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

	AffineEval<T> evalf(f);
	IntervalMatrix res = evalf.eval(IntervalVector(1,Interval::ZERO)).m();
	CPPUNIT_ASSERT(res == Matrix(3,3,_M1));


	// perform another evaluation
	double _M2[] = { 9, 10, 11, 15, 16, 17 };
	BitSet bitset(3);
	bitset.add(0);
	bitset.add(2);

	IntervalMatrix res2 = evalf.eval(IntervalVector(1,Interval(9,9)),bitset).m();
	CPPUNIT_ASSERT(res2 == Matrix(2,3,_M2));

}

