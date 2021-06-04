/* ============================================================================
 * I B E X - Test of the Affine forms
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENSE
 *
 * Author(s)   : Jordan NININ
 * Created     : Juin 20, 2021
 * ---------------------------------------------------------------------------- */

#include "TestAffineArith.h"
#include <string>
#include <stdio.h>
#include <iomanip>

template<class T>
bool TestAffineArith<T>::compare_results (Interval r, AffineMain<T>  a) {
	Interval ra = a.itv();
	if (! a.is_actif()) {
		if (r==ra) {
			return true;
		} else {
			std::cout<<std::endl<<"Affine:  "<<a<<std::endl;
			std::cout<<"AffineI: "<<ra<<std::endl;
			std::cout<<"Interval:"<<r<<std::endl;
			return false;
		}
	} else {
		bool e1 = ra == r;
		bool e2 = (r.is_subset(ra) && ra.is_subset(r.inflate(abs(r).ub()*1.e-14)) );
		bool e3 = fabs(r.mid()-a.mid()) < 1.e-14;
		bool e4 = fabs(r.rad()-a.val(0)) < 1.e-14;
		if (e1 || e2 || (e3 && e4)) {
			return true;
		}
		else {

			std::cout<<std::endl<<"Affine:  "<<a<<std::endl;
			std::cout<<"AffineI: "<<ra<<std::endl;
			std::cout<<"Interval:"<<r<<std::endl;
			return false;
		}
	}
}


template<class T>
bool TestAffineArith<T>::compare_results2 (Interval r, AffineMain<T>  a) {
	Interval ra = a.itv();
	if (! a.is_actif()) {
		return r==ra;
	} else {
		bool e1 = ra == r;
		bool e2 = (r.is_subset(ra) && ra.is_subset(r.inflate(abs(r).ub()*1.e-14)) );
		bool e3 = fabs(r.mid()-a.mid()) < 1.e-14;
		bool e4 = fabs(r.rad()-(a.val(0)+a.val(1))) < 1.e-14;
		if (e1 || e2 || (e3 && e4)) {
			return true;
		}
		else {

			std::cout<<std::endl<<"Affine:  "<<a<<std::endl;
			std::cout<<"AffineI: "<<ra<<std::endl;
			std::cout<<"Interval:"<<r<<std::endl;
			return false;
		}
	}
}





template<class T>
void TestAffineArith<T>::test01(){
	Interval x(-1, 1);
	AffineMainVector<T> ax(1, x);
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test02() {
	Interval x(-1, 2.13);
	AffineMainVector<T> ax(1, x);
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test03()  {
	Interval x(5, 7.43);
	AffineMainVector<T> ax (1,x);
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test04()  {
	Interval x(-15, -7.45);
	AffineMainVector<T> ax (1,x);
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

// Float addition
template<class T>
void TestAffineArith<T>::test05()  {
	Interval x(1.0, 2.0);
	AffineMainVector<T> ax (1,x);
	x = x + 1.0;
	ax[0] = ax[0] + 1.0;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test06()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x + 0.0;
	ax[0] = ax[0] + 0.0;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test07()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x + 3.14;
	ax[0] = ax[0] + 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test08()  {
	Interval x(-5, -2);
	AffineMainVector<T> ax (1,x);
	x = x + 3.14;
	ax[0] = ax[0] + 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test09()  {
	Interval x(-5, 4);
	AffineMainVector<T> ax (1,x);
	x = x + 3.14;
	ax[0] = ax[0] + 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

// Float subtraction

template<class T>
void TestAffineArith<T>::test10()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x - 1.0;
	ax[0] = ax[0] - 1.0;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test11()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x - 0.0;
	ax[0] = ax[0] - 0.0;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test12()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x - 3.14;
	ax[0] = ax[0] - 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test13() {
	Interval x(-5, -2);
	AffineMainVector<T> ax (1,x);
	x = x - 3.14;
	ax[0] = ax[0] - 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test14() {
	Interval x(-5, 4);
	AffineMainVector<T> ax (1,x);
	x = x - 3.14;
	ax[0] = ax[0] - 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}


// Scalar multiplication

template<class T>
void TestAffineArith<T>::test15()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x * 1.0;
	ax[0] = ax[0] * 1.0;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test16()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x * 0.0;
	ax[0] = ax[0] * 0.0;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}


template<class T>
void TestAffineArith<T>::test17()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x * 3.14;
	ax[0] = ax[0] * 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test18()   {
	Interval x(-5, -2);
	AffineMainVector<T> ax (1,x);
	x = x * 3.14;
	ax[0] = ax[0] * 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test19()   {
	Interval x(-5, 4);
	AffineMainVector<T> ax (1,x);
	x = x * 3.14;
	ax[0] = ax[0] * 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

// Float division

template<class T>
void TestAffineArith<T>::test20()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x / 1.0;
	ax[0] = ax[0] / 1.0;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test21()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	x = x / 3.14;
	ax[0] = ax[0] / 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test22()   {
	Interval x(-5, -2);
	AffineMainVector<T> ax (1,x);
	x = x / 3.14;
	ax[0] = ax[0] / 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test23()  {
	Interval x(-5, 4);
	AffineMainVector<T> ax (1,x);
	x = x / 3.14;
	ax[0] = ax[0] / 3.14;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

/* ****************************** */

// Interval addition

template<class T>
void TestAffineArith<T>::test24()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(1.0, 1.0);
	x = x + prod;
	ax[0] = ax[0] + prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test25()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(0.0, 0.0);
	x = x + prod;
	ax[0] = ax[0] + prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test26()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(-1.0, 1.0);
	x = x + prod;
	ax[0] = ax[0] + prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test27()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(-1.7, 1.7);
	x = x + prod;
	ax[0] = ax[0] + prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test28()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(-1.08, 1.7);
	x = x + prod;
	ax[0] = ax[0] + prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}


// Interval multiplication

template<class T>
void TestAffineArith<T>::test29()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(1.0, 1.0);
	x = x * prod;
	ax[0] = ax[0] * prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test30()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(0.0, 0.0);
	x = x * prod;
	ax[0] = ax[0] * prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test31()  {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(-1.0, 1.0);
	x = x * prod;
	ax[0] = ax[0] * prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test32()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(-1.7, 1.7);
	x = x * prod;
	ax[0] = ax[0] * prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test33()   {
	Interval x(1, 2);
	Interval prod(-1.08, 1.7);
	AffineMainVector<T> ax (2);
	std::cout << std::endl << "coucou0: "<< ax[0]<< std::endl;
	ax[0]=x;
	std::cout << std::endl << "coucou0: "<< ax[0]<< std::endl;
	std::cout << std::endl << "coucou1: "<< ax[1]<< std::endl;
	ax[1]=prod;
	std::cout << std::endl << "coucou1: "<< ax[1]<< std::endl;
	AffineMain<T> tmp;
	tmp=prod;
	std::cout << std::endl << "coucoutmp: "<< tmp<< std::endl;
	std::cout << std::endl << "coucou0: "<< ax[0]<< std::endl;
	x = x * prod;
	ax[0] = ax[0] * ax[1];
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test33_2()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(1.08, 1.7);
	x = x * prod;
	ax[0] = ax[0] * prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test33_3()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	Interval prod(-1.08, -0.7);
	x = x * prod;
	ax[0] = ax[0] * prod;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}


/* ********************************************** */
template<class T>
void TestAffineArith<T>::test34()   {
	Interval x(1, 2);
	AffineMainVector<T> ax (1,x);
	AffineMain<T> ay = ax[0];
	CPPUNIT_ASSERT(compare_results ( x, ay));
}


template<class T>
void TestAffineArith<T>::test35()   {
	Interval x(-11.3, -4.3);
	AffineMainVector<T> ax (1,x);
	AffineMain<T> ay = ax[0];
	CPPUNIT_ASSERT(compare_results ( x, ay));
}

template<class T>
void TestAffineArith<T>::test36()   {
	Interval x(-11.3, 4.3);
	AffineMainVector<T> ax (1,x);
	AffineMain<T> ay = ax[0];
	CPPUNIT_ASSERT(compare_results ( x, ay));
}

template<class T>
void TestAffineArith<T>::test37()   {
	Interval x(0, 4.3);
	AffineMainVector<T> ax (1,x);
	AffineMain<T> ay = ax[0];
	CPPUNIT_ASSERT(compare_results ( x, ay));
}

/* ********************************************** */

// Assign-addition with interval
template<class T>
void TestAffineArith<T>::test38()   {
	Interval x(1.57, 2.23);
	AffineMainVector<T> ax (1,x);
	ax[0] += x;
	x += x;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}


template<class T>
void TestAffineArith<T>::test39()  {
	Interval x(-11.3, -4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] += x;
	x += x;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test49()  {
	Interval x(-11.3, 4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] += x;
	x += x;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test50()  {
	Interval x(0, 4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] += x;
	x += x;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

// Assign-multiplication with interval
template<class T>
void TestAffineArith<T>::test51()   {
	Interval x(1.57, 2.23);
	AffineMainVector<T> ax (1,x);
	ax[0] *= x;
	x *= x;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}


template<class T>
void TestAffineArith<T>::test52()  {
	Interval x(-11.3, -4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] *= x;
	x *= x;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test53()  {
	Interval x(-11.3, 4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] *= x;
	x *= x;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test54()  {
	Interval x(0, 4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] *= x;
	x *= x;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

/* ********************************************* */


// Assign-addition with float
template<class T>
void TestAffineArith<T>::test55()  {
	Interval x(1,2);
	AffineMainVector<T> ax (1,x);
	ax[0] += 3.14159;
	x += 3.14159;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}


template<class T>
void TestAffineArith<T>::test56()  {
	Interval x(-11.3, -4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] += 3.14159;
	x += 3.14159;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test57()   {
	Interval x(-11.3, 4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] += 3.14159;
	x += 3.14159;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test58()  {
	Interval x(0, 4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] += 3.14159;
	x += 3.14159;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test59()  {
	Interval x(0, 4.31);
	AffineMainVector<T> ax (1,x);
	ax[0] += 0.;
	x += 0.;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

// Assign-multiplication with float
template<class T>
void TestAffineArith<T>::test60()  {
	Interval x(1.57, 2.23);
	AffineMainVector<T> ax (1,x);
	ax[0] *= 1.;
	x *= 1.;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test61()  {
	Interval x(1.57, 2.23);
	AffineMainVector<T> ax (1,x);
	ax[0] *= 0.;
	x *= 0.;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test62()  {
	Interval x(1.57, 2.23);
	AffineMainVector<T> ax (1,x);
	ax[0] *= 3.14159;
	x *= 3.14159;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}


template<class T>
void TestAffineArith<T>::test63()  {
	Interval x(-11.3, -4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] *= 3.14159;
	x *= 3.14159;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test64()   {
	Interval x(-11.3, 4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] *= 3.14159;
	x *= 3.14159;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test65()   {
	Interval x(0, 4.3);
	AffineMainVector<T> ax (1,x);
	ax[0] *= 3.14159;
	x *= 3.14159;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

/* ****************************************** */

// Addition
template<class T>
void TestAffineArith<T>::test66() {
	Interval x(1.2,4.5);
	Interval y(3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x+y, ax[0]+ax[1]));
}

template<class T>
void TestAffineArith<T>::test67()  {
	Interval x(1.2,4.5);
	Interval y(-7.2,-3.77);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x+y, ax[0]+ax[1]));
}

template<class T>
void TestAffineArith<T>::test68()   {
	Interval x(1.2,4.5);
	Interval y(-3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x+y, ax[0]+ax[1]));
}

template<class T>
void TestAffineArith<T>::test69()  {
	Interval x(1.2,4.5);
	Interval y(0,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x+y, ax[0]+ax[1]));
}

// Subtraction

template<class T>
void TestAffineArith<T>::test70()   {
	Interval x(1.2,4.5);
	Interval y(3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x-y, ax[0]-ax[1]));
}

template<class T>
void TestAffineArith<T>::test71()   {
	Interval x(1.2,4.5);
	Interval y(-7.2,-3.77);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x-y, ax[0]-ax[1]));
}

template<class T>
void TestAffineArith<T>::test72()  {
	Interval x(1.2,4.5);
	Interval y(-3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x-y, ax[0]-ax[1]));
}

template<class T>
void TestAffineArith<T>::test73()  {
	Interval x(1.2,4.5);
	Interval y(0,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x-y,  ax[0]-ax[1]));
}

// Multiplication


template<class T>
void TestAffineArith<T>::test74()   {
	Interval x(1.2,4.5);
	Interval y(3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x*y,  ax[0]*ax[1]));
}

template<class T>
void TestAffineArith<T>::test75()  {
	Interval x(1.2,4.5);
	Interval y(-7.2,-3.77);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x*y,  ax[0]*ax[1]));
}

template<class T>
void TestAffineArith<T>::test76()  {
	Interval x(1.2,4.5);
	Interval y(-3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x*y,  ax[0]*ax[1]));
}

template<class T>
void TestAffineArith<T>::test77()  {
	Interval x(1.2,4.5);
	Interval y(0,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x*y,  ax[0]*ax[1]));
}

// Division

template<class T>
void TestAffineArith<T>::test78()  {
	Interval x(1.2,4.5);
	Interval y(3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x/y,  ax[0]/ax[1]));
}

template<class T>
void TestAffineArith<T>::test79()  {
	Interval x(1.2,4.5);
	Interval y(-7.2,-3.77);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x/y,  ax[0]/ax[1]));
}

template<class T>
void TestAffineArith<T>::test80()   {
	Interval x(1.2,4.5);
	Interval y(-3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x/y,  ax[0]/ax[1]));
}

template<class T>
void TestAffineArith<T>::test81()  {
	Interval x(1.2,4.5);
	Interval y(0,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	CPPUNIT_ASSERT(compare_results2 ( x/y,  ax[0]/ax[1]));
}

/* ********************************************** */

// Assign-addition with affine forms
template<class T>
void TestAffineArith<T>::test82()  {
	Interval x(1.2,4.5);
	Interval y(3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	ax[0] += ax[1];
	x += y;
	CPPUNIT_ASSERT(compare_results2 (  x,  ax[0]));
}

template<class T>
void TestAffineArith<T>::test83()  {
	Interval x(1.2,4.5);
	Interval y(-7.2,-3.77);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	ax[0] += ax[1];
	x += y;
	CPPUNIT_ASSERT(compare_results2 (  x,  ax[0]));
}

template<class T>
void TestAffineArith<T>::test84()  {
	Interval x(1.2,4.5);
	Interval y(-3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	ax[0] += ax[1];
	x += y;
	CPPUNIT_ASSERT(compare_results2 (  x,  ax[0]));
}

template<class T>
void TestAffineArith<T>::test85()   {
	Interval x(1.2,4.5);
	Interval y(0,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	ax[0] += ax[1];
	x += y;
	CPPUNIT_ASSERT(compare_results2 (  x,  ax[0]));
}

// Assign-multiplication with affine forms
template<class T>
void TestAffineArith<T>::test86() {
	Interval x(1.2,4.5);
	Interval y(3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	ax[0] *= ax[1];
	x *= y;
	CPPUNIT_ASSERT(compare_results2 (  x,  ax[0]));;
}

template<class T>
void TestAffineArith<T>::test87() {
	Interval x(1.2,4.5);
	Interval y(-7.2,-3.77);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	ax[0] *= ax[1];
	x *= y;
	CPPUNIT_ASSERT(compare_results2 (  x,  ax[0]));
}

template<class T>
void TestAffineArith<T>::test88() {
	Interval x(1.2,4.5);
	Interval y(-3.77,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	ax[0] *= ax[1];
	x *= y;
	CPPUNIT_ASSERT(compare_results2 (  x,  ax[0]));
}

template<class T>
void TestAffineArith<T>::test89() {
	Interval x(1.2,4.5);
	Interval y(0,7.2);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	ax[0] *= ax[1];
	x *= y;
	CPPUNIT_ASSERT(compare_results2 (  x,  ax[0]));
}

/* **************************** */
template<class T>
void TestAffineArith<T>::test90() {
	Interval x(92.4,1909.3);
	Interval y(92.4,1909.3);
	IntervalVector v(2);
	v[0]=x;
	v[1]=y;
	AffineMainVector<T> ax(v);
	Interval res = x - y;
	AffineMain<T> resa = ax[0] - ax[1];
	CPPUNIT_ASSERT(compare_results2 ( res,  resa));
}


template<class T>
void TestAffineArith<T>::test91() {
	Interval x(92.4,1909.3);
	AffineMainVector<T> ax (1,x);
	Interval res = x - x;
	AffineMain<T> resa = ax[0] - ax[0];
	CPPUNIT_ASSERT(compare_results ( Interval(0),  resa));
}

template<class T>
void TestAffineArith<T>::test92() {
	Interval x(-92.4,-10.3);
	AffineMainVector<T> ax (1,x);
	Interval res = x - x;
	AffineMain<T> resa = ax[0] - ax[0];
	CPPUNIT_ASSERT(compare_results ( Interval(0),  resa));
}

template<class T>
void TestAffineArith<T>::test93() {
	Interval x(-92.4, 10.3);
	AffineMainVector<T> ax (1,x);
	Interval res = x - x;
	AffineMain<T> resa = ax[0] - ax[0];
	CPPUNIT_ASSERT(compare_results ( Interval(0),  resa));
}

template<class T>
void TestAffineArith<T>::test94() {
	Interval x(0.0, 10.3);
	AffineMainVector<T> ax (1,x);
	Interval res = x - x;
	AffineMain<T> resa = ax[0] - ax[0];
	CPPUNIT_ASSERT(compare_results ( Interval(0),  resa));
}

/* **************************** */
template<class T>
void TestAffineArith<T>::test95() {

	int iter = 1000;
	Interval x(1.3,78.4);
	AffineMainVector<T> ax (1,x);
	for (int i=0; i < iter; i++)
	{
		//cout << "ax : "<< ax << endl;
		ax[0] += -0.5*ax[0];//+0.01;
		x += -0.5*x;//+0.01;
	}
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));

}

template<class T>
void TestAffineArith<T>::test96()  {
	int iter = 1000;
	Interval x(-78.4,-1.3);
	AffineMainVector<T> ax (1,x);
	for (int i=0;i < iter;i++)
	{
		//cout << "ax : "<< ax << endl;
		ax[0] +=-0.5*ax[0];//+0.01;
		x +=-0.5*x;//+0.01;
	}
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test97() {
	int iter = 1000;
	Interval x(-1.3,78.4);
	AffineMainVector<T> ax (1,x);
	for (int i=0;i < iter; i++)
	{
		//cout << "ax : "<< ax << endl;
		ax[0] += -0.5*ax[0];//+0.01;
		x +=-0.5*x;//+0.01;
	}
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

template<class T>
void TestAffineArith<T>::test98()   {
	int iter = 1000;
	Interval x(0.,78.4);
	AffineMainVector<T> ax (1,x);
	for (int i=0;i < iter;i++)
	{
		//cout << "ax : "<< ax << endl;
		ax[0] += -0.5 * ax[0];// + 0.01;
		x += -0.5 * x;// + 0.01;
	}
	//    cout << ax << endl;
	//    cout << ax.itv() << endl;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

/* ********************* */
//test of sqrt
template<class T>
void TestAffineArith<T>::test99()   {
	Interval x(10,11);
	AffineMainVector<T> ax (1,x);
	for (int i=0;i < 1;i++)
	{
		ax[0] = sqrt(ax[0]);
		x = sqrt(x);
	}
	// cout << x << endl;
	// cout << ax.itv() << endl;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}

/* ********************* */
//test of 1/x
template<class T>
void TestAffineArith<T>::test100()   {
	Interval x(3,21);
	AffineMainVector<T> ax (1,x);
	for (int i=0;i < 1;i++)
	{
		ax[0] = 1/ax[0];
		x = 1/x;
	}
	//  cout << x << endl;
	//  cout << ax.itv() << endl;
	CPPUNIT_ASSERT(compare_results ( x, ax[0]));
}


/* ********************* */
template<class T>
void TestAffineArith<T>::test101()   {
	Variable yd(3);
	Interval sigma(10.0);
	Interval rho(28.0);
	Interval beta = Interval(8.0)/3.0;

	Function ydot = Function(yd,Return(sigma*(yd[1]-yd[0]),
			yd[0]*(rho-yd[2])-yd[1],
			yd[0]*yd[1]-beta*yd[2]));

	IntervalVector v(3);
	v[0] = Interval(-1, 1);
	v[1] = Interval(2, 3);
	v[2] = Interval(-4,-4);
	AffineMainVector<T> va(v);

	AffineEval<T> cal(ydot);
	IntervalVector res  = cal.eval(v).v();
	AffineMainVector<T> resa = cal.af2.top->v();

	for (int j = 0; j < 3; ++j) {
		CPPUNIT_ASSERT(compare_results ( res[j], resa[j]));
	}
}



