/* ============================================================================
 * I B E X - Affine Vector Tests
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENSE
 *
 * Author(s)   : Jordan NININ
 * Created     : Juin 20, 2021
 * ---------------------------------------------------------------------------- */

#include "TestAffineVector.h"
#include "ibex_Affine.h"
#include "utils.h"

using namespace std;


template <class T>
void TestAffineVector<T>::cons01() {
	IntervalVector x(2);
	x[0]=Interval::ALL_REALS;
	x[1]=Interval::ALL_REALS;
	check_af(AffineVarMainVector<T>(2),x);
	check_af(AffineVarMainVector<T>(x),x);
	check_af(AffineMainVector<T>(2),x);
	check_af(AffineMainVector<T>(2)=x,x);
}

template <class T>
void TestAffineVector<T>::cons02() {
	IntervalVector x(2);
	x[0]=(Interval(0,1));
	x[1]=(Interval(0,1));
	check_af(AffineVarMainVector<T>(2,Interval(0,1)),x);
	check_af(AffineVarMainVector<T>(x),x);
	AffineMainVector<T> tmp(2);
	tmp.init(Interval(0,1));
	check_af(tmp,x);
	AffineVarMainVector<T> ttt(x);
	check_af(AffineMainVector<T>(ttt),x);
	check_af(AffineMainVector<T>(2)=x,x);
	AffineVarMainVector<T> tmp2(2);
	tmp2.init(Interval(0,1));
	check_af(tmp2,x);
}

template <class T>
void TestAffineVector<T>::cons03() {
	IntervalVector x(2);
	x[0]=Interval(0,1);
	x[1]=Interval(2,3);
	AffineVarMainVector<T> af(x);
	check_af(af,x);
	check_af(af,AffineMainVector<T>(2)=x);
}

template <class T>
void TestAffineVector<T>::cons04() {
	double bounds[][2] = {{0,1},{2,3}};
	IntervalVector x(2);
	x[0]=Interval(0,1);
	x[1]=Interval(2,3);
	AffineVarMainVector<T> af(2,bounds);
	check_af(af,x);
	check_af(AffineMainVector<T>(2)=af,x);
}

template <class T>
void TestAffineVector<T>::cons05() {
	AffineVarMainVector<T>  x(2);
	x[0].set_empty();
	x[1].set_empty();
	check_af(AffineMainVector<T>::empty(2),IntervalVector::empty(2));
	check_af(x,IntervalVector::empty(2));
	CPPUNIT_ASSERT(x.is_empty());
	check_af(x,IntervalVector(2)=x.itv());
	CPPUNIT_ASSERT((AffineMainVector<T>(2)=x).is_empty());
}

template <class T>
void TestAffineVector<T>::set_empty01() {
	AffineVarMainVector<T> x(2);
	CPPUNIT_ASSERT(!x.is_empty());
	x.set_empty();
	CPPUNIT_ASSERT(x.is_empty());
	AffineMainVector<T> x2(2);
	CPPUNIT_ASSERT(!x2.is_empty());
	x2.set_empty();
	CPPUNIT_ASSERT(x2.is_empty());
}



template <class T>
void TestAffineVector<T>::is_empty02() {
	CPPUNIT_ASSERT(!AffineMainVector<T>(2).is_empty());
	CPPUNIT_ASSERT(!AffineVarMainVector<T>(2).is_empty());
}

template <class T>
void TestAffineVector<T>::resize01() {
	AffineVarMainVector<T> x(1,Interval(1,2));
	x.resize(3);
	CPPUNIT_ASSERT(x.size()==3);
	check_af(x[0],Interval(1,2));
	check_af(x[1],Interval::ALL_REALS);
	check_af(x[2],Interval::ALL_REALS);
	AffineMainVector<T> x2(1);
	x2[0] = Interval(1,2);
	x2.resize(3);
	CPPUNIT_ASSERT(x2.size()==3);
	check_af(x2[0],Interval(1,2));
	check_af(x2[1],Interval::ALL_REALS);
	check_af(x2[2],Interval::ALL_REALS);
}

template <class T>
void TestAffineVector<T>::resize02() {
	AffineVarMainVector<T> x(1,Interval(1,2));
	x.resize(1);
	CPPUNIT_ASSERT(x.size()==1);
	check_af(x[0],Interval(1,2));
	AffineMainVector<T> x2(1);
	x2[0] = Interval(1,2);
	x2.resize(1);
	CPPUNIT_ASSERT(x2.size()==1);
	check_af(x2[0],Interval(1,2));
}

template <class T>
void TestAffineVector<T>::resize03() {
	AffineVarMainVector<T> x(2,Interval(1,2));
	x.set_empty();
	x.resize(3);
	CPPUNIT_ASSERT(x.size()==3);
	CPPUNIT_ASSERT(x.is_empty());
	CPPUNIT_ASSERT(x[2]==Interval::ALL_REALS);
	AffineMainVector<T> x2(2);
	x2[0]= Interval(1,2);
	x2.set_empty();
	x2.resize(3);
	CPPUNIT_ASSERT(x2.size()==3);
	CPPUNIT_ASSERT(x2.is_empty());
	CPPUNIT_ASSERT(x2[2]==Interval::ALL_REALS);
}

template <class T>
void TestAffineVector<T>::resize04() {
	IntervalVector x(5);
	x[0]=Interval(1,2);
	x[1]=Interval(3,4);
	AffineVarMainVector<T> af(x);
	AffineMainVector<T> af2(af);

	af.resize(2);
	CPPUNIT_ASSERT(af.size()==2);
	check_af(af[0],Interval(1,2));
	check_af(af[1],Interval(3,4));
	af.resize(5);
	CPPUNIT_ASSERT(af.size()==5);
	check_af(af[0],Interval(1,2));
	check_af(af[1],Interval(3,4));
	af.resize(1);
	CPPUNIT_ASSERT(af.size()==1);
	check_af(af[0],Interval(1,2));

	CPPUNIT_ASSERT(af2.size()==5);
	af2.resize(2);
	CPPUNIT_ASSERT(af2.size()==2);
	check_af(af2[0],Interval(1,2));
	check_af(af2[1],Interval(3,4));
	af2.resize(5);
	CPPUNIT_ASSERT(af2.size()==5);
	check_af(af2[0],Interval(1,2));
	check_af(af2[1],Interval(3,4));
	af2.resize(1);
	CPPUNIT_ASSERT(af2.size()==1);
	check_af(af[0],Interval(1,2));
}

static double _x[][2]={{0,1},{2,3},{4,5}};

template <class T>
void TestAffineVector<T>::subvector01() {
	double _x01[][2]={{0,1},{2,3}};
	CPPUNIT_ASSERT(AffineVarMainVector<T>(3,_x).subvector(0,1)==IntervalVector(2,_x01));
}

template <class T>
void TestAffineVector<T>::subvector02() {
	double _x12[][2]={{2,3},{4,5}};
	CPPUNIT_ASSERT(AffineVarMainVector<T>(3,_x).subvector(1,2)==IntervalVector(2,_x12));
}

template <class T>
void TestAffineVector<T>::subvector03() {
	double _x11[][2]={{2,3}};
	CPPUNIT_ASSERT(AffineVarMainVector<T>(3,_x).subvector(1,1)==IntervalVector(1,_x11));
}

template <class T>
void TestAffineVector<T>::subvector04() {
	double _x22[][2]={{4,5}};
	CPPUNIT_ASSERT(AffineVarMainVector<T>(3,_x).subvector(2,2)==IntervalVector(1,_x22));
}

template <class T>
void TestAffineVector<T>::subvector05() {
	CPPUNIT_ASSERT(AffineVarMainVector<T>(3,_x).subvector(0,2)==IntervalVector(3,_x));
}


template <class T>
void TestAffineVector<T>::eq01() {
	AffineVarMainVector<T> x(3,_x);
	CPPUNIT_ASSERT(x==x);
	CPPUNIT_ASSERT(!(x!=x));
}

template <class T>
void TestAffineVector<T>::eq02() {
	AffineVarMainVector<T> x(3,_x);
	double _x01[][2]={{0,1},{2,3}};
	AffineVarMainVector<T> x1(2,_x01);
	CPPUNIT_ASSERT(!(x==x1));
	CPPUNIT_ASSERT(x!=x1);
}

template <class T>
void TestAffineVector<T>::eq03() {
	double _x1[][2]={{0,1},{4,5}};
	double _x2[][2]={{2,3},{6,7}};
	AffineVarMainVector<T> x1(2,_x1);
	AffineVarMainVector<T> x2(2,_x2);
	x1.set_empty();
	x2.set_empty();
	CPPUNIT_ASSERT(x1==x2);
	CPPUNIT_ASSERT(!(x1!=x2));
}

template <class T>
void TestAffineVector<T>::mid01() {
	AffineVarMainVector<T> x(3,_x);
	Vector m=x.mid();
	check(m[0],0.5);
	check(m[1],2.5);
	check(m[2],4.5);
}

template <class T>
void TestAffineVector<T>::is_unbounded02() {
	double _x1[][2]={{0,1},{0,2},{NEG_INFINITY,0}};
	CPPUNIT_ASSERT(AffineVarMainVector<T>(3,_x1).is_unbounded());
}

template <class T>
void TestAffineVector<T>::is_unbounded03() {
	double _x1[][2]={{0,1},{0,2}};
	CPPUNIT_ASSERT(!AffineVarMainVector<T>(2,_x1).is_unbounded());
}

template <class T>
void TestAffineVector<T>::is_unbounded04() {
	CPPUNIT_ASSERT(AffineVarMainVector<T>(1).is_unbounded());
}

template <class T>
void TestAffineVector<T>::minus01() {
	double _x1[][2]={{0,3},{0,2},{0,1}};
	double _x2[][2]={{-3,0},{-2,0},{-1,0}};
	CPPUNIT_ASSERT((-AffineVarMainVector<T>(3,_x1))==AffineVarMainVector<T>(3,_x2));
}

template <class T>
void TestAffineVector<T>::minus02() {
	double _x1[][2]={{0,1},{0,POS_INFINITY}};
	double _x2[][2]={{-1,0},{NEG_INFINITY,0}};
	check_af(-AffineVarMainVector<T>(2,_x1),IntervalVector(2,_x2));
}


template <class T>
void TestAffineVector<T>::add01() {
	double _x1[][2]={{0,3},{0,2},{0,1}};
	double _x2[][2]={{0,1},{0,1},{0,1}};
	double _x3[][2]={{0,4},{0,3},{0,2}};

	AffineVarMainVector<T> x1(3,_x1);
	AffineVarMainVector<T> x2(3,_x2);
	AffineVarMainVector<T> x3(3,_x3);
	AffineMainVector<T> e(3);
	e.set_empty();

	check_af(x1+x2,x3);
	check_af(x1+x2+x1-x1,x3);
	check_af(x1+x2,x3+x2-x2);
	//check_af(x1+e,e);
	CPPUNIT_ASSERT((x1+e).is_empty());
	//check_af(IntervalVector(x1)+=e,e);
	CPPUNIT_ASSERT((AffineMainVector<T>(x1)+=e).is_empty());

	//check_af(e+x1,e);
	CPPUNIT_ASSERT((e+x1).is_empty());
	//check_af(e+=x1,e);
	CPPUNIT_ASSERT((e+=x1).is_empty());
	//check_af(e+e,e);
	CPPUNIT_ASSERT((e+e).is_empty());
	//check_af(e+=e,e);
	CPPUNIT_ASSERT((e+=e).is_empty());

	check_af(AffineMainVector<T>(x1)+=x2,x3);
	//check_af(IntervalVector(x1)+=e,e);
	CPPUNIT_ASSERT((AffineMainVector<T>(x1)+=e).is_empty());

	check_af(AffineMainVector<T>(x2)+=x1,x3);
}

template <class T>
void TestAffineVector<T>::sub01() {
	double _x1[][2]={{0,3},{0,2},{0,1}};
	double _x2[][2]={{0,1},{0,1},{0,1}};
	double _x3[][2]={{0,2},{0,1},{0,0}};
	AffineVarMainVector<T> x1(3,_x1);
	AffineVarMainVector<T> x2(3,_x2);
	AffineVarMainVector<T> x3(3,_x3);
	AffineMainVector<T> e(3);
	e.set_empty();

	check_af(x1-x2,x3);
	check_af(x1-x2+x2-x2,x3);
	check_af(x1-x2,x3+x1-x1);
	check_af(x2-x1,-x3);
	//check_af(x1-e,e);
	CPPUNIT_ASSERT((x1-e).is_empty());
	//check_af(IntervalVector(x1)-=e,e);
	CPPUNIT_ASSERT((AffineMainVector<T>(x1)-=e).is_empty());

	//check_af(e-x1,e);
	CPPUNIT_ASSERT((e-x1).is_empty());
	//check_af(e-=x1,e);
	CPPUNIT_ASSERT((e-=x1).is_empty());
	//check_af(e-e,e);
	CPPUNIT_ASSERT((e-e).is_empty());
	//check_af(e-=e,e);
	CPPUNIT_ASSERT((e-=e).is_empty());

	check_af(AffineMainVector<T>(x1)-=x2,x3);
	check_af(AffineMainVector<T>(x2)-=x1,-x3);
}
