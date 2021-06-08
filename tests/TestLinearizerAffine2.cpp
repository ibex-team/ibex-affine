//============================================================================
//                                  I B E X                                   
// File        : TestLinearizerAffine2.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Apr 10, 2012
// Last Update : Apr 10, 2012
//============================================================================

#include "TestLinearizerAffine2.h"

#include "ibex_CtcFwdBwd.h"
#include "ibex_SystemFactory.h"
#include "ibex_System.h"
#include "ibex_CtcPolytopeHull.h"
#include "ibex_LinearizerXTaylor.h"
#include "ibex_Array.h"

using namespace std;

namespace ibex {



void TestLinearizerAffine2::fixbug01() {

	SystemFactory f;
	Variable x,extra;
	f.add_var(x); f.add_var(extra); f.add_ctr(x<=0);
	f.add_ctr(extra>=0); f.add_ctr(extra<= 15);
	System sys(f);
	// if I replace the domain of "extra" by [1,1.1]
	// there is no more problem.
	double _box[][2] = {{-100,100},{1,1}};
	IntervalVector box(2,_box);
	LinearizerAffine2 linear_relax(sys);

	CtcPolytopeHull polytope(linear_relax);
	polytope.contract(box);

	double _box2[][2] = {{-100,0},{1,1}};
	IntervalVector box2(2,_box2);
	check(box,box2);

	double _e[][2] = {{-100,100},{10,20}};
	IntervalVector e(2,_e);

	polytope.contract(e);

	double _ee[][2] = {{-100,0},{10,15}};
	IntervalVector ee2(2,_ee);
	check(e,ee2);
}


} // end namespace ibex
