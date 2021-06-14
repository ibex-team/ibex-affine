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

#include "TestCtcAffineFwdBwd.h"

#include "ibex_CtcAffineFwdBwd.h"
#include "Ponts30.h"

using namespace std;

namespace ibex {

void TestCtcAffineFwdBwd::sqrt_issue28() {
	Variable x;
	Function f(x,sqrt(x));
	NumConstraint c(f);

	CtcAffineFwdBwd ctc(c);

	IntervalVector box(1,Interval(-2,-1));

	ctc.contract(box);

	CPPUNIT_ASSERT(box.is_empty());
}




void TestCtcAffineFwdBwd::atan2_issue134() {
	Variable x, y, z;
	Function f(y, x, z, z - atan2(y, x));
	NumConstraint c(f);
	CtcAffineFwdBwd ctc(c);

	IntervalVector box(3);
	box[0] = Interval(1,1);
	box[1] = Interval(0,0);
	box[2] = Interval(-100,100);
	ctc.contract(box);
	check(box[2],Interval::half_pi());
}



void TestCtcAffineFwdBwd::func02() {
	try {
		System sys(SRCDIR_TESTS "/quimper/func02.qpr");

		//cout << "sys nb ctr=" << sys.nb_ctr << endl;
		//CPPUNIT_ASSERT(sys.nb_ctr==12);
		CtcFwdBwd* c[24];
		for (int i=0; i<sys.ctrs.size(); i++)
			c[i]=new CtcFwdBwd(sys.ctrs[i]);

		for (int i=0; i<sys.ctrs.size(); i++) {
			IntervalVector subbox(2);
			subbox[0]=11;
			subbox[1]=12;
			IntervalVector box(8);
			box.put(1,subbox); // load x1[1] and x1[2]
			box.put(4,subbox); // load z1[1] and z1[2]

			c[i]->contract(box);

			CPPUNIT_ASSERT(!box.is_empty());
			check(box.subvector(6,7),subbox); // check x2
		}
		for (int i=0; i<sys.ctrs.size(); i++)
			delete c[i];

	} catch(SyntaxError& e) {
		cout << e << endl;
		CPPUNIT_ASSERT(false);
	}
}


void TestCtcAffineFwdBwd::solve() {
    const ExprSymbol& x=ExprSymbol::new_("x");
	const ExprNode& y=ExprGenericUnaryOp::new_("sinc",x);

    SystemFactory fac;
    fac.add_var(x);
    fac.add_ctr(y=Interval(0.5,0.5));
    System sys(fac);

    CtcAffineFwdBwd c(sys.f_ctrs[0]);
    RoundRobin rr(1e-7);
    CellStack stack;

	Solver solver(sys, c, rr, stack, Vector(1,1e-7), Vector(1,1e-7));
	solver.solve(IntervalVector(1,Interval(-100,100)));

    CPPUNIT_ASSERT(solver.get_data().nb_solution()==2);
    CPPUNIT_ASSERT(solver.get_data().nb_unknown()==1);
    Interval sol=solver.get_data().solution(1)[0];
    CPPUNIT_ASSERT(almost_eq(sin(sol),0.5*sol,1e-6));
}


void TestCtcAffineFwdBwd::ponts30() {
	Ponts30 p30;
	IntervalVector box = p30.init_box;

	NumConstraint* ctr[30];
	CtcAffineFwdBwd* c[30];

	for (int i=0; i<30; i++) {
		Function* fi=dynamic_cast<Function*>(&((*p30.f)[i]));
		CPPUNIT_ASSERT(fi!=NULL);
		ctr[i]=new NumConstraint(*fi,EQ);
		c[i]=new CtcAffineFwdBwd(*ctr[i]);
	}
//	cout << "before="<< box << endl;

	for (int i=0; i<30; i++) {
		c[i]->contract(box);
		//cout << box << endl;
		//cout << p30.hc4r_box[i] << endl;
		//c[i]->hc4r.eval.f.cf.print<Domain>();
		CPPUNIT_ASSERT(almost_eq(box, p30.hc4r_box[i],1e-02));
	}
	//cout << "after="<< box << endl;
//	box = p30.init_box;
//
//	Array<NumConstraint> a(ctr,30);
//	CtcHC4 hc4(a,0.1);
//	hc4.accumulate=true;
//	box=p30.init_box;
//	hc4.contract(box);
//
//	CPPUNIT_ASSERT(almost_eq(box, p30.hc4_box,1e-04));

	for (int i=0; i<30; i++) {
		delete c[i];
		delete ctr[i];
	}
}

} // namespace ibex
