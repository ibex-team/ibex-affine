/* ============================================================================
 * I B E X - HC4 Revise with Affine arithmetic (forward-backward algorithm)
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENCE.
 *
 * Author(s)   : Jordan Ninin
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */

#include "ibex_CtcAffineFwdBwd.h"

namespace ibex {


CtcAffineFwdBwd::CtcAffineFwdBwd(const Function& f, CmpOp op): CtcFwdBwd(f,op),
		 _myaffineeval(f), _myhc4revise(_myaffineeval) { }

CtcAffineFwdBwd::CtcAffineFwdBwd(const Function& f, const Domain& y): CtcFwdBwd(f,y),
		 _myaffineeval(f), _myhc4revise(_myaffineeval) { }

CtcAffineFwdBwd::CtcAffineFwdBwd(const Function& f, const Interval& y): CtcFwdBwd(f,y),
		 _myaffineeval(f), _myhc4revise(_myaffineeval) { }

CtcAffineFwdBwd::CtcAffineFwdBwd(const Function& f, const IntervalVector& y): CtcFwdBwd(f,y),
		 _myaffineeval(f), _myhc4revise(_myaffineeval) { }

CtcAffineFwdBwd::CtcAffineFwdBwd(const Function& f, const IntervalMatrix& y): CtcFwdBwd(f,y),
		 _myaffineeval(f), _myhc4revise(_myaffineeval) { }

CtcAffineFwdBwd::CtcAffineFwdBwd(const NumConstraint& ctr): CtcFwdBwd(ctr),
		 _myaffineeval(ctr.f), _myhc4revise(_myaffineeval) { }

CtcAffineFwdBwd::CtcAffineFwdBwd(const System& sys, int i): CtcFwdBwd(sys,i),
		 _myaffineeval(sys.ctrs[i].f), _myhc4revise(_myaffineeval) { }


void CtcAffineFwdBwd::contract(IntervalVector& box) {
	ContractContext context(box);
	contract(box,context);
}

CtcAffineFwdBwd::~CtcAffineFwdBwd(){ }

void CtcAffineFwdBwd::contract(IntervalVector& box, ContractContext& context) {


	assert(box.size()==_myaffineeval._af_f.nb_var());

	BxpActiveCtr* p=NULL;
	BxpSystemCache* sp=NULL;

//
// Notes:
// Using properties here is counter-productive on many benchmarks just because of the cost
// of looking inside a map (operator[]). See issue #342
//
//	if (ctr_num==-1)
//		p=(BxpActiveCtr*) context.prop[active_prop_id];
//	else
//		sp=(BxpSystemCache*) context.prop[system_cache_id];
//

	if ((p && !p->active()) || (sp && sp->is_active_ctrs_uptodate() && !(sp->active_ctrs()[ctr_num]))) {
		context.output_flags.add(INACTIVE);
		context.output_flags.add(FIXPOINT);
		return;
	}

	//std::cout << " hc4 of " << f << "=" << d << " with box=" << box << std::endl;
	if (_myhc4revise.proj(d,box)) {
		if (p) p->set_inactive();
		if (sp) sp->active_ctrs().remove(ctr_num);
		context.output_flags.add(INACTIVE);
		context.output_flags.add(FIXPOINT);
	}
	else if (box.is_empty()) {
		context.output_flags.add(FIXPOINT);
	}
	//std::cout << " ---> " << box << std::endl;

	// Note: setting the FIXPOINT flag is incorrect when there
	// is no multiple occurrence because some operators
	// may be non-optimal in backward mode.

	//TODO: set used_vars in output impact
	context.prop.update(BoxEvent(box,BoxEvent::CONTRACT));
}

} // namespace ibex
