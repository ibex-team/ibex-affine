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

#ifndef __IBEX_CTC_AFFINEFWDBWD_H__
#define __IBEX_CTC_AFFINEFWDBWD_H__

#include "ibex_Ctc.h"
#include "ibex_HC4Revise.h"
#include "ibex_Affine.h"
#include "ibex_CtcFwdBwd.h"
#include "ibex_BxpActiveCtr.h"
#include "ibex_BxpSystemCache.h"

namespace ibex {

/**
 * \ingroup contractor
 * \brief Forward-backward contractor (HC4Revise).
 *
 */
class CtcAffineFwdBwd: public CtcFwdBwd {
public:
	/**
	 * \brief Build the contractor for "f(x)=0" or "f(x)<=0".
	 *
	 * \param op: by default: EQ.
	 *
	 */
	CtcAffineFwdBwd(const Function& f, CmpOp op=EQ);

	/**
	 * \brief Build the contractor for "f(x) in [y]".
	 */
	CtcAffineFwdBwd(const Function& f, const Domain& y);

	/**
	 * \brief Build the contractor for "f(x) in [y]".
	 */
	CtcAffineFwdBwd(const Function& f, const Interval& y);

	/**
	 * \brief Build the contractor for "f(x) in [y]".
	 */
	CtcAffineFwdBwd(const Function& f, const IntervalVector& y);

	/**
	 * \brief Build the contractor for "f(x) in [y]".
	 */
	CtcAffineFwdBwd(const Function& f, const IntervalMatrix& y);

	/**
	 * \remark ctr is not kept by reference.
	 */
	CtcAffineFwdBwd(const NumConstraint& ctr);

	/**
	 * \brief Build the contrator for the ith constraint
	 *
	 * Allow to benefit from the system cache associated
	 * to the system.
	 */
	CtcAffineFwdBwd(const System& sys, int i);

	/**
	 * \brief Delete this.
	 */
	~CtcAffineFwdBwd();

	/**
	 * \brief Contract a box.
	 */
	void contract(IntervalVector& box);

	/**
	 * \brief Contract the box.
	 */
	virtual void contract(IntervalVector& box, ContractContext& context);

protected:
	Affine2Eval _myaffineeval;
	HC4Revise 	_myhc4revise;
};

} // namespace ibex
#endif // __IBEX_CTC_FWDBWD_H__
