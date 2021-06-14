/* ============================================================================
 * I B E X - Definition of the Linearizer of Affine2 forms
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENCE.
 *
 * Author(s)   : Jordan Ninin
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */


#ifndef __IBEX_LINEARIZER_AFFINE2_H__
#define __IBEX_LINEARIZER_AFFINE2_H__

#include "ibex_Affine.h"
#include "ibex_System.h"
#include "ibex_Linearizer.h"

//#include <vector>

namespace ibex {

/**
 * \ingroup numeric
 * \brief Affine-based linearization
 *
 * This class is an implementation of the ART algorithm
 * \author Jordan Ninin
 * \date May 2013
 */

class LinearizerAffine2 : public Linearizer {

public:

	LinearizerAffine2 (const System& sys);

	~LinearizerAffine2 ();

	/**
	 * \biref  ART iteration.
	 *
	 *  Linearize the system and performs 2n calls to Simplex in order to reduce the 2 bounds of each variable
	 */
	virtual int linearize(const IntervalVector& box, LPSolver& lp_solver);

	/**
	 * \brief Generation of a linear approximation of the inner region
	 *
	 */
	int inlinearization(const IntervalVector& box, LPSolver& lp_solver);

	/**
	 * \brief Generation of a linear approximation of the linear objective function
	 *
	 */
	bool goal_linearization(const IntervalVector& box, LPSolver& lp_solver);

private:
	/**
	 * \brief The system
	 */
	const System& sys;

	/**
	 * \brief Affine evaluator for the goal function (if any)
	 */
	AffineEval<AF_Default>* goal_af_evl;

	/**
	 * \brief Affine evaluators for the constraints functions
	 */
	AffineEval<AF_Default>** ctr_af_evl;
};

} // end namespace ibex

#endif /* __IBEX_LINEAR_RELAX_AFFINE2_H__ */

