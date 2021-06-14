/* ============================================================================
 * I B E X - Definition of the Affine3 class based on fAFFull version 1 from DynIbex
 * ============================================================================
 * Copyright   : ENSTA Bretagne
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENSE.
 *
 * Author(s)   : Jordan Ninin, Julien Alexandre dit Sandretto and Alexandre Chapoutot
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */

#ifndef IBEX_AFFINE3_FAFFULLI_H_
#define IBEX_AFFINE3_FAFFULLI_H_

#include "ibex_Interval.h"

#include <list>
#include <ostream>
#include <utility>


namespace ibex {


typedef unsigned long int af3_int;
class AF_fAFFullI;


//typedef AF_fAF1  AF_Other;
//typedef AF_fAF2  AF_Other;
//typedef AF_fAF2_fma  AF_Other;
//typedef AF_iAF  AF_Other;
//typedef AF_sAF  AF_Other;
//typedef AF_No  AF_Other;
typedef AF_fAFFullI AF_Other;



template<class T>  class AffineMain;
template<class T>  class AffineVarMain;

class AF_fAFFullI {


private:
	friend class AffineVarMain<AF_fAFFullI>;
	friend class AffineMain<AF_fAFFullI>;
	template<class A>
	friend std::ostream& operator<<(std::ostream& os, const AffineMain<A>& x);

	static af3_int _counter;
	static double maTol;
	/**
	 * Code for the particular case:
	 * if the affine form is actif, _actif=1  and _n is the size of the affine form
	 * if the set is degenerate, _actif = 0 and itv().diam()< AF_EC
	 * if the set is empty, _actif = -1
	 * if the set is ]-oo,+oo[, _actif = -2 and _err =]-oo,+oo[
	 * if the set is [a, +oo[ , _actif = -3 and _err = [a, +oo[
	 * if the set is ]-oo, a] , _actif = -4 and _err = ]-oo, a]
	 *
	 */
	double _center;
	std::list< std::pair<af3_int, double> > _rays;
	Interval _garbage;


	AF_fAFFullI (double center, std::list<std::pair<af3_int,double> > rays, Interval garbage);

public:
	/** \brief  Delete the affine form */
	virtual ~AF_fAFFullI();


};

inline AF_fAFFullI::AF_fAFFullI (double center, std::list<std::pair<af3_int,double> > rays, Interval garbage) :
    		_center(center), _rays(rays), _garbage(garbage) {

}

inline AF_fAFFullI::~AF_fAFFullI() {
	if (!_rays.empty()) { _rays.clear(); }
}




}



#endif /* IBEX_AFFINE3_FAFFULL_H_ */
