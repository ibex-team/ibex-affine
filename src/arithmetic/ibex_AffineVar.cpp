/* ============================================================================
 * I B E X - Implementation of the AffineVarMain class
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENCE.
 *
 * Author(s)   : Jordan Ninin
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */
#include "ibex_AffineVar.h"
#include "ibex_AffineMain.h"
#include <iostream>
#include <cassert>

namespace ibex {


/**
 * Code for the particular case:
 * if the affine form is actif, _actif=1  and _n is the size of the affine form
 * if the set is degenerate, _actif = 0 and itv().diam()< AF_EC
 * if the set is empty, _actif = -1
 * if the set is ]-oo,+oo[, _actif = -2 and _err = ]-oo,+oo[
 * if the set is [a, +oo[ , _actif = -3 and _err = [a, +oo[
 * if the set is ]-oo, a] , _actif = -4 and _err = ]-oo, a]
 *
 */

//===========================================================================================
//===========================================================================================




template<>
AffineVarMain<AF_fAF2>& AffineVarMain<AF_fAF2>::operator=(const Interval& x) {

	if (x.is_empty()) {
		_actif = -1;
		_elt._err = 0.0;
	} else if (x.ub()>= POS_INFINITY && x.lb()<= NEG_INFINITY ) {
		_actif = -2;
		_elt._err = 0.0;
	} else if (x.ub()>= POS_INFINITY ) {
		_actif = -3;
		_elt._err = x.lb();
	} else if (x.lb()<= NEG_INFINITY ) {
		_actif = -4;
		_elt._err = x.ub();
	} else  {
		assert((unsigned long int)_n > _var);
		_elt._err = 0.0;
		if (_elt._val==NULL) _elt._val = new double[_n+1];
		_elt._val[0] = x.mid();
		for (int i = 1; i <= size(); i++){
			_elt._val[i] = 0.0;
		}
		if (x.is_degenerated()){
			_actif=0;
		} else {
			_actif=1;
			_elt._val[_var+1] = x.rad();
		}
	}
	return *this;

}


template<>
AffineVarMain<AF_fAF2>& AffineVarMain<AF_fAF2>::operator=(const AffineVarMain<AF_fAF2>& x) {
	if (this != &x) {
		_var = x._var;
		_elt._err = x._elt._err;
		_actif = x._actif;
		_n =x._n;
		if (x.is_actif()) {
			if (_elt._val!=NULL) { delete[] _elt._val; }
			_elt._val = new double[_n+1];
			int i = 0;
			for (; i <= x.size(); i++) {
				_elt._val[i] = x._elt._val[i];
			}
		}
	}
	return *this;
}



//template<>
//AffineVarMain<AF_fAF2>::AffineVarMain(const Interval & itv) :
//		AffineMain<AF_fAF2>(_count+1, _count, itv),
//		var		(_count) {
//	_count++;
//}
//
//template<>
//AffineVarMain<AF_fAF2>::AffineVarMain(double d) :
//		AffineMain<AF_fAF2>(_count+1, _count, Interval(d)),
//		var		(_count) {
//	_count++;
//}







//===========================================================================================
//===========================================================================================

template<>
AffineVarMain<AF_fAFFullI>& AffineVarMain<AF_fAFFullI>::operator=(const Interval& x) {
	_elt._garbage = Interval(0.0); //
	if (x.is_empty()) {
		_actif = -1;
		_elt._center = 0.0;
		// Note the empty Affine form is an empty list
		if (!_elt._rays.empty()) {
			_elt._rays.clear();
		}
	}
	else if (x.ub()>= POS_INFINITY && x.lb()<= NEG_INFINITY ) {
		_actif = -2;
		_elt._center = 0.0;
		// Note the entire set in Affine form is an empty list
		if (!_elt._rays.empty()) {
			_elt._rays.clear();
		}
	}
	else if (x.ub()>= POS_INFINITY ) {
		_actif = -3;
		_elt._center = x.lb();
		if (!_elt._rays.empty()) {
			_elt._rays.clear();
		}
	}
	else if (x.lb()<= NEG_INFINITY ) {
		_actif = -4;
		_elt._center = x.ub();
		if (!_elt._rays.empty()) {
			_elt._rays.clear();
		}
	}
	else  {
		if (!_elt._rays.empty()) {
			_elt._rays.clear();
		}
		_elt._center = x.mid();
		if ( x.is_degenerated()) {
			_actif = 0;
		} else {
			_actif = 1;
			std::pair<af3_int,double> p(_var, x.rad());
			_elt._rays.push_back(p);
		}
		if (AF_fAFFullI::_counter <= _var ) {AF_fAFFullI::_counter = _var+1;}

	}
	return *this;
}

template<>
AffineVarMain<AF_fAFFullI>& AffineVarMain<AF_fAFFullI>::operator=(const AffineVarMain<AF_fAFFullI>& x) {
	if (this != &x) {
		_var = x._var;
		_n = x._n;
		_actif = x._actif;
		_elt._center = x._elt._center;
		_elt._garbage = x._elt._garbage;
		_elt._rays.clear();

		if ((x.is_actif()) && (!x._elt._rays.empty()))	{
			std::list<std::pair<af3_int,double> >::const_iterator it = x._elt._rays.begin();
			for (; it != x._elt._rays.end(); ++it) {
				_elt._rays.push_back(std::pair<af3_int,double>(it->first,it->second));
			}
		}
	}
	return *this;
}


//template<>
//AffineVarMain<AF_fAFFullI>::AffineVarMain(const Interval & itv) :
//		AffineMain<AF_fAFFullI>(),
//		var		(_count) {
//	*this = itv;
//	_count++;
//}
//
//template<>
//AffineVarMain<AF_fAFFullI>::AffineVarMain(double d) :
//		AffineMain<AF_fAFFullI>(),
//		var		(_count) {
//	*this = Interval(d);
//	_count++;
//}
//
//
//template<>
//AffineVarMain<AF_fAFFullI>::AffineVarMain(int size, int var1, const Interval& itv) :
//		AffineMain<AF_fAFFullI>(),
//		var		(var1) {
//	*this = itv;
//}
//
//




//===========================================================================================
//===========================================================================================



}// end namespace ibex




