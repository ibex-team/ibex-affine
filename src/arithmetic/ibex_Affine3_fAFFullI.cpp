/* ============================================================================
 * I B E X - Definition of the AffineMain<AF_fAFFull> class based on fAFFull version 1 of DynIbex
 * ============================================================================
 * Copyright   : ENSTA ParisTech
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Julien Alexandre dit Sandretto, Alexandre Chapoutot and Jordan Ninin
 * Created     : Jul 18, 2014
 * Sponsored   : This research benefited from the support of the "Chair Complex Systems Engineering - Ecole Polytechnique, THALES, DGA, FX, DASSAULT AVIATION, DCNS Research, ENSTA ParisTech, Telecom ParisTech, Fondation ParisTech and FDO ENSTA"
 * ---------------------------------------------------------------------------- */
#include "ibex_Affine3_fAFFullI.h"
#include "ibex_Affine.h"
//#include <iomanip>
#include <cassert>

namespace ibex {

double AF_fAFFullI::maTol = 2.5e-15;

unsigned long int AF_fAFFullI::_counter = 0;

bool noise_null (const std::pair<int,double> value) { return (value.second >= 0)&&(value.second <= 0.0); }


//TODO  Changer la list en vector pour faire bien le AffineMain::val(i)
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
template<>
AffineMain<AF_fAFFullI>::AffineMain() :
	_actif (-2     ),
	_n		(0		),
	_elt	(0.0, std::list<std::pair<int,double> >(), Interval(0.0))	{
}

template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator=(const Interval& x) {
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
			std::pair<int,double> p(AF_fAFFullI::_counter++, x.rad());
			_elt._rays.push_back(p);
		}
	}
	return *this;
}


template<>
AffineMain<AF_fAFFullI>::AffineMain(int size, int var, const Interval& itv) :
	_actif	(0),
	_n 		(size),
	_elt	(0, std::list<std::pair<int,double> >(), Interval(0.0))  {
	assert((size>=0) && (var>=0) && (var<=size));
	if (!(itv.is_unbounded()||itv.is_empty())) {
		_elt._center = itv.mid();
		if (! itv.is_degenerated()) {
			_actif =1;
			std::pair<int,double> p(var, itv.rad());
			_elt._rays.push_back(p);
			if (AF_fAFFullI::_counter <= (unsigned long int)var ) {AF_fAFFullI::_counter = var+1;}

		}
	} else {
		*this = itv;
	}
}



template<>
AffineMain<AF_fAFFullI>::AffineMain(const AffineMain<AF_fAFFullI>& x) :
		_actif	(x._actif),
		_n		(x._n),
		_elt	(x._elt._center, std::list<std::pair<int,double> >(),x._elt._garbage) {
	if (!x._elt._rays.empty())	{
		std::list<std::pair<int,double> >::const_iterator it = x._elt._rays.begin();
		for (; it != x._elt._rays.end(); ++it) {
			_elt._rays.push_back(std::pair<int,double>(it->first,it->second));
		}
	}

}


template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator=(const AffineMain<AF_fAFFullI>& x) {
	if (this != &x) {
		_n = x._n;
		_actif = x._actif;
		_elt._center = x._elt._center;
		_elt._garbage = x._elt._garbage;
		_elt._rays.clear();

		if ((x.is_actif()) && (!x._elt._rays.empty()))	{
			std::list<std::pair<int,double> >::const_iterator it = x._elt._rays.begin();
			for (; it != x._elt._rays.end(); ++it) {
				_elt._rays.push_back(std::pair<int,double>(it->first,it->second));
			}
		}
	}
	return *this;
}

template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator=(double d) {
	_elt._garbage = Interval(0.0);
	if (!_elt._rays.empty()) {
		_elt._rays.clear();
	}

	if (fabs(d)<POS_INFINITY) {
		_actif = 0;
		_elt._center = d;

	} else {
		if (d>0) {
			_actif = -3;
			_elt._center = d;
		} else {
			_actif = -4;
			_elt._center = d;
		}
	}
	return *this;
}



template<>
double AffineMain<AF_fAFFullI>::val(int i) const{  // TODO changer la std::list par une autre structure pour ameliorer
	assert((0<=i) && (i<=size()));
	if (!_elt._rays.empty()) {
		std::list<std::pair<int,double> >::const_iterator iter = _elt._rays.begin();
		for (; iter != _elt._rays.end(); ++iter) {
			if (iter -> first == i) { return iter -> second; }
			if (iter -> first > i) { return 0.0; }
		}
	}
	return 0.0;
}

template<>
double AffineMain<AF_fAFFullI>::err() const{
	return _elt._garbage.rad();
}



template<>
int AffineMain<AF_fAFFullI>::size() const {
	if (!_elt._rays.empty())	{
		std::pair<int,double> p = _elt._rays.back();
		if (_n>p.first) {
			return _n;
		} else {
			return p.first;
		}
	}
	else {
		return _n;
	}
}

template<>
const Interval AffineMain<AF_fAFFullI>::itv() const {

	switch(_actif) {
	case -1 : {
		return Interval::empty_set();
		break;
	}
	case -2 : {
		return Interval::all_reals();
		break;
	}
	case -3 : {
		return Interval(_elt._center,POS_INFINITY);
		break;
	}
	case -4: {
		return Interval(NEG_INFINITY,_elt._center);
		break;
	}
	case 0: {
		return Interval(_elt._center);
		break;
	}
	default: { // _actif==1
		Interval res(_elt._center);
		Interval pmOne(-1.0, 1.0);
		if (!_elt._rays.empty()) {
			std::list<std::pair<int,double> >::const_iterator it = _elt._rays.begin();
			for (; it != _elt._rays.end(); ++it) {
				res += (it -> second * pmOne);
			}
		}
		res += _elt._garbage;
		return res;
		break;
	}
	}
}

template<>
double AffineMain<AF_fAFFullI>::mid() const{
	return (is_actif()) ? _elt._center : itv().mid();
}






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

template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::Aneg() {
	switch(_actif) {
	case -3 : {
		_elt._center = -_elt._center;
		_actif    = -4;
		break;
	}
	case -4 : {
		_elt._center = -_elt._center;
		_actif    = -3;
		break;
	}
	case 0 :{
		_elt._center = (-_elt._center);
		break;
	}
	case 1 : {
		_elt._center = -_elt._center;
		_elt._garbage = -_elt._garbage;
		if (!_elt._rays.empty()) {
			std::list<std::pair<int,double> >::iterator it = _elt._rays.begin();
			for (; it != _elt._rays.end(); ++it) {
				it->second = -(it->second);
			}
		}
		break;
	}
	default :
		break;
	}

	return *this;
}




template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator*=(double alpha) {
	if (_actif==1) {  // multiply by a scalar alpha
		if (alpha==0.0) {
			_actif = 0;
			_elt._center = 0.0;
			_elt._rays.clear();
			_elt._garbage = Interval(0.0);
		} else if ((fabs(alpha)) < POS_INFINITY) {
			Interval roundoff_error(0.,0.);
			Interval intermediate(0.,0.);
			// Computation step for the center
			intermediate = Interval(_elt._center) * alpha;
			_elt._center = intermediate.mid();
			roundoff_error += intermediate.rad();

			_elt._garbage *= alpha;

			// Computation step for the rays
			if (!_elt._rays.empty()) {
				std::list<std::pair<int,double> >::iterator it =  _elt._rays.begin();
				for (; it != _elt._rays.end(); ++it) {
					intermediate = Interval (it -> second) * alpha;
					it -> second = intermediate.ub(); // Check if it true
				}
			}

			_elt._garbage += roundoff_error * Interval(-1,1);
			if (_elt._garbage.rad() > AF_fAFFullI::maTol) {
				std::pair<int,double> pcumul(AF_fAFFullI::_counter++, _elt._garbage.rad());
				_elt._rays.push_back(pcumul);
				_elt._garbage = Interval(0.0);
			}

			_elt._rays.remove_if(noise_null);
		} else {
			*this = itv()*alpha;

		}
	} else {  //scalar alpha
		*this = itv()* alpha;
	}
	return *this;
}



template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator+=(double beta) {

	if ((_actif==1) && (fabs(beta))<POS_INFINITY) {
		Interval roundoff_error(0.,0.);
		Interval intermediate(0.,0.);
		intermediate = Interval(_elt._center) + beta;
		_elt._center = intermediate.mid();
		roundoff_error += intermediate.rad();

		_elt._garbage += roundoff_error * Interval(-1,1);
		if (_elt._garbage.rad() > AF_fAFFullI::maTol) {
			std::pair<int,double> pcumul(AF_fAFFullI::_counter++, _elt._garbage.rad());
			_elt._rays.push_back(pcumul);
			_elt._garbage = Interval(0.0);
		}

		_elt._rays.remove_if(noise_null);
	} else {
		*this = itv()+ beta;
	}
	return *this;
}




template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::inflate(double ddelta) {
	assert(ddelta>=0);
	if (ddelta>0) {
		if (is_actif()) {
			if ((ddelta)<POS_INFINITY) {
				_actif=1;
				std::pair<int,double> pdelta(AF_fAFFullI::_counter++, ddelta);
				_elt._rays.push_back(pdelta);
				_elt._rays.remove_if(noise_null);
			}
			else {
				*this = Interval::all_reals();
			}
		} else {
			*this = itv()+Interval(-1,1)*ddelta;
		}
	}
	return *this;
}



template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator+=( const AffineMain<AF_fAFFullI>& y) {

	if (is_actif() && y.is_actif()) {
		if (y.is_degenerated()) {
			*this += y.mid();
		} else if (is_degenerated()) {
			double tmp = _elt._center;
			*this = y;
			*this += tmp;
		} else {
			Interval roundoff_error(0.,0.);
			Interval intermediate(0.,0.);
			// Computation step for the center
			intermediate = Interval(_elt._center) + y._elt._center;
			_elt._center = intermediate.mid();
			roundoff_error += intermediate.rad();

			_elt._garbage += y._elt._garbage;

			// Computation step for the rays
			if (_elt._rays.empty())	{
				if (!y._elt._rays.empty())	{
					std::list<std::pair<int,double> >::const_iterator ity = y._elt._rays.begin();
					for (; ity != y._elt._rays.end(); ++ity)	{
						_elt._rays.push_back(std::pair<int,double>(ity->first,ity->second));
					}
				}
			} else if (/*!_elt._rays.empty() &&*/ !y._elt._rays.empty()) {
				std::list<std::pair<int,double> >::iterator it =  _elt._rays.begin();
				std::list<std::pair<int,double> >::const_iterator ity =  y._elt._rays.begin();

				while ((ity != y._elt._rays.end()) || (it != _elt._rays.end()))		{
					if (ity == y._elt._rays.end()) { //y is finished : stop
						break;
					}
					else if (it == _elt._rays.end()) { //x is finished : we push y
						std::pair<int,double> py(ity->first,ity->second );
						_elt._rays.insert(it, py);
						ity++;
					}
					else if (it -> first == ity -> first) { //same noise term : add
						intermediate = Interval(it -> second) + ity -> second;
						it -> second = intermediate.ub();
						it++;
						ity++;
					}
					else if (it -> first < ity -> first) { //noise of y after current x noise : x++
						it++;
					}
					else  { //noise of y before current x noise : add y before x
						std::pair<int,double> py(ity->first,ity->second );
						_elt._rays.insert(it, py);
						ity++;
					}
				}
			}

			_elt._garbage += roundoff_error * Interval(-1,1);
			if (_elt._garbage.rad() > AF_fAFFullI::maTol) {
				std::pair<int,double> pcumul(AF_fAFFullI::_counter++, _elt._garbage.rad());
				_elt._rays.push_back(pcumul);
				_elt._garbage = Interval(0.0);
			}

			_elt._rays.remove_if(noise_null);
		}
	} else if (is_actif()) { // y is not a valid Affine Form. So we add y.itv() such as an interval
		*this += y.itv();
	} else if (y.is_actif()) { // this case is perhap not necessary
		Interval tmp = itv();
		*this = y;
		*this += tmp;
	} else {
		*this = itv() + y.itv();
	}


	return *this;
}





// classical form
template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator*=(const AffineMain<AF_fAFFullI>& y) {

	if (is_actif() && (y.is_actif())) {
		if (y.is_degenerated()) {
			*this *= y._elt._center;
		} else if (is_degenerated()) {
			double tmp = _elt._center;
			*this = y;
			*this *= tmp;
		} else  {
			_actif =1;
			AffineMain<AF_fAFFullI> ax(*this);

			double y0 = y._elt._center;
			ax._elt._center = 0.0;

			*this = (((y*(this->_elt._center)) + (ax*y0) )+ ( ax.itv() * (y.itv() - y._elt._center)));
		}
	} else {// y is not a valid Affine Form. So we add y.itv() such as an interval
		*this = itv() * y.itv();
	}

	return *this;
}

template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator*=(const Interval& y) {

	//*this *= Affine2Main<AF_fAFFullI>(y);
	if (	(!is_actif())||
			y.is_empty()||
			y.is_unbounded() ) {
		*this = itv()*y;

	} else if (y.is_degenerated()) {
		*this *= y.mid();
	} else {
		 AffineMain<AF_fAFFullI> tmp;
		 tmp = y;
		*this *= tmp;
	}
	return *this;
}





template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator+=(const Interval& y) {

	if (	(!is_actif())||
			y.is_empty()||
			y.is_unbounded() ) {
		*this = (itv()+y);

	} else {
		 AffineMain<AF_fAFFullI> tmp;
		 tmp = y;
		*this += tmp;

	}
	return *this;
}




template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::Asqr(const Interval& itv) {

	if (	(!is_actif())||
			itv.is_empty()||
			itv.is_unbounded()||
			(itv.diam() < AF_EC)  ) {
		*this = pow(itv,2);

	} else  {
		// TODO il faut faire this->Apow(int n, const Interval itv);
		// TODO ou tester AffineMain<AF_fAFFullI> y(*this); *this *=y;
		if (itv.lb() >0 ) {
			this->Apow(Interval(2),itv);
		} else if (itv.ub()<0) {
			this->Aneg();
			this->Apow(Interval(2),-itv);
		} else {
			this->Aabs(abs(itv));
			this->Apow(Interval(2),abs(itv));
		}
	}

	return *this;
}




template<>
void AffineMain<AF_fAFFullI>::compact(double tol){
	Interval pmOne(-1.0, 1.0);

	if (!_elt._rays.empty()) {
		std::list<std::pair<int,double> >::iterator iter = _elt._rays.begin();
		Interval cumul(0.0);

		int list_size = _elt._rays.size();

		double seuil = tol * list_size * 0.0001;

		while (list_size > 10)
		{

			for (iter = _elt._rays.begin(); iter != _elt._rays.end(); ++iter) {
				if (std::abs(iter -> second) < seuil) {
					cumul += (iter->second * pmOne);

					iter = _elt._rays.erase(iter);
					iter--;
				}
			}
			seuil*=10;
			list_size = _elt._rays.size();
		}

		_elt._garbage += cumul;

		if (_elt._garbage.rad() > AF_fAFFullI::maTol)
		{
			std::pair<int,double> pcumul(AF_fAFFullI::_counter++, _elt._garbage.rad());
			_elt._rays.push_back(pcumul);
			_elt._garbage = Interval(0.0);
		}
	}
	return;
}


template<>
std::ostream& operator<<(std::ostream& os, const AffineMain<AF_fAFFullI>& x) {
	os << x.itv() << " : ";
	if (x.is_actif()) {
		os << x.mid();

		if (!x._elt._rays.empty()) {
			std::list<std::pair<int,double> >::const_iterator iter = x._elt._rays.begin();
			for (; iter != x._elt._rays.end(); ++iter) {
				double v = iter -> second;
				if (v!=0)
					os << " + " << v  << " eps_" << iter -> first;
			}
		}
		os << " + " << x.err() << "[-1,1]";
	} else {
		os << "Affine3 Form not activate ";
	}
	return os;
}


template<>
void AffineMain<AF_fAFFullI>::resize(int n) { assert(n>=1); }





//===========================================================================================
//===========================================================================================






}// end namespace ibex









