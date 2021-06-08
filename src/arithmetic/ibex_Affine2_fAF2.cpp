/* ============================================================================
 * I B E X - Implementation of the AffineMain<AF_fAF2> class based on fAF version 2
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Jordan Ninin
 * Created     : Jul 16, 2013
 * ---------------------------------------------------------------------------- */


#ifndef IBEX_AFFINE_AF2_H_
#define IBEX_AFFINE_AF2_H_


#include "ibex_Affine.h"
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
template<>
AffineMain<AF_fAF2>::AffineMain() :
		 _actif (-2     ),
		 _n		(0		),
		 _elt	(NULL	,POS_INFINITY)	{
 }


template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator=(const Interval& x) {

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
		if (_elt._val==NULL) { _elt._val = new double[_n+1]; }
		_elt._val[0] = x.mid();
		for (int i=1; i<=_n;i++) {
			_elt._val[i] =0;
		}
		if ( x.is_degenerated()) {
			_actif=0;
			_elt._err	= 0;
		} else {
			_actif = 1;
			_elt._err	= x.rad();
		}
	}
	return *this;
}


template<>
AffineMain<AF_fAF2>::AffineMain(int size, int var, const Interval& itv) :
			_actif	(0),
			_n 		(size),
			_elt	(NULL,0.0)
{
	assert((size>=0) && (var>=0) && (var<=size));
	if (!(itv.is_unbounded()||itv.is_empty())) {
		_elt._val	=new double[size + 1];
		_elt._val[0] = itv.mid();
		for (int i = 1; i <= size; i++){
			_elt._val[i] = 0.0;
		}
		if (! itv.is_degenerated()) {
			_actif =1;
			_elt._val[var+1] = itv.rad();
		}
	} else {
		*this = itv;
	}
}



template<>
AffineMain<AF_fAF2>::AffineMain(const AffineMain<AF_fAF2>& x) :
		_actif	(x._actif),
		_n		(x.size()),
		_elt	(NULL	,x._elt._err ) {
	if (is_actif()) {
		_elt._val =new double[x.size() + 1];
		for (int i = 0; i <= x.size(); i++){
			_elt._val[i] = x._elt._val[i];
		}
	}
}



template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator=(const AffineMain<AF_fAF2>& x) {
	if (this != &x) {
		_elt._err = x._elt._err;
		_actif = x._actif;
		_n =x._n;
		if (x.is_actif()) {
			if (_elt._val!=NULL) { delete[] _elt._val; }
			_elt._val = new double[size()+1];
			int i = 0;
			for (; i <= x.size(); i++) {
				_elt._val[i] = x._elt._val[i];
			}
		}
	}
	return *this;
}

template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator=(double d) {

	if (fabs(d)<POS_INFINITY) {
		_actif = 0;
		if (_elt._val==NULL) _elt._val = new double[size()+1];
		_elt._err = 0.0; //abs(d)*AF_EE;
		_elt._val[0] = d;
		for (int i = 1; i <= size(); i++){
			_elt._val[i] = 0.0;
		}
	} else {
		if (d>0) {
			_actif = -3;
			_elt._err = d;
		} else {
			_actif = -4;
			_elt._err = d;
		}
	}
	return *this;
}





template<>
double AffineMain<AF_fAF2>::val(int i) const{
	assert(is_actif() &&(0<=i) && (i<size()));
	return _elt._val[i+1];
}

template<>
double AffineMain<AF_fAF2>::err() const{
	return _elt._err;
}


template<>
const Interval AffineMain<AF_fAF2>::itv() const {

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
		return Interval(_elt._err,POS_INFINITY);
		break;
	}
	case -4: {
		return Interval(NEG_INFINITY,_elt._err);
		break;
	}
	case 0: {
		return Interval(_elt._val[0]);
		break;
	}
	default: { // _actif==1
		Interval res(_elt._val[0]);
		Interval pmOne(-1.0, 1.0);
		for (int i = 1; i <= size(); i++){
			res += (_elt._val[i] * pmOne);
		}
		res += _elt._err * pmOne;
		return res;
		break;
	}
	}

}


template<>
double AffineMain<AF_fAF2>::mid() const{
	return (is_actif())? _elt._val[0] : itv().mid();
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
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::Aneg() {
	switch(_actif) {
	case -3 : {
		_elt._err = -_elt._err;
		_actif    = -4;
		break;
	}
	case -4 : {
		_elt._err = -_elt._err;
		_actif    = -3;
		break;
	}
	case 0 :{
		_elt._val[0] = (-_elt._val[0]);
		break;
	}
	case 1 : {
		for (int i = 0; i <= size(); i++) {
			_elt._val[i] = (-_elt._val[i]);
		}
		break;
	}
	default :
		break;
	}

	return *this;
}




template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator*=(double alpha) {
	if (_actif==1) {  // multiply by a scalar alpha
		if (alpha==0.0) {
			_actif = 0;
			for (int i=0; i<=_n;i++) {
				_elt._val[i]=0;
			}
			_elt._err = 0;
		} else if ( fabs(alpha) < POS_INFINITY) {
			double temp, ttt, sss, eee;
			ttt= 0.0;
			sss= 0.0;
			for (int i=0; i<=size();i++) {
				eee = _elt.twoProd(_elt._val[i], alpha, &temp);
				_elt._val[i] = temp;
				ttt = (1+2*AF_EM)*(ttt+fabs(eee));
				if (fabs(_elt._val[i])<AF_EC) {
					sss = (1+2*AF_EM)*(sss+ fabs(_elt._val[i]));
					_elt._val[i] = 0.0;
				}
			}
			//	_elt._err = (1+2*AF_EM)*((1+2*AF_EM)*fabs(alpha)*_elt._err+AF_EE*AF_EM*ttt + AF_EE*sss);
			_elt._err = (1+2*AF_EM)*( ((1+2*AF_EM)*fabs(alpha)*_elt._err) +	((AF_EE*ttt) +	(AF_EE*sss)) );

			bool b = (_elt._err<POS_INFINITY);
			for (int i=0;i<=size();i++) {
				b &= (fabs(_elt._val[i])<POS_INFINITY);
			}
			if (!b) { *this = Interval::all_reals(); }

		} else {
			*this = itv()*alpha;
		}
	} else {  //scalar alpha
		*this = itv()* alpha;
	}
	return *this;
}



template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator+=(double beta) {

	if ((_actif==1) && (fabs(beta)<POS_INFINITY)) {
		double temp, ttt, sss, eee;
		ttt=0.0;
		sss=0.0;
		eee = _elt.twoSum(_elt._val[0],beta,&temp);
		ttt = (1+2*AF_EM)*(ttt+fabs(eee));
		if (fabs(temp)<AF_EC) {
			sss = (1+2*AF_EM)*(sss+fabs(temp));
			_elt._val[0] = 0.0;
		} else {
			_elt._val[0]=temp;
		}
		//				_elt._err = (1+2*AF_EM)*(_elt._err+ (AF_EE*(AF_EM*ttt)+AF_EE*sss));
		_elt._err = (1+2*AF_EM)*(_elt._err +	(AF_EE*(ttt)+ AF_EE*sss) );

		if (!(_elt._err<POS_INFINITY && (fabs(_elt._val[0])<POS_INFINITY))) { *this = Interval::ALL_REALS; }

	} else {
		*this = itv()+ beta;
	}
	return *this;

}




template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::inflate(double ddelta) {
	assert(ddelta>=0);
	if (ddelta>0) {
		if (is_actif()) {
			if ((ddelta)<POS_INFINITY) {
				_actif=1;
				double temp, ttt, sss, eee;
				ttt=0.0;
				sss=0.0;
				eee = _elt.twoSum(_elt._err,fabs(ddelta), &temp);
				ttt = (1+2*AF_EM)*(fabs(eee));
				if (fabs(temp)<AF_EC) {
					sss = (1+2*AF_EM)*(fabs(temp));
					temp =0;
				}
				//				_elt._err = (1+2*AF_EM)*(temp+ (AF_EE*(AF_EM*ttt)));;
				_elt._err = (1+2*AF_EM)*( temp + (AF_EE*(ttt) + AF_EE*sss) );

				if (!(_elt._err<POS_INFINITY)) { *this = Interval::all_reals(); }
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
void AffineMain<AF_fAF2>::resize(int n) {
	assert(n>=1);
	if (n>_n) {
		double * tmp= new double[n+1];
		int i=0;
		for (;i<=_n;i++) {
			tmp[i] = _elt._val[i];
		}
		for (;i<=n;i++) {
			tmp[i] = 0;
		}
		_n = n;
		delete[] _elt._val;
		_elt._val = tmp;
	} else if (n<_n) {
		double * tmp= new double[n+1];
		int i=0;
		for (;i<=n;i++) {
			tmp[i] = _elt._val[i];
		}
		_n = n;
		delete[] _elt._val;
		_elt._val = tmp;
	}
}





template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator+=(const AffineMain<AF_fAF2>& y) {

	if (is_actif() && y.is_actif()) {
		if (y.is_degenerated()) {
			*this += y._elt._val[0];
		} else if (is_degenerated()) {
			double tmp = _elt._val[0];
			*this = y;
			*this += tmp;
		} else {
			if (_n < y.size()) {
				resize(y.size());
			}
			double temp, ttt, sss, eee;
			ttt=0.0;
			sss=0.0;
			for(int i=0;i<=y.size();i++) {
				eee = _elt.twoSum(_elt._val[i], y._elt._val[i], &temp);
				ttt = (1+2*AF_EM)*(ttt+fabs(eee));
				if (fabs(temp)<AF_EC) {
					sss = (1+2*AF_EM)*(sss+ fabs(temp));
					_elt._val[i] = 0.0;
				}
				else {
					_elt._val[i]=temp;
				}
			}
			// _elt._err = (1+2*AF_EM)*((_elt._err+y._elt._err+ (AF_EE*(AF_EM*ttt)+AF_EE*sss));
			_elt._err = (1+2*AF_EM)*( (_elt._err+y._elt._err) + ((AF_EE*(ttt)) + (AF_EE*sss)) );

			bool b = (_elt._err<POS_INFINITY);
			for (int i=0;i<=_n;i++) {
				b &= (fabs(_elt._val[i])<POS_INFINITY);
			}
			if (!b) {*this = Interval::all_reals(); }

		}
	} else if (is_actif()) { // y is not a valid affine form. So we add y.itv() such as an interval
		*this += y.itv();
	} else if (y.is_actif()) {// *this is not a valid affine form. We inverse *this and y.
		Interval tmp = itv();
		*this = y;
		*this += tmp;
	} else {
		*this = itv() + y.itv();
	}
	//std::cout << "OUT += "<<std::endl<< *this  << std::endl;
	return *this;
}


/**
 * see  Equation(17)  of
 * X.-H. Vu, D. Sam-Haroud, and B. Faltings. Combining multiple inclusion representa-
tions in numerical constraint propagation. In Tools with Artificial Intelligence, IEEE
International Conference on, pages 458–467, Los Alamitos, CA, USA, 2004. IEEE Com-
puter Society.
 */
template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator*=(const AffineMain<AF_fAF2>& y) {

	if (is_actif() && (y.is_actif())) {
		if (y.is_degenerated()) {
			*this *= y._elt._val[0];
		}	else if (is_degenerated()) {
			double tmp = _elt._val[0];
			*this = y;
			*this *= tmp;
		} else 	 {
			if (_n < y.size()) {
				resize(y.size());
			}
			double Sx, Sy, Sxy, Sz, ttt, sss, ppp, tmp, xVal0, eee;
			double * xTmp;

			xTmp = new double[_n + 1];
			Sx=0.0; Sy=0.0; Sxy=0.0; Sz=0.0; ttt=0.0; sss=0.0; ppp=0.0; tmp=0.0; xVal0=0.0; eee=0.0;

			for (int i = 1; i <= _n; i++) {
				ppp=0;
				if (i<=y.size()) {
					eee = _elt.twoProd(_elt._val[i],y._elt._val[i], &ppp);
					ttt = (1+2*AF_EM)*(ttt+fabs(eee));

					eee = _elt.twoSum(Sz,ppp, &tmp);
					ttt = (1+2*AF_EM)*(ttt+fabs(eee));
					Sz = tmp;

					if (fabs(Sz) < AF_EC) {
						sss = (1+2*AF_EM)*(sss+ fabs(Sz));
						Sz = 0.0;
					}

					eee = _elt.twoSum(Sxy,fabs(ppp), &tmp);
					ttt = (1+2*AF_EM)*(ttt+fabs(eee));
					Sxy = tmp;

					if (fabs(Sxy) < AF_EC) {
						sss = (1+2*AF_EM)*(sss+ fabs(Sxy));
						Sxy = 0.0;
					}
				}
				eee = _elt.twoSum(Sx,fabs(_elt._val[i]), &tmp);
				ttt = (1+2*AF_EM)*(ttt+fabs(eee));
				Sx = tmp;

				if (fabs(Sx) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ fabs(Sx));
					Sx = 0.0;
				}

				if (i<=y.size()) {
					eee = _elt.twoSum(Sy,fabs(y._elt._val[i]), &tmp);
					ttt = (1+2*AF_EM)*(ttt+fabs(eee));
					Sy = tmp;

					if (fabs(Sy) < AF_EC) {
						sss = (1+2*AF_EM)*(sss+ fabs(Sy));
						Sy = 0.0;
					}
				}

			}

			xVal0 = _elt._val[0];
			// RES = X%T(0) * res
			for (int i = 0; i <= _n; i++) {
				eee = _elt.twoProd(_elt._val[i],y._elt._val[0], &ppp);
				ttt = (1+2*AF_EM)*(ttt+fabs(eee));
				_elt._val[i] = ppp;

				if (fabs(_elt._val[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ fabs(_elt._val[i]));
					_elt._val[i] = 0.0;
				}
			}

			// Xtmp = X%T(0) * Y
			xTmp[0] = 0.0;
			for (int i = 1; i <= y.size(); i++) {
				eee = _elt.twoProd(xVal0,y._elt._val[i], &ppp);
				ttt = (1+2*AF_EM)*(ttt+fabs(eee));
				xTmp[i] = ppp;

				if (fabs(xTmp[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ fabs(xTmp[i]));
					xTmp[i] = 0.0;
				}

			}

			//RES =  RES + Xtmp = ( Y%(0) * X ) + ( X%T(0) * Y - X%T(0)*Y%(0) )
			for (int i = 0; i <= y.size(); i++) {

				eee = _elt.twoSum(_elt._val[i],xTmp[i], &tmp);
				ttt = (1+2*AF_EM)*(ttt+fabs(eee));
				_elt._val[i] = tmp;

				if (fabs(_elt._val[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ fabs(_elt._val[i]));
					_elt._val[i] = 0.0;
				}

			}

			eee = _elt.twoProd(0.5,Sz, &ppp);
			ttt = (1+2*AF_EM)*(ttt+fabs(eee));

			eee = _elt.twoSum(_elt._val[0],ppp, &tmp);
			ttt = (1+2*AF_EM)*(ttt+fabs(eee));
			_elt._val[0] = tmp;

			if (fabs(_elt._val[0]) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ fabs(_elt._val[0]));
				_elt._val[0] = 0.0;
			}

			eee = _elt.twoSum(_elt._err,Sx, &tmp);
			ttt = (1+2*AF_EM)*(ttt+fabs(eee));

			eee = _elt.twoSum(y._elt._err,Sy, &ppp);
			ttt = (1+2*AF_EM)*(ttt+fabs(eee));


			_elt._err = (1+ 2*AF_EM) * (
					((1+ 2*AF_EM) *fabs(y._elt._val[0]) * _elt._err)  +
					((1+ 2*AF_EM) *fabs(xVal0) * y._elt._err)  +
					((1+ 2*AF_EM) *(tmp * ppp)) +
					((1- 2*AF_EM) *(-0.5) *  Sxy)  +
					//					(AF_EE * (AF_EM * ttt))  +
					(AF_EE * (ttt))  +
					(AF_EE * sss)
			);


			bool b = (_elt._err<POS_INFINITY);
			for (int i=0;i<=_n;i++) {
				b &= (fabs(_elt._val[i])<POS_INFINITY);
			}
			if (!b) *this = Interval::all_reals();

			delete[] xTmp;
		}

	} else { // y or x is not a valid affine form. So we add y.itv() such as an interval
		*this = (itv() * y.itv());
	}

	//std::cout << "OUT *= "<<std::endl<< *this << std::endl;
	return *this;
}



template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator*=(const Interval& y) {
	if (	(!is_actif())||
			y.is_empty()||
			y.is_unbounded() ) {
		*this = itv()*y;
	} else {
		if (y.is_degenerated()) {
			*this *= y.mid();
		} else {
			AffineMain<AF_fAF2> tmp;
			tmp= y.mid();	 // to check if it is the best way to do it
			tmp.inflate(y.rad());
			*this *= tmp;
		}
	}
	return *this;
}




template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::Asqr(const Interval& itv) {

	if (	(!is_actif())||
			itv.is_empty()||
			itv.is_unbounded()||
			(itv.diam() < AF_EC)  ) {
		*this = pow(itv,2);

	} else  {

		double Sx, Sx2, ttt, sss, ppp, x0, eee,tmp;
		Sx = 0; Sx2 = 0; ttt = 0; sss = 0; ppp = 0; x0 = 0; eee =0.0; tmp =0.0;

		// compute the error
		for (int i = 1; i <= _n; i++) {

			eee = _elt.twoProd(_elt._val[i],_elt._val[i], &ppp);
			ttt = (1+2*AF_EM)*(ttt+fabs(eee));

			eee = _elt.twoSum(Sx2,ppp, &tmp);
			ttt = (1+2*AF_EM)*(ttt+fabs(eee));
			Sx2 = tmp;

			if (fabs(Sx2) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ fabs(Sx2));
				Sx2 = 0.0;
			}

			eee = _elt.twoSum(Sx,fabs(_elt._val[i]), &tmp);
			ttt = (1+2*AF_EM)*(ttt+fabs(eee));
			Sx = tmp;

			if (fabs(Sx) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ fabs(Sx));
				Sx = 0.0;
			}

		}
		// compute 2*_elt._val[0]*(*this)
		x0 = _elt._val[0];

		eee = _elt.twoProd(x0,x0, &ppp);
		ttt = (1+2*AF_EM)*(ttt+fabs(eee));
		_elt._val[0] = ppp;

		if (fabs(_elt._val[0]) < AF_EC) {
			sss = (1+2*AF_EM)*(sss+ fabs(_elt._val[0]));
			_elt._val[0] = 0.0;
		}

		// compute 2*_elt._val[0]*(*this)
		for (int i = 1; i <= _n; i++) {

			eee = _elt.twoProd((2*x0),_elt._val[i], &ppp);
			ttt = (1+2*AF_EM)*(ttt+fabs(eee));
			_elt._val[i] = ppp;

			if (fabs(_elt._val[i]) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ fabs(_elt._val[i]));
				_elt._val[i] = 0.0;
			}

		}

		eee = _elt.twoProd(0.5,Sx2, &ppp);
		ttt = (1+2*AF_EM)*(ttt+fabs(eee));

		eee = _elt.twoSum(_elt._val[0],ppp, &tmp);
		ttt = (1+2*AF_EM)*(ttt+fabs(eee));
		_elt._val[0] = tmp;

		if (fabs(_elt._val[0]) < AF_EC) {
			sss = (1+2*AF_EM)*(sss+ fabs(_elt._val[0]));
			_elt._val[0] = 0.0;
		}

		eee = _elt.twoSum(_elt._err,Sx, &tmp);
		ttt = (1+2*AF_EM)*(ttt+fabs(eee));

		_elt._err = (1+ 2*AF_EM) * (
				((1+ 2*AF_EM) *2*fabs(x0) * _elt._err)  +
				((1+ 2*AF_EM) *(tmp * tmp)) +
				((1- 2*AF_EM) *(-0.5) *  Sx2)  +
//					(AF_EE * (AF_EM * ttt))  +
				(AF_EE * (ttt))  +
				(AF_EE * sss)
				);

		{
			bool b = (_elt._err<POS_INFINITY);
			for (int i=0;i<=_n;i++) {
				b &= (fabs(_elt._val[i])<POS_INFINITY);
			}
			if (!b) {
				*this = Interval::all_reals();
			}
		}

	}

//	std::cout << "out sqr "<<std::endl;
	return *this;
}


template<>
void AffineMain<AF_fAF2>::compact(double tol){
	for (int i=1;i<=_n;i++) {
		if (fabs(_elt._val[i])<tol) {
			double temp=0.0;
			double sss=0.0;
			double eee = _elt.twoSum(_elt._err,fabs(_elt._val[i]), &temp);
			double ttt = (1+2*AF_EM)*(fabs(eee));
			if (fabs(temp)<AF_EC) {
				sss = (1+2*AF_EM)*(fabs(temp));
				temp =0;
			}
//			_elt._err = (1+2*AF_EM)*(temp+ (AF_EE*(AF_EM*ttt)));;
			_elt._err = (1+2*AF_EM)*( temp + (AF_EE*(ttt) + AF_EE*sss) );

			_elt._val[i] =0;
		}
	}
}






//===========================================================================================
//===========================================================================================



}// end namespace ibex


#endif /* IBEX_Affine_AF2_H_ */


