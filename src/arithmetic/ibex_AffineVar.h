/* ============================================================================
 * I B E X - AffineVar definition
 * ============================================================================
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Jordan Ninin
 * Bug fixes   :
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */

#ifndef IBEX_AFFINEVAR_H_
#define IBEX_AFFINEVAR_H_

#include "ibex_IntervalVector.h"
#include "ibex_Affine.h"
#include <math.h>
#include <ostream>
#include <cassert>
#include "ibex_Exception.h"

namespace ibex {

template<class T>  class AffineVarMainVector;

/**
 * \ingroup arithmetic
 *
 * \brief Affine Arithmetic AF2
 *
 * Code for the particular case:
 * if the affine form is actif, _actif=1  and _n is the size of the affine form
 * if the set is degenerate, _actif = 0 and itv().diam()< AF_EC
 * if the set is empty, _actif = -1
 * if the set is ]-oo,+oo[, _actif = -2 and _err =]-oo,+oo[
 * if the set is [a, +oo[ , _actif = -3 and _err = [a, +oo[
 * if the set is ]-oo, a] , _actif = -4 and _err = ]-oo, a]
 *
 */


typedef AffineVarMainVector<AF_Default> Affine2Vars;
typedef AffineVarMainVector<AF_Other>  Affine3Vars;



template<class T=AF_Default>
class AffineVarMain : public AffineMain<T> {


public:
	/** \brief Set *this to itv.
	 */
	AffineVarMain& operator=(const Interval& itv);

private:
	friend class AffineVarMainVector<T>;
//	friend class AffineMainMatrix<T>;


    unsigned long int _var;

    AffineVarMain() : _var(0) {};

	/** \brief Create an affine form by copy. */
    AffineVarMain(const AffineVarMain<T>& x);

	/** \brief Create an affine form with at most \size variables and  initialized the \var^th variable with  itv. */
	AffineVarMain(int size, int var, const Interval& itv);

	/** \brief Copy an affine Var form. */
	AffineVarMain& operator=(const AffineVarMain<T>& x) ;

	AffineVarMain& operator+=(const Vector& x2) {ibex_error(" AffineVarMain : operator+= non valid");};
	AffineVarMain& operator+=(const Interval& x2) {ibex_error(" AffineVarMain : operator+= non valid");};
	AffineVarMain& operator+=(const AffineMain<T>& x2) {ibex_error(" AffineVarMain : operator+= non valid");};
	AffineVarMain& operator-=(const Vector& x2) {ibex_error(" AffineVarMain : operator-= non valid");};
	AffineVarMain& operator-=(const Interval& x2) {ibex_error(" AffineVarMain : operator-= non valid");};
	AffineVarMain& operator-=(const AffineMain<T>& x2) {ibex_error(" AffineVarMain : operator-= non valid");};
	AffineVarMain& operator*=(double d) {ibex_error(" AffineVarMain : operator*= non valid");};
	AffineVarMain& operator*=(const Interval& x1) {ibex_error(" AffineVarMain : operator*= non valid");};
	AffineVarMain& operator*=(const AffineMain<T>& x1) {ibex_error(" AffineVarMain : operator*= non valid");};

	AffineVarMain& Asqr(const Interval& itv)  {ibex_error(" AffineVarMain : operator Asqr non valid");};
	AffineVarMain&  Aneg()  {ibex_error(" AffineVarMain : operator Aneg non valid");};
	AffineVarMain&  Ainv(const Interval& itv)  {ibex_error(" AffineVarMain : operator Ainv non valid");};
	AffineVarMain&  Asqrt(const Interval& itv)  {ibex_error(" AffineVarMain : operator non Asqrt valid");};
	AffineVarMain&  Aexp(const Interval& itv)  {ibex_error(" AffineVarMain : operator Aexp non valid");};
	AffineVarMain&  Alog(const Interval& itv)  {ibex_error(" AffineVarMain : operator Alog non valid");};
	AffineVarMain&  Apow(int n, const Interval& itv)  {ibex_error(" AffineVarMain : operator Apow non valid");};
	AffineVarMain&  Apow(double d, const Interval& itv)  {ibex_error(" AffineVarMain : operator Apow non valid");};
	AffineVarMain&  Apow(const Interval &y, const Interval& itvx)  {ibex_error(" AffineVarMain : operator Apow non valid");};
	AffineVarMain&  Aroot(int n, const Interval& itv)  {ibex_error(" AffineVarMain : operator Aroot non valid");};
	AffineVarMain&  Acos(const Interval& itv)  {ibex_error(" AffineVarMain : operator Acos non valid");};
	AffineVarMain&  Asin(const Interval& itv)  {ibex_error(" AffineVarMain : operator Asin non valid");};
	AffineVarMain&  Atan(const Interval& itv)  {ibex_error(" AffineVarMain : operator Atan non valid");};
	AffineVarMain&  Aacos(const Interval& itv)  {ibex_error(" AffineVarMain : operator Aacos non valid");};
	AffineVarMain&  Aasin(const Interval& itv)  {ibex_error(" AffineVarMain : operator Aasin non valid");};
	AffineVarMain&  Aatan(const Interval& itv)  {ibex_error(" AffineVarMain : operator Aatan non valid");};
	AffineVarMain&  Acosh(const Interval& itv)  {ibex_error(" AffineVarMain : operator Acosh non valid");};
	AffineVarMain&  Asinh(const Interval& itv)  {ibex_error(" AffineVarMain : operator Asinh non valid");};
	AffineVarMain&  Atanh(const Interval& itv)  {ibex_error(" AffineVarMain : operator Atanh non valid");};
	AffineVarMain&  Aabs(const Interval& itv)  {ibex_error(" AffineVarMain : operator Aabs non valid");};
	AffineVarMain&  Ainv_CH(const Interval& itv)  {ibex_error(" AffineVarMain : operator Ainv_CH non valid");};
	AffineVarMain&  Asqrt_CH(const Interval& itv)  {ibex_error(" AffineVarMain : operator Asqrt_CH non valid");};
	AffineVarMain&  Aexp_CH(const Interval& itv)  {ibex_error(" AffineVarMain : operator Aexp_CH non valid");};
	AffineVarMain&  Alog_CH(const Interval& itv)  {ibex_error(" AffineVarMain : operator Alog_CH non valid");};
	AffineVarMain&  Ainv_MR(const Interval& itv)  {ibex_error(" AffineVarMain : operator Ainv_CH non valid");};
	AffineVarMain&  Asqrt_MR(const Interval& itv)  {ibex_error(" AffineVarMain : operator Asqrt_MR non valid");};
	AffineVarMain&  Aexp_MR(const Interval& itv)  {ibex_error(" AffineVarMain : operator Aexp_MR non valid");};
	AffineVarMain&  Alog_MR(const Interval& itv)  {ibex_error(" AffineVarMain : operator Alog_MR non valid");};
};



template<class T>
AffineVarMain<T>::AffineVarMain(int size, int var1, const Interval& itv) :
		AffineMain<T>(size, var1, itv),
		_var		(var1) {
}



template<class T>
AffineVarMain<T>::AffineVarMain(const AffineVarMain& x) :
		AffineMain<T>(x.size(), x._var, x.itv()),
		_var		(x._var) {
}


//===============================================================================================





template<class T=AF_Default>
class AffineVarMainVector : public AffineMainVector<T> {


//private:

//	AffineVarMainVector() : //_n(0),
//		AffineMainVector<T>(0) { }
//	Array<AffineVarMain<T> > _vec;	   // vector of elements

public:

	/** \brief  Create \a n Affine form . All the components are Affine([-oo,+oo])
	 * \pre n>0
	 */
	explicit AffineVarMainVector(int n);

	/**
	 * \brief  Create an AffineVarMainVector of dimension \a n with
	 * all the components initialized to Affine(\a n,i,\a x) }.
	 * \pre n>0
	 */
	explicit AffineVarMainVector(int n, const Interval& x);

	/**
	 * \brief  Create \a n AffineVarMainVector of dimension \a n with
	 * all the components initialized to \a x.
	 * \pre n>0
	 */
//	explicit AffineVarMainVector(int n, const AffineMain<T>& x);

	/**
	 * \brief Create  a copy of \a x.
	 */
	AffineVarMainVector(const AffineVarMainVector& x);

	/**
	 * \brief Create \a n AffineVarMainVector  initialized by
	 * AffineMain(\a n, i,Interval(bounds[i][0],bounds[i][1]) )
	 * \param bounds an nx2 array of doubles
	 * \pre n>0
	 */
	explicit AffineVarMainVector(int n, double  bounds[][2]);

	/**
	 * \brief Create \a x.size AffineVarMainVector of dimension \a x.size with
	 * the [i] component initialized to
	 * if !(\a b) Affine(x[i])
	 * else  Affine(x.size(), i,x[i])
	 */
	explicit AffineVarMainVector(const IntervalVector& x);

	/**
	 * \brief Create the degenerated AffineVarMainVector x
	 *
	 */
	explicit AffineVarMainVector(const Vector& x);


	/**
	 * \brief Delete this vector
	 */
	//virtual ~AffineVarMainVector();

	/**
	 * \brief Resize this AffineMainVector.
	 *
	 * If the size is increased, the existing components are not
	 * modified and the new ones are set to (-inf,+inf), even if
	 * (*this) is the empty Interval (however, in this case, the status of
	 * (*this) remains "empty").
	 */
	void resize(int n2);

private:
	void put(int start_index, const AffineMainVector<T>& subvec) {ibex_error(" AffineVarMainVector : operator put non valid");};
	AffineVarMainVector& operator=(const AffineMainVector<T>& x) {ibex_error(" AffineVarMainVector : operator= non valid");};
	AffineVarMainVector& operator+=(const Vector& x2) {ibex_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator+=(const IntervalVector& x2) {ibex_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator+=(const AffineMainVector<T>& x2) {ibex_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator-=(const Vector& x2) {ibex_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator-=(const IntervalVector& x2) {ibex_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator-=(const AffineMainVector<T>& x2) {ibex_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator*=(double d) {ibex_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator*=(const Interval& x1) {ibex_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator*=(const AffineMain<T>& x1) {ibex_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator/=(double d) {ibex_error(" AffineVarMainVector : operator/= non valid");};
	AffineVarMainVector& operator/=(const Interval& x1) {ibex_error(" AffineVarMainVector : operator/= non valid");};
	AffineVarMainVector& operator/=(const AffineMain<T>& x1) {ibex_error(" AffineVarMainVector : operator/= non valid");};

};




//===============================================================================================

//template<class T>
//inline AffineVarMainVector<T>::~AffineVarMainVector<T>() {
//	for (int i=0; i<_vec.size(); i++) {
//		delete _vec[i];
//	}
//}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n) {
	assert(n>=1);
	this->_n   = n;
	this->_vec = new AffineVarMain<T>[n];
	assert(n>=1);
	for (int i = 0; i < n; i++){
		(this->_vec)[i] = AffineVarMain<T>(n, i, Interval::all_reals());
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n, const Interval& x) {
	assert(n>=1);
	this->_n   = n;
	this->_vec  = new AffineVarMain<T>[n];
	for (int i = 0; i < n; i++) {
		(this->_vec)[i] = AffineVarMain<T>(n, i, x);
	}
}


template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const AffineVarMainVector<T>& x) {
	this->_n   = x.size();
	this->_vec  = new AffineVarMain<T>[x.size()];
	for (int i = 0; i < x.size(); i++){
		(this->_vec)[i] = AffineVarMain<T>(x.size(),i,(x[i]).itv());
	}

}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n, double bounds[][2])  {
	assert(n>=1);
	this->_n   = n;
	this->_vec  = new AffineVarMain<T>[n];
	if (bounds == 0){ // probably, the user called AffineVarMainVector<T>(n,0) and 0 is interpreted as NULL!
		for (int i = 0; i < n; i++){
			(this->_vec)[i] = AffineVarMain<T>(n, i, Interval(0));
		}
	}
	else {
		for (int i = 0; i < n; i++){
			(this->_vec)[i] = AffineVarMain<T>(n, i, Interval(bounds[i][0], bounds[i][1]));
		}

	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const IntervalVector& x) {
	this->_n   = x.size();
	this->_vec  = new AffineVarMain<T>[x.size()];
	for (int i = 0; i < x.size(); i++){
		(this->_vec)[i] = AffineVarMain<T>(x.size(), i, x[i]);
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const Vector& x) {
	this->_n   = x.size();
	this->_vec.resize(x.size());
	for (int i = 0; i < x.size(); i++){
		(this->_vec)[i] = AffineVarMain<T>(x.size(), i, Interval(x[i]));
	}
}

template<class T>
void AffineVarMainVector<T>::resize(int n2) {
	assert(n2>=1);
	assert((this->_vec==NULL && this->_n==0) || (this->_n!=0 && this->_vec!=NULL));

	if (n2==this->size()) return;

	AffineVarMain<T>* newVec=new AffineVarMain<T>[n2];
	int i=0;
	for (; i<this->size() && i<n2; i++){
		newVec[i]=AffineVarMain<T>(n2, i, (this->_vec[i]).itv());
	}
	for (; i<n2; i++) {
		newVec[i]=AffineVarMain<T>(n2, i, Interval::all_reals());
	}

	if (this->_vec!=NULL) {// vec==NULL happens when default constructor is used (n==0)
		delete[] this->_vec;
	}

	this->_n   = n2;
	this->_vec = newVec;
}



}


#endif /* IBEX_Affine_H_ */


/** \brief atan2(AF[y],AF[x]). */
//Affine2 atan2(const Affine2& y, const Affine2& x);
/** \brief atan2([y],AF[x]). */
//Affine2 atan2(const Interval& y, const Affine2& x);
/** \brief atan2(AF[y],[x]). */
//Affine2 atan2(const Affine2& y, const Interval& x);
/** \brief cosh(AF[x]). */
