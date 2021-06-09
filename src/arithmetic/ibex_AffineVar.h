/* ============================================================================
 * I B E X - AffineVar definition
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENCE.
 *
 * Author(s)   : Jordan Ninin
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */

#ifndef IBEX_AFFINEVAR_H_
#define IBEX_AFFINEVAR_H_

#include "ibex_TemplateVector.h"
#include "ibex_IntervalVector.h"
#include "ibex_AffineMain.h"
#include "ibex_AffineVector.h"
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


typedef AffineVarMainVector<AF_Default> Affine2Variables;
typedef AffineVarMainVector<AF_Other>  Affine3Variables;



template<class T=AF_Default>
class AffineVarMain : public AffineMain<T> {


public:
	/** \brief Set *this to itv.
	 */
	AffineVarMain& operator=(const Interval& itv);

	/** \brief Return -*this. */
	AffineMain<T> operator-() const;


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

/** \brief Return (-x) */
template<class T>
inline AffineMain<T> AffineVarMain<T>::operator-() const {
	AffineMain<T> res(*this);
	res.Aneg();
	return res;
}


//===============================================================================================





template<class T=AF_Default>
class AffineVarMainVector {

	int _n;						// size of _vec
	AffineVarMain<T>* _vec;	   // vector of elements

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
	virtual ~AffineVarMainVector() {delete [] _vec;};


	/**
	 * \brief Set this AffineMainVector to the empty AffineMainVector
	 *
	 * The dimension remains the same.
	 */
	void set_empty();

	/**
	 * \brief Set all the elements to 0 (even if empty).
	 *
	 * \note Emptiness is "overridden".
	 */
	void clear();

	/**
	 * \brief Set all the elements to Affine(n,i,x)
	 *
	 * \note Emptiness is "overridden".
	 */
	void init(const Interval& x);
	void init(const IntervalVector& x);

	/**
	 * \brief Create [empty; ...; empty]
	 *
	 * Create an empty AffineMainVector of dimension \a n
	 * (all the components being empty Intervals)
	 *
	 * \pre n>0
	 */
	static AffineVarMainVector empty(int n);
	/**
	 * \brief Resize this AffineMainVector.
	 *
	 * If the size is increased, the existing components are not
	 * modified and the new ones are set to (-inf,+inf), even if
	 * (*this) is the empty Interval (however, in this case, the status of
	 * (*this) remains "empty").
	 */
	void resize(int n2);

	/**
	 * \brief Return the ith AffineVar
	 *
	 * A return a const reference to the
	 * i^th component (i starts from 0)
	 */
	const AffineVarMain<T>& operator[](int i) const;

	/**
	 * \brief Return the ith Interval
	 *
	 * A return a non-const reference to the
	 * i^th component (i starts from 0)
	 */
	AffineVarMain<T>& operator[](int i);

	/**
	 * \brief Return a subvector.
	 *
	 * \pre (*this) must not be empty
	 * \return [ (*this)[start_index]; ...; (*this)[end_index] ].
	 */
	AffineVarMainVector<T> subvector(int start_index, int end_index) const;

	/**
	 * \brief Put a subvector into *this at a given position.
	 *
	 * \param start_index - the position where the subvector
	 * \param subvec - the subvector
	 *
	 * \pre (*this) must not be empty
	 */
	void put(int start_index, const IntervalVector& subvec);

	/**
	 * \brief Assign this IntervalVector to x.
	 *
	 * \pre Dimensions of this and x must match.
	 * \note Emptiness is overridden.
	 */
	AffineVarMainVector<T>& operator=(const IntervalVector& x);

	/**
	 * \brief Return the IntervalVector compose by the interval of each Affine form
	 * \pre (*this) must be nonempty
	 */
	IntervalVector itv() const;

	/**
	 * \brief The dimension (number of components)
	 */
	int size() const;

	/**
	 * \brief Return the lower bound vector
	 * \pre (*this) must be nonempty
	 */
	Vector lb() const;

	/**
	 * \brief Return the upper bound vector
	 * \pre (*this) must be nonempty
	 */
	Vector ub() const;

	/**
	 * \brief Return the midpoint
	 * \pre (*this) must be nonempty
	 */
	Vector mid() const;

	/**
	 * \brief Return the mignitude vector.
	 * \pre (*this) must be nonempty
	 */
	Vector mig() const;

	/**
	 * \brief Return the magnitude vector.
	 * \pre (*this) must be nonempty
	 */
	Vector mag() const;

	/**
	 * \brief Vector of radii.
	 */
	Vector rad() const;

	/**
	 * \brief Return the vector of diameters.
	 */
	Vector diam() const;

	/**
	 * \brief Return true iff this AffineMainVector is empty
	 */
	bool is_empty() const;

	/**
	 * \brief true iff this interval vector contains an infinite bound.
	 *
	 * \note An empty interval vector is always bounded.
	 */
	bool is_unbounded() const;


	/**
	 * \brief Return true if the bounds of this AffineMainVector match that of \a x.
	 */
	bool operator==(const AffineMainVector<T>& x) const;
	bool operator==(const AffineVarMainVector<T>& x) const;
	bool operator==(const IntervalVector& x) const;

	/**
	 * \brief Return true if one bounds of one component of *this differs from \a x.
	 */
	bool operator!=(const IntervalVector& x) const;
	bool operator!=(const AffineMainVector<T>& x) const;
	bool operator!=(const AffineVarMainVector<T>& x) const;

private:
	void put(int start_index, const AffineMainVector<T>& subvec) {ibex_error(" AffineVarMainVector : operator put non valid");};
	AffineVarMainVector& operator=(const AffineMainVector<T>& x) {ibex_error(" AffineVarMainVector : operator= non valid");};
	AffineVarMainVector& operator+=(const Vector& x2) {ibex_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator+=(const IntervalVector& x2) {ibex_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator+=(const AffineMainVector<T>& x2) {ibex_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator+=(const AffineVarMainVector<T>& x2) {ibex_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator-=(const Vector& x2) {ibex_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator-=(const IntervalVector& x2) {ibex_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator-=(const AffineMainVector<T>& x2) {ibex_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator-=(const AffineVarMainVector<T>& x2) {ibex_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator*=(double d) {ibex_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator*=(const Interval& x1) {ibex_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator*=(const AffineMain<T>& x1) {ibex_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator*=(const AffineVarMain<T>& x1) {ibex_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator/=(double d) {ibex_error(" AffineVarMainVector : operator/= non valid");};
	AffineVarMainVector& operator/=(const Interval& x1) {ibex_error(" AffineVarMainVector : operator/= non valid");};
	AffineVarMainVector& operator/=(const AffineMain<T>& x1) {ibex_error(" AffineVarMainVector : operator/= non valid");};
	AffineVarMainVector& operator/=(const AffineVarMain<T>& x1) {ibex_error(" AffineVarMainVector : operator/= non valid");};

};




/**
 * \brief Return the intersection of x and y.
 */
template<class T>
IntervalVector operator&(const AffineVarMainVector<T>& x, const AffineVarMainVector<T>& y);
template<class T>
IntervalVector operator&(const IntervalVector& x, const AffineVarMainVector<T>& y);
template<class T>
IntervalVector operator&(const AffineVarMainVector<T>& x, const IntervalVector& y);
template<class T>
IntervalVector operator&(const AffineMainVector<T>& x, const AffineVarMainVector<T>& y);
template<class T>
IntervalVector operator&(const AffineVarMainVector<T>& x, const AffineMainVector<T>& y);

/**
 * \brief Return the hull of x & y.
 */
template<class T>
IntervalVector operator|(const AffineVarMainVector<T>& x, const AffineVarMainVector<T>& y);
template<class T>
IntervalVector operator|(const IntervalVector& x, const AffineVarMainVector<T>& y);
template<class T>
IntervalVector operator|(const AffineVarMainVector<T>& x, const IntervalVector& y);
template<class T>
IntervalVector operator|(const AffineMainVector<T>& x, const AffineVarMainVector<T>& y);
template<class T>
IntervalVector operator|(const AffineVarMainVector<T>& x, const AffineMainVector<T>& y);

/**
 * \brief -x.
 */
template<class T>
AffineMainVector<T> operator-(const AffineVarMainVector<T>& x);

/**
 * \brief x1+x2.
 */
template<class T>
AffineMainVector<T> operator+(const Vector& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainVector<T> operator+(const AffineVarMainVector<T>& x1, const Vector& x2);
template<class T>
AffineMainVector<T> operator+(const AffineVarMainVector<T>& x1, const IntervalVector& x2);
template<class T>
AffineMainVector<T> operator+(const IntervalVector& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainVector<T> operator+(const AffineVarMainVector<T>& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainVector<T> operator+(const AffineMainVector<T>& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainVector<T> operator+(const AffineVarMainVector<T>& x1, const AffineMainVector<T>& x2);

/**
 * \brief x1-x2.
 */
template<class T>
AffineMainVector<T> operator-(const Vector& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainVector<T> operator-(const AffineVarMainVector<T>& x1, const Vector& x2);
template<class T>
AffineMainVector<T> operator-(const AffineVarMainVector<T>& x1, const IntervalVector& x2);
template<class T>
AffineMainVector<T> operator-(const IntervalVector& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainVector<T> operator-(const AffineVarMainVector<T>& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainVector<T> operator-(const AffineMainVector<T>& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainVector<T> operator-(const AffineVarMainVector<T>& x1, const AffineMainVector<T>& x2);

/**
 * \brief x1*x2.
 */
template<class T>
AffineMain<T> operator*(const Vector& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMain<T> operator*(const AffineVarMainVector<T>& x1, const Vector& x2);
template<class T>
AffineMain<T> operator*(const AffineVarMainVector<T>& x1, const IntervalVector& x2);
template<class T>
AffineMain<T> operator*(const IntervalVector& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMain<T> operator*(const AffineVarMainVector<T>& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMain<T> operator*(const AffineMainVector<T>& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMain<T> operator*(const AffineVarMainVector<T>& x1, const AffineMainVector<T>& x2);

/**
 * \brief d*x
 */
template<class T>
AffineMainVector<T> operator*(double d, const AffineVarMainVector<T>& x);
template<class T>
AffineMainVector<T> operator*(const AffineMain<T>& d, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainVector<T> operator*(const Interval& d, const AffineVarMainVector<T>& x2);


// Matrix operator

/**
 * \brief $[m]*[x]$.
 */
template<class T>
AffineMainVector<T> operator*(const Matrix& m, const AffineVarMainVector<T>& x);
template<class T>
AffineMainVector<T> operator*(const AffineVarMainVector<T>& x, const Matrix& m);

/**
 * \brief $[m]*[x]$.
 */
template<class T>
AffineMainVector<T> operator*(const AffineMainMatrix<T>& m, const AffineVarMainVector<T>& x);
template<class T>
AffineMainVector<T> operator*(const IntervalMatrix& m, const AffineVarMainVector<T>& x);

/**
 * \brief $[x]*[m]$.
 */
template<class T>
AffineMainVector<T> operator*(const AffineVarMainVector<T>& x, const AffineMainMatrix<T>& m);
template<class T>
AffineMainVector<T> operator*(const AffineVarMainVector<T>& x, const IntervalMatrix& m);



/**
 * \brief Outer product (multiplication of a column vector by a row vector).
 */
template<class T>
AffineMainMatrix<T> outer_product(const AffineVarMainVector<T>& x1, const AffineMainVector<T>& x2);
template<class T>
AffineMainMatrix<T> outer_product(const AffineVarMainVector<T>& x1, const AffineVarMainVector<T>& x2);
template<class T>
AffineMainMatrix<T> outer_product(const AffineVarMainVector<T>& x1, const Vector& x2);
template<class T>
AffineMainMatrix<T> outer_product(const Vector& x1, const AffineVarMainVector<T>& x2);


/**
 * \brief |x|.
 */
template<class T>
AffineMainVector<T> abs(const AffineVarMainVector<T>& x);


/**
 * \brief Display the AffineVarMainVector<T> \a x
 */
template<class T>
std::ostream& operator<<(std::ostream& os, const AffineVarMainVector<T>& x);


template<class T>
AffineMainVector<T> operator-(const AffineVarMainVector<T>& x);
/**
 * \brief Cartesian product of x and y.
 *
 */
template<class T>
AffineVarMainVector<T> cart_prod(const AffineVarMainVector<T>& x, const AffineVarMainVector<T>& y);



// the following functions are
// introduced to allow genericity
template<class T> inline bool ___is_empty(const AffineVarMainVector<T>& v) { return v.is_empty(); }
template<class T> inline void ___set_empty(AffineVarMainVector<T>& v)      { v.set_empty(); }


} // end namespace ibex

#include "ibex_LinearArith.h"

namespace ibex {

template<class T>
AffineVarMainVector<T>& AffineVarMainVector<T>::operator=(const IntervalVector& x)  { resize(x.size()); // see issue #10
                                                                                    return _assignV(*this,x); }
template<class T>
AffineVarMainVector<T>  AffineVarMainVector<T>::subvector(int start_index, int end_index) const   { //return _subvector(*this,start_index,end_index);
	//assert(!___is_empty(v));
	assert(end_index>=0 && start_index>=0);
	assert(end_index<size() && start_index<=end_index);
	//throw InvalidIntervalVectorOp("Invalid indices for IntervalVector::subvector");

	AffineVarMainVector<T> v2(end_index-start_index+1);
	int j=0;
	for (int i=start_index; i<=end_index; i++) {
		v2[j++]=(*this)[i];
	}
	return v2;
}
template<class T>
void           AffineVarMainVector<T>::put(int start_index, const IntervalVector& subvec) {
	assert(!is_empty());
	assert(!subvec.is_empty());
	int end_index=start_index+subvec.size()-1;
	assert(start_index>=0 && end_index<size());

	int j=0;
	for (int i=start_index; i<=end_index; i++) {
		(*this)[i]=subvec[j++];
	}
}
template<class T>
bool           AffineVarMainVector<T>::operator==(const AffineMainVector<T>& x) const          { return _equalsV(*this,x); }
template<class T>
bool           AffineVarMainVector<T>::operator==(const AffineVarMainVector<T>& x) const          { return _equalsV(*this,x); }
template<class T>
bool           AffineVarMainVector<T>::operator==(const IntervalVector& x) const         { return _equalsV(*this,x); }
template<class T>
Vector         AffineVarMainVector<T>::lb() const                                        { return _lb(*this); }
template<class T>
Vector         AffineVarMainVector<T>::ub() const                                        { return _ub(*this); }
template<class T>
Vector         AffineVarMainVector<T>::mid() const                                       { return _mid(*this); }
template<class T>
Vector         AffineVarMainVector<T>::mig() const                                       { return _mig(*this); }
template<class T>
Vector         AffineVarMainVector<T>::mag() const                                       { return _mag(*this); }
//bool           AffineMainVector<T>::is_flat() const                                   { return _is_flat(*this); }
//bool           AffineMainVector<T>::contains(const Vector& x) const                   { return _contains(*this,x); }
template<class T>
bool           AffineVarMainVector<T>::is_unbounded() const                              { return _is_unbounded(*this); }
//bool           AffineMainVector<T>::is_subset(const AffineMainVector<T>& x) const           { return _is_subset(*this,x); }
//bool           AffineMainVector<T>::is_subset(const IntervalVector& x) const          { return _is_subset(*this,x); }
//bool           AffineMainVector<T>::is_strict_subset(const AffineMainVector<T>& x) const    { return _is_strict_subset(*this,x); }
//bool           AffineMainVector<T>::is_strict_subset(const IntervalVector& x) const   { return _is_strict_subset(*this,x); }
//bool           AffineMainVector<T>::is_zero() const                                   { return _is_zero(*this); }
//bool           AffineMainVector<T>::is_bisectable() const                             { return _is_bisectable(*this); }
template<class T>
Vector         AffineVarMainVector<T>::rad() const                                       { return _rad(*this); }
template<class T>
Vector         AffineVarMainVector<T>::diam() const                                      { return _diam(*this); }
//int            AffineMainVector<T>::extr_diam_index(bool min) const                   { return _extr_diam_index(*this,min); }
template<class T>
std::ostream& operator<<(std::ostream& os, const AffineVarMainVector<T>& x)              { return _displayV(os,x); }
//double         AffineMainVector<T>::volume() const                                    { return _volume(*this); }
//double         AffineMainVector<T>::perimeter() const                                 { return _perimeter(*this); }
//double         AffineMainVector<T>::rel_distance(const AffineMainVector<T>& x) const        { return _rel_distance(*this,x); }
//double         AffineMainVector<T>::rel_distance(const IntervalVector& x) const       { return _rel_distance(*this,x); }
//Vector         AffineMainVector<T>::random() const                                    { return _random<AffineMainVector<T>,Affine>(*this); }

//int          AffineMainVector<T>::diff(const AffineMainVector<T>& y, IntervalVector*& result) const   { return _diff(itv(), y.itv(), result); }
//int          AffineMainVector<T>::diff(const IntervalVector& y, IntervalVector*& result) const  { return _diff(itv(), y, result); }
//int          AffineMainVector<T>::complementary(IntervalVector*& result) const                  { return _complementary(itv(), result); }
//std::pair<IntervalVector,IntervalVector> AffineMainVector<T>::bisect(int i, double ratio) const { return _bisect(*this, i, ratio); }


template<class T>
void AffineVarMainVector<T>::init(const Interval& x) {
	for (int i = 0; i < size(); i++) {
		(*this)[i] =  x;
	}
}

template<class T>
void AffineVarMainVector<T>::init(const IntervalVector& x) {
	for (int i = 0; i < size(); i++) {
		(*this)[i] = x[i];
	}
}


template<class T>
inline AffineVarMainVector<T> AffineVarMainVector<T>::empty(int n) {
	return AffineVarMainVector<T>(n, Interval::empty_set());
}


template<class T>
inline void AffineVarMainVector<T>::set_empty() {
	(*this)[0] = Interval::empty_set();
}

template<class T>
inline void AffineVarMainVector<T>::clear() {
	init(Interval::zero());
}

template<class T>
inline int AffineVarMainVector<T>::size() const {
//	return _vec.size();
	return _n;
}

template<class T>
inline bool AffineVarMainVector<T>::is_empty() const {
	return (*this)[0].is_empty();
}


template<class T>
IntervalVector AffineVarMainVector<T>::itv() const {
	IntervalVector tmp(_n);
	for (int i = 0; i < _n; i++) {
		tmp[i] = (*this)[i].itv();
	}
	return tmp;
}

template<class T>
inline bool AffineVarMainVector<T>::operator!=(const IntervalVector& x) const {
	return !(*this==x);
}
template<class T>
inline bool AffineVarMainVector<T>::operator!=(const AffineMainVector<T>& x) const {
	return !(*this==x);
}
template<class T>
inline bool AffineVarMainVector<T>::operator!=(const AffineVarMainVector<T>& x) const {
	return !(*this==x);
}


template<class T>
inline AffineVarMainVector<T> cart_prod(const AffineVarMainVector<T>& x, const AffineVarMainVector<T>& y) {
	AffineVarMainVector<T> z(x.size()+y.size());
	z.put(0,x);
	z.put(x.size(),y);
	return z;
}

//===============================================================================================

//template<class T>
//inline AffineVarMainVector<T>::~AffineVarMainVector<T>() {
//	for (int i=0; i<_vec.size(); i++) {
//		delete _vec[i];
//	}
//}


template<class T>
const AffineVarMain<T>& AffineVarMainVector<T>::operator[](int i) const {
	assert(i>=0 && i<_n);
	return (_vec)[i];
}

template<class T>
AffineVarMain<T>& AffineVarMainVector<T>::operator[](int i) {
	assert(i>=0 && i<_n);
	return (_vec)[i];
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n) {
	assert(n>=1);
	_n   = n;
	_vec = new AffineVarMain<T>[n];
	assert(n>=1);
	for (int i = 0; i < n; i++){
		(*this)[i] = AffineVarMain<T>(n, i, Interval::all_reals());
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n, const Interval& x) {
	assert(n>=1);
	_n   = n;
	_vec  = new AffineVarMain<T>[n];
	for (int i = 0; i < n; i++) {
		(*this)[i] = AffineVarMain<T>(n, i, x);
	}
}


template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const AffineVarMainVector<T>& x) {
	_n   = x.size();
	_vec  = new AffineVarMain<T>[x.size()];
	for (int i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(),i,(x[i]).itv());
	}

}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n, double bounds[][2])  {
	assert(n>=1);
	_n   = n;
	_vec  = new AffineVarMain<T>[n];
	if (bounds == 0){ // probably, the user called AffineVarMainVector<T>(n,0) and 0 is interpreted as NULL!
		for (int i = 0; i < n; i++){
			(*this)[i] = AffineVarMain<T>(n, i, Interval(0));
		}
	}
	else {
		for (int i = 0; i < n; i++){
			(*this)[i] = AffineVarMain<T>(n, i, Interval(bounds[i][0], bounds[i][1]));
		}

	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const IntervalVector& x) {
	_n   = x.size();
	_vec  = new AffineVarMain<T>[x.size()];
	for (int i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(), i, x[i]);
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const Vector& x) {
	_n   = x.size();
	_vec = new AffineVarMain<T>[x.size()];
	for (int i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(), i, Interval(x[i]));
	}
}

template<class T>
void AffineVarMainVector<T>::resize(int n2) {
	assert(n2>=1);
	assert((_vec==NULL && _n==0) || (_n!=0 && _vec!=NULL));

	if (n2==size()) return;

	AffineVarMain<T>* newVec=new AffineVarMain<T>[n2];
	int i=0;
	for (; i<size() && i<n2; i++){
		newVec[i]=AffineVarMain<T>(n2, i, ((*this)[i]).itv());
	}
	for (; i<n2; i++) {
		newVec[i]=AffineVarMain<T>(n2, i, Interval::all_reals());
	}

	if (_vec!=NULL) {// vec==NULL happens when default constructor is used (n==0)
		delete[] _vec;
	}

	_n   = n2;
	_vec = newVec;
}


template<class T>
inline AffineMainVector<T> operator-(const AffineVarMainVector<T>& x) {
	const int n=x.size();

	AffineMainVector<T> y(n);

	if (x.is_empty()) {y.set_empty(); return y; }

	for (int i=0; i<n; i++) {
		y[i]= (-x[i]);
	}
	return y;
}


template<class T>
inline AffineMainVector<T> abs( const AffineVarMainVector<T>& x) {
	return _abs(x);
}


template<class T>
inline AffineMainVector<T> operator+(const Vector& x1, const AffineVarMainVector<T>& x2) {
	return AffineMainVector<T>(x2)+=x1;
}

template<class T>
inline AffineMainVector<T> operator+(const AffineVarMainVector<T>& x1, const Vector& x2) {
	return AffineMainVector<T>(x1)+=x2;
}

template<class T>
inline AffineMainVector<T> operator+(const IntervalVector& x1, const AffineVarMainVector<T>& x2) {
	return AffineMainVector<T>(x2)+=x1;
}

template<class T>
inline AffineMainVector<T> operator+(const AffineVarMainVector<T>& x1, const IntervalVector& x2) {
	return AffineMainVector<T>(x1)+=x2;
}

template<class T>
AffineMainVector<T> operator+(const AffineMainVector<T>& x1, const AffineVarMainVector<T>& x2) {
	return AffineMainVector<T>(x1)+=x2;
}

template<class T>
AffineMainVector<T> operator+(const AffineVarMainVector<T>& x1, const AffineMainVector<T>& x2) {
	return AffineMainVector<T>(x1)+=x2;
}

template<class T>
AffineMainVector<T> operator+(const AffineVarMainVector<T>& x1, const AffineVarMainVector<T>& x2) {
	return AffineMainVector<T>(x1)+=x2;
}


template<class T>
inline AffineMainVector<T> operator-(const Vector& x1, const AffineVarMainVector<T>& x2) {
	AffineMainVector<T> res(x2.size());
	res = (-x2);
	return res += x1;
}

template<class T>
inline AffineMainVector<T> operator-(const AffineVarMainVector<T>& x1, const Vector& x2) {
	return AffineMainVector<T>(x1)-=x2;
}

template<class T>
inline AffineMainVector<T> operator-(const AffineVarMainVector<T>& x1, const IntervalVector& x2) {
	return AffineMainVector<T>(x1)-=x2;
}

template<class T>
inline AffineMainVector<T> operator-(const IntervalVector& x1, const AffineVarMainVector<T>& x2) {
	AffineMainVector<T> res(x2.size());
	res = (-x2);
	return res += x1;
}

template<class T>
inline AffineMainVector<T> operator-(const AffineVarMainVector<T>& x1, const AffineMainVector<T>& x2) {
	return AffineMainVector<T>(x1) += (-x2);
}

template<class T>
inline AffineMainVector<T> operator-(const AffineMainVector<T>& x1, const AffineVarMainVector<T>& x2) {
	return AffineMainVector<T>(x1) += (-x2);
}

template<class T>
inline AffineMainVector<T> operator-(const AffineVarMainVector<T>& x1, const AffineVarMainVector<T>& x2) {
	return AffineMainVector<T>(x1) += (-x2);
}


//Prod Scalar


template<class T>
AffineMain<T> operator*(const Vector& v1, const AffineVarMainVector<T>& v2) {
	assert(v1.size()==v2.size());

	int n=v1.size();
	AffineMain<T> y(0);

	if (v2.is_empty()) {
		y.set_empty();
		return y;
	}

	for (int i=0; i<n; i++) {
		y+=v1[i]*v2[i];
	}
	return y;
}

template<class T>
AffineMain<T> operator*(const AffineVarMainVector<T>& v1, const Vector& v2) {
	assert(v1.size()==v2.size());

	int n=v1.size();
	AffineMain<T> y(0);

	if (v1.is_empty()) {
		y.set_empty();
		return y;
	}

	for (int i=0; i<n; i++) {
		y+=v1[i]*v2[i];
	}
	return y;
}


template<class T>
AffineMain<T> operator*(const IntervalVector& x1, const AffineVarMainVector<T>& x2){
	return x2*x1;
}

template<class T>
AffineMain<T> operator*(const AffineVarMainVector<T>& v1, const IntervalVector& v2) {
	assert(v1.size()==v2.size());

	int n=v1.size();
	AffineMain<T> y(0);

	if (v1.is_empty() || v2.is_empty()) {
		y.set_empty();
		return y;
	}

	for (int i=0; i<n; i++) {
		y+=v1[i] * v2[i];
	}
	return y;
}

template<class T>
AffineMain<T> operator*(const AffineMainVector<T>& v1, const AffineVarMainVector<T>& v2) {
	assert(v1.size()==v2.size());
	assert(v1.size()==v2.size());
	return v2*v1;
}

template<class T>
AffineMain<T> operator*(const AffineVarMainVector<T>& v1, const AffineMainVector<T>& v2) {
	assert(v1.size()==v2.size());
	assert(v1.size()==v2.size());

	int n=v1.size();
	AffineMain<T> y;

	if (v1.is_empty() || v2.is_empty()) {
		y.set_empty();
		return y;
	}
	y = 0;
	for (int i=0; i<n; i++) {
		y+=v1[i] * v2[i];
	}
	return y;
}

template<class T>
AffineMain<T> operator*(const AffineVarMainVector<T>& v1, const AffineVarMainVector<T>& v2) {
	assert(v1.size()==v2.size());
	assert(v1.size()==v2.size());

	int n=v1.size();
	AffineMain<T> y;

	if (v1.is_empty() || v2.is_empty()) {
		y.set_empty();
		return y;
	}
	y = 0;
	for (int i=0; i<n; i++) {
		y+=v1[i] * v2[i];
	}
	return y;
}


// Scalar * Vector

template<class T>
inline AffineMainVector<T> operator*(double d, const AffineVarMainVector<T>& x) {
	return AffineMainVector<T>(x)*=d;
}

template<class T>
inline AffineMainVector<T> operator*(const AffineMain<T>& x1, const AffineVarMainVector<T>& x2) {
	return AffineMainVector<T>(x2)*=x1;
}

template<class T>
inline AffineMainVector<T> operator*(const Interval& x1, const AffineVarMainVector<T>& x2) {
	return AffineMainVector<T>(x2)*=x1;
}





// Matrix

template<class T>
inline AffineMainVector<T> operator*(const IntervalMatrix& m, const AffineVarMainVector<T>& x) {
	return mulMV<IntervalMatrix,AffineVarMainVector<T>,AffineMainVector<T> >(m,x);
}

template<class T>
inline AffineMainVector<T> operator*(const AffineVarMainVector<T>& x, const Matrix& m) {
	return mulVM<AffineVarMainVector<T>,Matrix,AffineMainVector<T> >(x,m);
}


template<class T>
inline AffineMainVector<T> operator*(const AffineVarMainVector<T>& x, const IntervalMatrix& m) {
	return mulVM<AffineVarMainVector<T>,IntervalMatrix,AffineMainVector<T> >(x,m);
}


template<class T>
inline AffineMainVector<T> operator*(const AffineVarMainVector<T>& x, const AffineMainMatrix<T>& m) {
	return mulVM<AffineVarMainVector<T>,AffineMainMatrix<T>,AffineMainVector<T> >(x,m);
}

}


#endif /* IBEX_AffineVar_H_ */


/** \brief atan2(AF[y],AF[x]). */
//Affine2 atan2(const Affine2& y, const Affine2& x);
/** \brief atan2([y],AF[x]). */
//Affine2 atan2(const Interval& y, const Affine2& x);
/** \brief atan2(AF[y],[x]). */
//Affine2 atan2(const Affine2& y, const Interval& x);
/** \brief cosh(AF[x]). */
