/* ============================================================================
 * I B E X - Affine definition
 * ============================================================================
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Jordan Ninin
 * Bug fixes   :
 * Created     : Nov 12, 2012
 * ---------------------------------------------------------------------------- */

#ifndef IBEX_AFFINE_H_
#define IBEX_AFFINE_H_

#include "ibex_Interval.h"
#include <math.h>
#include <ostream>
#include <cassert>
#include "ibex_Exception.h"
#include "ibex_Array.h"
#include "ibex_TemplateDomain.h"
#include "ibex_Domain.h"


#include "ibex_Affine2_fAF2.h"
#include "ibex_Affine3_fAFFullI.h"


#ifdef _IBEX_WITH_AFFINE_EXTENDED_

#include "ibex_Affine2_fAF1.h"
#include "ibex_Affine2_iAF.h"
#include "ibex_Affine2_sAF.h"
#include "ibex_Affine2_No.h"

#endif

namespace ibex {

template<class T>  class AffineVariableMain;
template<class T>  class AffineVarMain;
template<class T>  class AffineVarMainVector;
template<class T>  class AffineMain;
template<class T> class AffineMainVector;
template<class T> class AffineMainMatrix;
template<class T> class AffineEval;


//template<class T> class TestAffineBase;
//template<class T> class TestAffineArih;

/**
 * \ingroup arithmetic
 *
 * \brief Affine Arithmetic AF2
 *
 */


//typedef AF_fAF1  AF_Default;
typedef AF_fAF2  AF_Default;
//typedef AF_fAF2_fma  AF_Default;
//typedef AF_iAF  AF_Default;
//typedef AF_sAF  AF_Default;
//typedef AF_No  AF_Default;
typedef AF_fAFFullI AF_Other;


typedef AffineMain<AF_Default> Affine2;
typedef AffineMain<AF_Other>  Affine3;

typedef AffineVarMainVector<AF_Default> Affine2Vars;
typedef AffineVarMainVector<AF_Other>  Affine3Vars;


//template<class T>
//Array< TemplateDomain< AffineMain<T> > > convert_to_affinedomain(const Array<Domain>& d);
//template<class T>
//Array< TemplateDomain< AffineMain<T> > > convert_to_affinedomain(const Array<const Domain>& d);



template<class T=AF_Default>
class AffineVarMain : public AffineMain<T> {


public:
	/** \brief Set *this to itv.
	 */
	AffineVarMain& operator=(const Interval& itv);

private:
	friend class AffineVarMainVector<T>;
//	friend class AffineMainMatrix<T>;


    const int var;

    AffineVarMain(const AffineVarMain<T>& x);
	/** \brief Create an affine form with n variables and  initialized val[0] with d. */
//	explicit AffineVarMain(double d);


	/** \brief Create an affine form with n variables and  initialized val[0] with  itv. */
//	explicit AffineVarMain(const Interval& itv);

	/** \brief Create an affine form with at most \size variables and  initialized the \var^th variable with  itv. */
	AffineVarMain(int size, int var, const Interval& itv);

	AffineVarMain& operator=(const AffineMain<T>& x) {ibex_error(" AffineVarMain : operator= non valid");};
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
		var		(var1) {
}



template<class T>
AffineVarMain<T>::AffineVarMain(const AffineVarMain& x) :
		AffineMain<T>(x.size(), x.var, x.itv()),
		var		(x.var) {
}

//===============================================================================================





template<class T=AF_Default>
class AffineVarMainVector : public AffineMainVector<T> {


private:

	AffineVarMainVector() : //_n(0),
		AffineMainVector<T>(0) { }

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
//	virtual ~AffineVarMainVector();

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
	void put(int start_index, const AffineMainVector<T>& subvec) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator=(const AffineMainVector<T>& x) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator+=(const Vector& x2) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator+=(const IntervalVector& x2) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator+=(const AffineMainVector<T>& x2) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator-=(const Vector& x2) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator-=(const IntervalVector& x2) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator-=(const AffineMainVector<T>& x2) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator*=(double d) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator*=(const Interval& x1) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator*=(const AffineMain<T>& x1) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator/=(double d) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator/=(const Interval& x1) {ibex_error(" AffineVarMainVector : operator non valid");};
	AffineVarMainVector& operator/=(const AffineMain<T>& x1) {ibex_error(" AffineVarMainVector : operator non valid");};

};




//===============================================================================================

//template<class T>
//inline AffineVarMainVector<T>::~AffineVarMainVector<T>() {
//	for (int i=0; i<_vec.size(); i++) {
//		delete _vec[i];
//	}
//}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n)  {
	this->_vec.resize(n);
	assert(n>=1);
	for (int i = 0; i < n; i++){
		this->_vec.set_ref(i,*new AffineVarMain<T>(n, i, Interval::all_reals()));
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n, const Interval& x) {
	assert(n>=1);
	this->_vec.resize(n);
	for (int i = 0; i < n; i++) {
		this->_vec.set_ref(i,*new AffineVarMain<T>(n, i, x));
	}

}


template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const AffineVarMainVector<T>& x) {
	this->_vec.resize(x.size());
	for (int i = 0; i < x.size(); i++){
		this->_vec.set_ref(i,*new AffineVarMain<T>(x.size(),i,(x[i]).itv()));
	}

}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n, double bounds[][2])  {
	this->_vec.resize(n);
	if (bounds == 0){ // probably, the user called AffineVarMainVector<T>(n,0) and 0 is interpreted as NULL!
		for (int i = 0; i < n; i++){
			this->_vec.set_ref(i,*new AffineVarMain<T>(n, i, Interval(0)));
		}
	}
	else {
		for (int i = 0; i < n; i++){
			this->_vec.set_ref(i,*new AffineVarMain<T>(n, i, Interval(bounds[i][0], bounds[i][1])));
		}

	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const IntervalVector& x) {
	this->_vec.resize(x.size());
	for (int i = 0; i < x.size(); i++){
		this->_vec.set_ref(i,*new AffineVarMain<T>(x.size(), i, x[i]));
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const Vector& x) {
	this->_vec.resize(x.size());
	for (int i = 0; i < x.size(); i++){
		this->_vec.set_ref(i,*new AffineVarMain<T>(x.size(), i, Interval(x[i])));
	}
}

template<class T>
void AffineVarMainVector<T>::resize(int n) {
	assert(n>=1);
	assert((this->_vec.is_empty() && this->_vec.size()==0) || (this->_vec.size()!=0 && !this->_vec.is_empty()));

	if (n==this->size()) return;
	else {
		int old_size= this->size();
		this->_vec.resize(n);
		Interval tmp;
		int i = 0;
		for (; i < old_size; i++){
			tmp = this->_vec[i].itv();
			this->_vec.set_ref(i,*new AffineVarMain<T>(n, i, tmp));
		}
		for (; i < n; i++){
			this->_vec.set_ref(i,*new AffineVarMain<T>(n, i, Interval::all_reals()));
		}
	}
}


template<class T=AF_Default>
class AffineMain {

protected:
//	friend class AffineVarMain<T>;
//	friend class AffineMainVector<T>;
//	friend class AffineMainMatrix<T>;
//	friend 	Array< TemplateDomain< AffineMain > > convert_to_affinedomain<T>(const Array<Domain>& d);
//	friend 	Array< TemplateDomain< AffineMain > > convert_to_affinedomain<T>(const Array<const Domain>& d);
	template<class A>
	friend std::ostream& operator<<(std::ostream& os, const AffineMain<A>& x);
	/** \brief tolerance for default compact procedure  */
	static const double AF_COMPAC_Tol;
	static const double AF_EM;
	static const double AF_EC;
	static const double AF_EE;

	static bool mode;

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
	int _actif;
	int _n; 		// dimension (size of _elt._val)-1  , ie number of variable

	T _elt;			// core of the affine form

	/** \brief Create an affine form with n variables and  initialized val[0] with d. */
	//explicit AffineMain(double d);

	/** \brief Create an affine form with n variables and  initialized the m^th variable with  itv. */
	AffineMain(int size, int var, const Interval& itv);

	/**
	 * \brief Change the number of affine variables
	 */
	void resize(int n);

public:

	typedef enum {
		AF_Default=0, AF_Chebyshev=1, AF_MinRange=2
	} Affine_Mode; // ...etc...

	/** \brief change the linearisation approximation of all the affine form: Chebyshev (by default), or Min-Range	 */
	static void change_mode(Affine_Mode tt=AF_Default);

	/** \brief Create an all_reals affine form (like Interval()). */
	AffineMain();

	/** \brief Create an affine form with n variables, initialized with x  */
	AffineMain(const AffineMain& x);

	/** \brief  Delete the affine form */
	virtual ~AffineMain() { };

	/** \brief Return -*this. */
	AffineMain operator-() const;

	/** \brief True iff *this and x are exactly the same intervals. */
	bool operator==(const AffineMain& x) const;

	/** \brief True iff *this and x are exactly the same intervals. */
	bool operator==(const Interval& x) const;

	/** \brief True iff *this and x are not exactly the same intervals. */
	bool operator!=(const AffineMain& x) const;

	/** \brief True iff *this and x are not exactly the same intervals. */
	bool operator!=(const Interval& x) const ;

	/** \brief Set this interval to the empty set. */
	void set_empty();

	/** \brief Set *this to x.
	 */
	virtual AffineMain& operator=(const AffineMain& x);

	/** \brief Set *this to d.
	 */
	AffineMain& operator=(double x);

	/** \brief Set *this to itv.
	 */
	virtual AffineMain& operator=(const Interval& itv);

	/* Union and Intersection of two Affine form must not be implemented
	 * That could produce to much confusion.
	 */
	/** \brief Intersection of *this and x.
	 * \param x - the interval to compute the intersection with.*/
	//	AffineMain& operator&=(const AffineMain& x);
	/** \brief Union of *this and I.
	 * \param x - the interval to compute the hull with.*/
	//	AffineMain& operator|=(const AffineMain& x);

	/**
	 * \brief Add [-rad,+rad] to *this.
	 *
	 * Return a reference to *this.
	 */
	AffineMain& inflate(double radd);

    /** \brief Lower bound.
     *
     * Return the lower bound of *this. */
    double lb() const;

    /** \brief Upper bound.
     *
     * Return the upper bound of *this. */
    double ub() const;

    /** \brief Midpoint.
     *
     * Returns the midpoint of *this.*/
    double mid() const;

    /**
     * \brief Radius.
     *
     * Return the diameter of *this.
     * By convention, 0 if *this is empty.*/
    double rad() const;

    /**
     * \brief Diameter.
     *
     * Return the diameter of *this.
     * By convention, 0 if *this is empty.*/
    double diam() const;

    /**
     * \brief Mignitude.
	 *
     * Returns the mignitude of *this:
     * <lu>
     * <li> +(lower bound)  if *this > 0
     * <li> -(upper bound) if *this < 0
     * <li> 0 otherwise.
     * </lu> */
    double mig() const;

    /**
     * \brief Magnitude.
	 *
     * Returns the magnitude of *this:
     * mag(*this)=max(|lower bound|, |upper bound|). */
    double mag() const;

	/**
	 * \brief number of variable represented
	 */
	int size() const;

	/**
	 * \brief Range of the affine form
	 */
	const Interval itv() const ;

	/**
	 * \brief return _val[i]
	 */
	double val(int i) const;

	/**
	 * \brief return _err
	 */
	double err() const;
	/**
	 * \brief return 1 if the affine form is actif and valid
	 */
	bool is_actif() const;

	/**
	 * \brief True iff *this is empty.
	 */
	bool is_empty() const;

	/**
	 * \brief True iff *this is empty.
	 */
	bool is_degenerated() const;

	/**
	 * \brief True if one bound of *this is infinite.
	 *
	 * \note An empty affine form is always bounded.
	 */
	bool is_unbounded() const;

	/**
	 * \brief reduce the number of noise variable if the value is inferior to \param tol
	 */
	void compact(double tol);
	void compact();

    /**
     * \brief True iff this interval is a subset of \a x.
     *
     * \note Always return true if *this is empty.
     */
    bool is_subset(const Interval& x) const;

    /**
     * \brief True iff this interval is a subset of \a x and not \a x itself.
     *
     * \note In particular, (-oo,oo) is not a strict subset of (-oo,oo)
     * and the empty set is not a strict subset of the empty set although
     * in both cases, the first is inside the interior of the second.
     */
    bool is_strict_subset(const Interval& x) const;

    /**
     * \brief True iff this interval is in the interior of \a x.
     *
    * \note In particular, (-oo,oo) is in the interior of (-oo,oo)
     * and the empty set is in the interior of the empty set.
     * \note Always return true if *this is empty.
     */
    bool is_interior_subset(const Interval& x) const;

    /**
     * \brief True iff this interval is in the relative interior of \a x.
     *
     * When x is degenerated, the relative interior of x is x itself
     * (it is not an open set). Otherwise, the relative interior of x
     * is the interior (in the usual meaning).
     */
    bool is_relative_interior_subset(const Interval& x) const;

    /**
     * \brief True iff this interval is in the interior of \a x and different from x.
     *
     * \note In particular, (-oo,oo) is not "strictly" in the interior of (-oo,oo)
     * and the empty set is not "strictly" in the interior of the empty set.
     */
    bool is_strict_interior_subset(const Interval& x) const;

    /**
     * \brief True iff this interval is a superset of \a x.
     *
     * \note Always return true if x is empty.
     */
    bool is_superset(const Interval& x) const;

    /**
     * \brief True iff this interval is a superset of \a x different from x.
     *
     * \see #is_strict_subset(const Interval&) const.
     */
    bool is_strict_superset(const Interval& x) const;

    /**
     * \brief True iff *this contains \a d.
     *
     * \note d can also be an "open bound", i.e., infinity.
     * So this function is not restricted to a set-membership
     * interpretation.
     */
    bool contains(const double& d) const;

    /**
     * \brief True iff the interior of *this contains \a d.
     *
     */
    bool interior_contains(const double& d) const;

    /**
     * \brief True iff *this and \a x intersect.
     */
    bool intersects(const Interval &x) const;

    /**
     * \brief True iff *this and \a x intersect and the intersection has a non-null volume.
     *
     * Equivalently, some interior points (of this or x) must belong to the intersection.
     */
    bool overlaps(const Interval &x) const;

    /**
     * \brief True iff *this and \a x do not intersect.
     *
     */
    bool is_disjoint(const Interval &x) const;


	/** \brief Add \a d to *this and return the result.  */
	virtual AffineMain& operator+=(double d);

	/** \brief Subtract \a d to *this and return the result. */
	virtual AffineMain& operator-=(double d);

	/** \brief Multiply *this by \a d and return the result. */
	virtual AffineMain& operator*=(double d);

	/** \brief Divide *this by \a d and return the result. */
	virtual AffineMain& operator/=(double d) ;

	/** \brief Add \a x to *this and return the result. */
	virtual AffineMain& operator+=(const Interval& x);

	/** \brief Subtract \a x to *this and return the result. */
	virtual AffineMain& operator-=(const Interval& x);

	/** \brief Multiply *this by \a x and return the result. */
	virtual AffineMain& operator*=(const Interval& x);

	/** \brief Divide *this by \a x and return the result.*/
	virtual AffineMain& operator/=(const Interval& x);

	/** \brief Add \a x to *this and return the result. */
	virtual AffineMain& operator+=(const AffineMain& x);

	/** \brief Subtract \a x to *this and return the result. */
	virtual AffineMain& operator-=(const AffineMain& x);

	/** \brief Multiply *this by \a x and return the result. */
	virtual AffineMain& operator*=(const AffineMain& x);

	/** \brief Divide *this by \a x and return the result. */
	virtual AffineMain& operator/=(const AffineMain& x);



	typedef AffineMain<T> SCALAR;
	typedef AffineMainVector<T> VECTOR;
	typedef AffineMainMatrix<T> MATRIX;



private:
	/** \brief Return sqr(*this) */
	virtual AffineMain& Asqr(const Interval& itv);
	virtual AffineMain&  Aneg();
	virtual AffineMain&  Ainv(const Interval& itv);
	virtual AffineMain&  Asqrt(const Interval& itv);
	virtual AffineMain&  Aexp(const Interval& itv);
	virtual AffineMain&  Alog(const Interval& itv);
	virtual AffineMain&  Apow(int n, const Interval& itv);
	virtual AffineMain&  Apow(double d, const Interval& itv);
	virtual AffineMain&  Apow(const Interval &y, const Interval& itvx);
	virtual AffineMain&  Aroot(int n, const Interval& itv);
	virtual AffineMain&  Acos(const Interval& itv);
	virtual AffineMain&  Asin(const Interval& itv);
	virtual AffineMain&  Atan(const Interval& itv);
	virtual AffineMain&  Aacos(const Interval& itv);
	virtual AffineMain&  Aasin(const Interval& itv);
	virtual AffineMain&  Aatan(const Interval& itv);
	virtual AffineMain&  Acosh(const Interval& itv);
	virtual AffineMain&  Asinh(const Interval& itv);
	virtual AffineMain&  Atanh(const Interval& itv);
	virtual AffineMain&  Aabs(const Interval& itv);
	virtual AffineMain&  Ainv_CH(const Interval& itv);
	virtual AffineMain&  Asqrt_CH(const Interval& itv);
	virtual AffineMain&  Aexp_CH(const Interval& itv);
	virtual AffineMain&  Alog_CH(const Interval& itv);
	virtual AffineMain&  Ainv_MR(const Interval& itv);
	virtual AffineMain&  Asqrt_MR(const Interval& itv);
	virtual AffineMain&  Aexp_MR(const Interval& itv);
	virtual AffineMain&  Alog_MR(const Interval& itv);



	template<class A>   friend AffineMain<A> operator/(double d, const AffineMain<A>& x);
	template<class A>	friend AffineMain<A> inv(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> sqr(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> sqrt(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> exp(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> log(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> pow(const AffineMain<A>&  x, int n);
	template<class A>	friend AffineMain<A> pow(const AffineMain<A>&  x, double d);
	template<class A>	friend AffineMain<A> pow(const AffineMain<A>& x, const Interval &y);
	template<class A>	friend AffineMain<A> pow(const AffineMain<A>& x, const AffineMain<A>& y);
	template<class A>	friend AffineMain<A> root(const AffineMain<A>&  x, int n);
	template<class A>	friend AffineMain<A> cos(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> sin(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> tan(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> acos(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> asin(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> atan(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> cosh(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> sinh(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> tanh(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> abs(const AffineMain<A>& x);


	template<class A> friend void AffineEval<A>::div_fwd   (int x1, int x2, int y);
	template<class A> friend void AffineEval<A>::abs_fwd    (int x, int y);
	template<class A> friend void AffineEval<A>::power_fwd  (int x, int y, int p);
	template<class A> friend void AffineEval<A>::sqr_fwd    (int x, int y);
	template<class A> friend void AffineEval<A>::sqrt_fwd   (int x, int y);
	template<class A> friend void AffineEval<A>::exp_fwd    (int x, int y);
	template<class A> friend void AffineEval<A>::log_fwd    (int x, int y);
	template<class A> friend void AffineEval<A>::cos_fwd    (int x, int y);
	template<class A> friend void AffineEval<A>::sin_fwd    (int x, int y);
	template<class A> friend void AffineEval<A>::tan_fwd    (int x, int y);
	template<class A> friend void AffineEval<A>::cosh_fwd   (int x, int y);
	template<class A> friend void AffineEval<A>::sinh_fwd   (int x, int y);
	template<class A> friend void AffineEval<A>::tanh_fwd   (int x, int y);
	template<class A> friend void AffineEval<A>::acos_fwd   (int x, int y);
	template<class A> friend void AffineEval<A>::asin_fwd   (int x, int y);
	template<class A> friend void AffineEval<A>::atan_fwd   (int x, int y);

};

// the following functions are
// introduced to allow genericity
template<class T> inline bool ___is_empty(const AffineMain<T>& x)       { return x.is_empty(); }
template<class T> inline void ___set_empty(AffineMain<T>& x)               { x.set_empty(); }


/** \ingroup arithmetic */
/*@{*/

/** \brief Stream out \a x. */
template<class T>
std::ostream& operator<<(std::ostream& os, const AffineMain<T>&  x);

template<class T>
inline double AffineMain<T>::lb() const { return this->itv().lb(); }

template<class T>
inline double AffineMain<T>::ub() const { return this->itv().ub(); }

template<class T>
inline double AffineMain<T>::rad() const { return this->itv().rad(); }

template<class T>
inline double AffineMain<T>::diam() const { return this->itv().diam(); }

template<class T>
inline double AffineMain<T>::mig() const { return this->itv().mig(); }

template<class T>
inline double AffineMain<T>::mag() const { return this->itv().mag(); }

template<class T>
inline bool AffineMain<T>::is_subset(const Interval& x) const { return this->itv().is_subset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_subset(const Interval& x) const { return this->itv().is_strict_subset(x); }

template<class T>
inline bool AffineMain<T>::is_interior_subset(const Interval& x) const { return this->itv().is_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_relative_interior_subset(const Interval& x) const { return this->itv().is_relative_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_interior_subset(const Interval& x) const { return this->itv().is_strict_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_superset(const Interval& x) const { return this->itv().is_superset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_superset(const Interval& x) const { return this->itv().is_strict_superset(x); }

template<class T>
inline bool AffineMain<T>::contains(const double& d) const { return this->itv().contains(d); }

template<class T>
inline bool AffineMain<T>::interior_contains(const double& d) const { return this->itv().interior_contains(d); }

template<class T>
inline bool AffineMain<T>::intersects(const Interval &x) const { return this->itv().intersects(x); }

template<class T>
inline bool AffineMain<T>::overlaps(const Interval &x) const { return this->itv().overlaps(x); }

template<class T>
inline bool AffineMain<T>::is_disjoint(const Interval &x) const { return this->itv().is_disjoint(x); }


/** \brief Return (-x) */
template<class T>
inline AffineMain<T> AffineMain<T>::operator-() const {
	AffineMain<T> res(*this);
	res.Aneg();
	return res;
}

/** \brief AF[x]_1+AF[x]_2. */
template<class T>
AffineMain<T> operator+(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/** \brief AF[x]+d. */
template<class T>
AffineMain<T> operator+(const AffineMain<T>&  x, double d);

/** \brief d+AF[x]. */
template<class T>
AffineMain<T> operator+(double d, const AffineMain<T>&  x);

/** \brief AF[x]_1+[x]_2. */
template<class T>
AffineMain<T> operator+(const AffineMain<T>&  x1, const Interval& x2);

/** \brief [x]_1+AF[x]_2. */
template<class T>
AffineMain<T> operator+(const Interval& x1, const AffineMain<T>&  x2);

/** \brief AF[x]_1-AF[x]_2. */
template<class T>
AffineMain<T> operator-(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/** \brief AF[x]-d. */
template<class T>
AffineMain<T> operator-(const AffineMain<T>&  x, double d);

/** \brief d-AF[x]. */
template<class T>
AffineMain<T> operator-(double d, const AffineMain<T>&  x);

/** \brief AF[x]_1-[x]_2. */
template<class T>
AffineMain<T> operator-(const AffineMain<T>&  x1, const Interval& x2);

/** \brief [x]_1-AF[x]_2. */
template<class T>
AffineMain<T> operator-(const Interval& x1, const AffineMain<T>&  x2);

/** \brief AF[x]_1*AF[x]_2. */
template<class T>
AffineMain<T> operator*(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/** \brief AF[x]*d. */
template<class T>
AffineMain<T> operator*(const AffineMain<T>&  x, double d);

/** \brief d*AF[x]. */
template<class T>
AffineMain<T> operator*(double d, const AffineMain<T>&  x);

/** \brief AF[x]_1*[x]_2. */
template<class T>
AffineMain<T> operator*(const AffineMain<T>&  x1, const Interval& x2);

/** \brief [x]_1*AF[x]_2. */
template<class T>
AffineMain<T> operator*(const Interval& x1, const AffineMain<T>&  x2);

/** \brief AF[x]_1/AF[x]_2. */
template<class T>
AffineMain<T> operator/(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/** \brief AF[x]/d. */
template<class T>
AffineMain<T> operator/(const AffineMain<T>&  x, double d);

/** \brief d/AF[x]. */
template<class T>
AffineMain<T> operator/(double d, const AffineMain<T>&  x);

/** \brief AF[x]_1/[x]_2. */
template<class T>
AffineMain<T> operator/(const AffineMain<T>&  x1, const Interval& x2);

/** \brief [x]_1/AF[x]_2. */
template<class T>
AffineMain<T> operator/(const Interval& x1, const AffineMain<T>&  x2);


/** \brief Hausdorff distance of AF[x]_1 and AF[x]_2. */
template<class T>
double distance(const AffineMain<T>& x1, const AffineMain<T>& x2);

/** \brief Hausdorff distance of [x]_1 and AF[x]_2. */
template<class T>
double distance(const Interval &x1, const AffineMain<T>& x2);

/** \brief Hausdorff distance of AF[x]_1 and [x]_2. */
template<class T>
double distance(const AffineMain<T>& x1, const Interval &x2);


/** \brief 1/AF[x] */
template<class T>
AffineMain<T> inv(const AffineMain<T>&  x);

/** \brief AF[x]^2 */
template<class T>
AffineMain<T> sqr(const AffineMain<T>&  x);

/** \brief sqrt{AF[x]} */
template<class T>
AffineMain<T> sqrt(const AffineMain<T>&  x);

/** \brief exp(AF[x]). */
template<class T>
AffineMain<T> exp(const AffineMain<T>&  x);

/** \brief log(AF[x]). */
template<class T>
AffineMain<T> log(const AffineMain<T>&  x);

/** \brief AF[x]^n. */
template<class T>
AffineMain<T> pow(const AffineMain<T>&  x, int n);

/** \brief AF[x]^d. */
template<class T>
AffineMain<T> pow(const AffineMain<T>&  x, double d);

/** \brief AF[x]^[y]. */
template<class T>
AffineMain<T> pow(const AffineMain<T>& x, const Interval &y);

/** \brief AF[x]^AF[y]. */
template<class T>
AffineMain<T> pow(const AffineMain<T>& x, const AffineMain<T>& y);

/** \brief nth root of AF[x]. */
template<class T>
AffineMain<T> root(const AffineMain<T>&  x, int n);

/** \brief cos(AF[x]). */
template<class T>
AffineMain<T> cos(const AffineMain<T>&  x);

/** \brief sin(AF[x]). */
template<class T>
AffineMain<T> sin(const AffineMain<T>&  x);

/** \brief tan(AF[x]). */
template<class T>
AffineMain<T> tan(const AffineMain<T>&  x);

/** \brief acos(AF[x]). */
template<class T>
AffineMain<T> acos(const AffineMain<T>&  x);

/** \brief asin(AF[x]). */
template<class T>
AffineMain<T> asin(const AffineMain<T>&  x);

/** \brief atan(AF[x]). */
template<class T>
AffineMain<T> atan(const AffineMain<T>&  x);

/** \brief cosh(AF[x]). */
template<class T>
AffineMain<T> cosh(const AffineMain<T>&  x);

/** \brief sinh(AF[x]). */
template<class T>
AffineMain<T> sinh(const AffineMain<T>&  x);

/** \brief tanh(AF[x]). */
template<class T>
AffineMain<T> tanh(const AffineMain<T>&  x);

/** \brief abs(AF[x]) = sqrt(sqr(AF[x])). */
template<class T>
AffineMain<T> abs(const AffineMain<T>& x);

template<class T>
Interval max(const AffineMain<T>&  x, const AffineMain<T>&  y);
template<class T>
Interval max(const Interval& x, const AffineMain<T>&  y);
template<class T>
Interval max(const AffineMain<T>&  x, const Interval& y);

template<class T>
Interval min(const AffineMain<T>&  x, const AffineMain<T>&  y);
template<class T>
Interval min(const Interval& x, const AffineMain<T>&  y) ;
template<class T>
Interval min(const AffineMain<T>&  x, const Interval& y);


/** \brief \f$[x]_1\cap [x]_2\f$.
 * \return Interval::EMPTY if the intersection is empty. */
template<class T>
Interval operator&(const AffineMain<T>&  x1, const AffineMain<T>&  x2);
template<class T>
Interval operator&(const Interval& x1, const AffineMain<T>&  x2);
template<class T>
Interval operator&(const AffineMain<T>&  x1, const Interval& x2);


/** \brief \f$\square([x]_1\cup [x]_2)\f$. */
template<class T>
Interval operator|(const AffineMain<T>&  x1, const AffineMain<T>&  x2);
template<class T>
Interval operator|(const Interval& x1, const AffineMain<T>&  x2);
template<class T>
Interval operator|(const AffineMain<T>&  x1, const Interval& x2);


/**
 * \brief Return the largest integer interval included in x.
 *  TODO Affine2::integer not yet implemented integer
 */
template<class T>
AffineMain<T> integer(const AffineMain<T>&  x);

/** \brief Sign of AF[x].
 *
 *  Return \f$sign(AF[x]) = hull \{ sign{x}, x\inAF[x] \}\f$.
 * \remark By convention, \f$ 0\inAF[x] \Longrightarrow sign(AF[x])=[-1,1]\f$. */
template<class T>
AffineMain<T> sign(const AffineMain<T>&  x);



/** \brief Chi of AF[a], AF[b] and AF[c].
 *
 *  Return \f$chi(AF[a],AF[b],AF[c]) = 0.5*(1-\sign(AF[a]))*AF[b] + 0.5*(\sign(AF[a])+1)*AF[c]. \f$
 * \remark  chi(AF[a],AF[b],AF[c]) =AF[b] if [a]<=0, AF[c] if AF[a]>0, hull \{AF[b], AF[c]\} else.  */
template<class T>
AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const AffineMain<T>&  c);
template<class T>
AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const AffineMain<T>&  c);
template<class T>
AffineMain<T> chi(const Interval&  a,const Interval&  b,const AffineMain<T>&  c);
template<class T>
AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const Interval&  c);
template<class T>
AffineMain<T> chi(const AffineMain<T>&  a,const Interval&  b,const AffineMain<T>&  c);
template<class T>
AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const Interval&  c);





}


/*@}*/

/*============================================ inline implementation ============================================ */



namespace ibex {


template<class T> const double AffineMain<T>::AF_COMPAC_Tol = 1.e-6;
template<class T> const double AffineMain<T>::AF_EM = __builtin_powi(2.0, -51);
template<class T> const double AffineMain<T>::AF_EC = __builtin_powi(2.0, -55);
template<class T> const double AffineMain<T>::AF_EE = 2.0;
template<class T> bool AffineMain<T>::mode=true;



template<class T>
inline void AffineMain<T>::change_mode(Affine_Mode tt) {
	switch(tt) {
	case AF_Default:
	case AF_Chebyshev: {
		mode =true;
		break;
	}
	case AF_MinRange:
		mode =false;
		break;
	}
}

template<class T>
inline void AffineMain<T>::compact(){	compact(AF_COMPAC_Tol); }

template<class T>
inline Interval operator&(const AffineMain<T>& x1, const AffineMain<T>& x2) {
	return ((x1.itv())&(x2.itv()));
}

template<class T>
inline Interval operator&(const Interval& x1, const AffineMain<T>& x2) {
	return ((x1)&(x2.itv()));
}

template<class T>
inline Interval operator&(const AffineMain<T>& x1, const Interval& x2) {
	return ((x1.itv())&(x2));
}

template<class T>
inline Interval operator|(const AffineMain<T>& x1, const AffineMain<T>& x2) {
	return ((x1.itv())|(x2.itv()));
}

template<class T>
inline Interval operator|(const Interval& x1, const AffineMain<T>& x2) {
	return ((x1)|(x2.itv()));
}

template<class T>
inline Interval operator|(const AffineMain<T>& x1, const Interval& x2) {
	return ((x1.itv())|(x2));
}

template<class T>
inline bool AffineMain<T>::operator==(const AffineMain<T>& x) const{
	return (this->itv() == x.itv());
}

template<class T>
inline bool AffineMain<T>::operator==(const Interval& x) const{
	return (this->itv() == x);
}

template<class T>
inline bool AffineMain<T>::operator!=(const AffineMain<T>& x) const{
	return (this->itv() != x.itv());
}

template<class T>
inline bool AffineMain<T>::operator!=(const Interval& x) const {
	return (this->itv() != x);
}

template<class T>
inline void AffineMain<T>::set_empty(){
	*this = Interval::empty_set();
}

//template<class T>
//inline AffineMain<T>& AffineMain<T>::inflate(double radd){
//	if (fabs(radd)>= POS_INFINITY) {
//		*this = Interval::all_reals();
//	} else {
//		saxpy(1.0, AffineMain<T>(),0.0, radd, false, false, false, true);
//	}
//	return *this;
//}

template<class T>
inline int AffineMain<T>::size() const{
	return _n;
}


template<class T>
inline bool AffineMain<T>::is_actif() const{
	return (_actif>-1);
}

template<class T>
inline bool AffineMain<T>::is_empty() const{
	return (_actif==-1);
}

template<class T>
inline bool AffineMain<T>::is_degenerated() const {
	//return (itv().diam() <	AF_EC);
	return (_actif==0);
}

template<class T>
inline bool AffineMain<T>::is_unbounded() const{
	return ((-1>_actif)&&(_actif>-5));
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(double d){
	*this += (-d);
	return *this;
}

template<class T>
inline 	AffineMain<T>& AffineMain<T>::operator/=(double d) {
	return *this *= (1.0 / Interval(d)) ;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator+=(const Interval& x){
	if ( x.is_empty()) {
		*this = Interval::empty_set();
	} else if (x.is_unbounded() || (!is_actif())) {
		*this = this->itv() + x;
	} else {
		*this += x.mid();
		this->inflate(x.rad());
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(const Interval& x){
	return *this += (-x);
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator/=(const Interval& x){
	return *this *= (1.0/x) ;
}

//template<class T>
//inline AffineMain<T>& AffineMain<T>::operator+=(const AffineMain<T>& x){
//	return saxpy(1.0, x, 0.0, 0.0, false, true, false, false);
//}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(const AffineMain<T>& x){
	return *this += (-AffineMain<T>(x));
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator/=(const AffineMain<T>& x){
	return *this *= (AffineMain<T>(x).Ainv(x.itv()));
}

template<class T>
inline AffineMain<T> operator+(const AffineMain<T>& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) += x2;
}

template<class T>
inline AffineMain<T> operator+(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) += d;
}

template<class T>
inline AffineMain<T> operator+(double d, const AffineMain<T>& x){
	return AffineMain<T>(x) += d;
}

template<class T>
inline AffineMain<T> operator+(const AffineMain<T>& x1, const Interval& x2){
	return AffineMain<T>(x1) += x2;
}

template<class T>
inline AffineMain<T> operator+(const Interval& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x2) += x1;
}

template<class T>
inline AffineMain<T> operator-(const AffineMain<T>& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) += (-x2);
}

template<class T>
inline AffineMain<T> operator-(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) -= d;
}

template<class T>
inline AffineMain<T> operator-(double d, const AffineMain<T>& x){
	AffineMain<T> res = (-x);
	return res += d;
}

template<class T>
inline AffineMain<T> operator-(const AffineMain<T>& x1, const Interval& x2) {
	return AffineMain<T>(x1) -= x2;
}

template<class T>
inline AffineMain<T> operator-(const Interval& x1, const AffineMain<T>& x2) {
	AffineMain<T> res = (- x2);
	return res += x1;
}

template<class T>
inline AffineMain<T> operator*(const AffineMain<T>& x1, const AffineMain<T>& x2) {
	return AffineMain<T>(x1) *= x2;
}

template<class T>
inline AffineMain<T> operator*(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) *= d;
}

template<class T>
inline AffineMain<T> operator*(double d, const AffineMain<T>& x){
	return AffineMain<T>(x) *= d;
}

template<class T>
inline AffineMain<T> operator*(const AffineMain<T>& x1, const Interval& x2){
	return AffineMain<T>(x1) *= x2;
}

template<class T>
inline AffineMain<T> operator*(const Interval& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x2) *= x1;
}

template<class T>
inline AffineMain<T> operator/(const AffineMain<T>& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) /= x2;
}

template<class T>
inline AffineMain<T> operator/(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) /= d;
}

template<class T>
inline AffineMain<T> operator/(double d, const AffineMain<T>& x){
	AffineMain<T> out;
	out = d;
	return out *= (AffineMain<T>(x).Ainv(x.itv()));
}

template<class T>
inline AffineMain<T> operator/(const AffineMain<T>& x1, const Interval& x2){
	return AffineMain<T>(x1) /= x2;
}

template<class T>
inline AffineMain<T> operator/(const Interval& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) *= (AffineMain<T>(x2).Ainv(x2.itv()));
}

template<class T>
inline double distance(const AffineMain<T> &x1, const AffineMain<T> &x2){
	return distance(x1.itv(), x2.itv());
}

template<class T>
inline double distance(const Interval &x1, const AffineMain<T> &x2){
	return distance(x1, x2.itv());
}

template<class T>
inline double distance(const AffineMain<T> &x1, const Interval &x2){
	return distance(x1.itv(), x2);
}

template<class T>
inline AffineMain<T> inv(const AffineMain<T>& x){
	return AffineMain<T>(x).Ainv(x.itv());
}

template<class T>
inline AffineMain<T> sqr(const AffineMain<T>& x){
	return AffineMain<T>(x).Asqr(x.itv());
}

template<class T>
inline AffineMain<T> sqrt(const AffineMain<T>& x){
	return AffineMain<T>(x).Asqrt(x.itv());
}

template<class T>
inline AffineMain<T> exp(const AffineMain<T>& x){
	return AffineMain<T>(x).Aexp(x.itv());
}

template<class T>
inline AffineMain<T> log(const AffineMain<T>& x){
	return AffineMain<T>(x).Alog(x.itv());
}

template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, int n) {
	return AffineMain<T>(x).Apow(n,x.itv());
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, double d){
	return AffineMain<T>(x).Apow(d,x.itv());
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T> &x, const Interval &y){
	// return exp(y * log(x));
	return AffineMain<T>(x).Apow(y,x.itv());
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T> &x, const AffineMain<T> &y){
	return AffineMain<T>(x).Apow(y.itv(),x.itv());
}

template<class T>
inline AffineMain<T> root(const AffineMain<T>& x, int n) {
	return AffineMain<T>(x).Aroot(n,x.itv());
}

template<class T>
inline AffineMain<T> cos(const AffineMain<T>& x){
	return AffineMain<T>(x).Acos(x.itv());
}

template<class T>
inline AffineMain<T> sin(const AffineMain<T>& x){
	return AffineMain<T>(x).Asin(x.itv());
}

template<class T>
inline AffineMain<T> tan(const AffineMain<T>& x){
	return AffineMain<T>(x).Atan(x.itv());
}

template<class T>
inline AffineMain<T> acos(const AffineMain<T>& x){
	return AffineMain<T>(x).Aacos(x.itv());
}

template<class T>
inline AffineMain<T> asin(const AffineMain<T>& x){
	return AffineMain<T>(x).Aasin(x.itv());
}

template<class T>
inline AffineMain<T> atan(const AffineMain<T>& x){
	return AffineMain<T>(x).Aatan(x.itv());
}

template<class T>
inline AffineMain<T> cosh(const AffineMain<T>& x){
	return AffineMain<T>(x).Acosh(x.itv());
}

template<class T>
inline AffineMain<T> sinh(const AffineMain<T>& x){
	return AffineMain<T>(x).Asinh(x.itv());
}

template<class T>
inline AffineMain<T> tanh(const AffineMain<T>& x){
	return AffineMain<T>(x).Atanh(x.itv());
}

template<class T>
inline AffineMain<T> abs(const AffineMain<T> &x){
	return AffineMain<T>(x).Aabs(x.itv());
}


template<class T>
inline Interval max(const AffineMain<T>& x, const AffineMain<T>& y) {
	return max(x.itv(), y.itv());
}
template<class T>
inline Interval max(const Interval& x, const AffineMain<T>& y) {
	return max(x, y.itv());
}
template<class T>
inline Interval max(const AffineMain<T>& x, const Interval& y) {
	return max(x.itv(), y);
}


template<class T>
inline Interval min(const AffineMain<T>& x, const AffineMain<T>& y) {
	return min(x.itv(), y.itv());
}
template<class T>
inline Interval min(const Interval& x, const AffineMain<T>& y) {
	return min(x, y.itv());
}
template<class T>
inline Interval min(const AffineMain<T>& x, const Interval& y) {
	return min(x.itv(), y);
}


template<class T>
inline AffineMain<T> integer(const AffineMain<T>& x){
	Interval tmp= x.itv();
	if ((tmp.lb()==(int) tmp.lb())&&(tmp.ub()==(int) tmp.ub())) {
		return AffineMain<T>(x);
	} else  {
		return AffineMain<T>(integer(tmp));
	}
}

template<class T>
inline AffineMain<T> sign(const AffineMain<T>& x) {
	return AffineMain<T>(sign(x.itv()));
}


template<class T>
inline AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const AffineMain<T>&  c) {
	return chi(a.itv(),b,c);
}
template<class T>
inline AffineMain<T> chi(const Interval&  a,const Interval&  b,const AffineMain<T>&  c) {
	return chi(a,AffineMain<T>(b),c);
}
template<class T>
inline AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const Interval&  c) {
	return chi(a,b,AffineMain<T>(c));
}
template<class T>
inline AffineMain<T> chi(const AffineMain<T>&  a,const Interval&  b,const AffineMain<T>&  c) {
	return chi(a.itv(),AffineMain<T>(b),c);
}
template<class T>
inline AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const Interval&  c){
	return chi(a.itv(),b,AffineMain<T>(c));
}
template<class T>
inline AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const AffineMain<T>&  c){
	if (a.ub()<=0) {
		return AffineMain<T>(b);
	} else if (a.lb()>0) {
		return AffineMain<T>(c);
	} else {
		AffineMain<T> out;
		out = b|c;
		return  out;
	}
}

template<class T>
inline std::ostream& operator<<(std::ostream& os, const AffineMain<T>& x) {
	{
		os << x.itv() << " : ";
		if (x.is_actif()) {
			os << x.mid();
			for (int i = 0; i < x.size(); i++) {
				os << " + " << x.val(i) << " eps_" << i;
			}
			os << " + " << x.err() << " [-1,1] ";
		} else {
			os << "AffineMain form not Activate ";
		}
		return os;

	}
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Ainv(const Interval& itv){
	if (mode) {
		this->Ainv_CH(itv);
	} else {
		this->Ainv_MR(itv);
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aexp(const Interval& itv){
	if (mode) {
		this->Aexp_CH(itv);
	} else {
		this->Aexp_MR(itv);
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Alog(const Interval& itv){
	if (mode) {
		this->Alog_CH(itv);
	} else {
		this->Alog_MR(itv);
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Asqrt(const Interval& itv){
	if (mode) {
		this->Asqrt_CH(itv);
	} else {
		this->Asqrt_MR(itv);
	}
	return *this;
}


// Chebyshev formula
// alpha = (F(sup(x)) - F(inf(x)))/diam(X)
// u = (f')^{-1}(alpha)
// d_a = f(inf(x)) -alpha*inf(X)
// d_b = f(sup(x)) -alpha*sup(x)
// d_min = min(d_a,d_b)
// d_max = f(u) - alpha*u
// beta = Interval(d_min,d_max).mid()
// zeta = Interval(d_min,d_max).rad()

template<class T>
inline AffineMain<T>& AffineMain<T>::Ainv_CH(const Interval& itv){
	Interval res_itv = 1.0/(itv);;

	// Particular case
	if ((itv.is_unbounded()) || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);
		Interval itv2 =abs(itv);
		dmm = (1.0/(itv2));
		alpha = -(dmm.diam()/itv2.diam());

		//u = 1/sqrt(-alpha);
		TEMP1 = (1.0/Interval((itv2).lb()))-alpha*Interval((itv2).lb());
		TEMP2 = (1.0/Interval((itv2).ub()))-alpha*Interval((itv2).ub());
		if (TEMP1.ub()>TEMP2.ub()) {
			band = Interval((2*sqrt(-Interval(alpha))).lb(),TEMP1.ub());
		}
		else {
			band = Interval((2*sqrt(-Interval(alpha))).lb(),TEMP2.ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		if (itv.lb()<0.0) beta = -beta;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}


	return *this;
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Asqrt_CH(const Interval& itv){
	Interval itv2 = itv & Interval::pos_reals();
	Interval res_itv = sqrt(itv2);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);

		alpha = res_itv.diam()/itv2.diam();

		//u = 1/(4*alpha^2);
		TEMP1 = res_itv.lb()-alpha*Interval(itv2.lb());
		TEMP2 = res_itv.ub()-alpha*Interval(itv2.ub());
		if (TEMP1.lb()>TEMP2.lb()) {
			band = Interval(TEMP2.lb(),(1.0/(4*Interval(alpha))).ub());
		} else {
			band = Interval(TEMP1.lb(),(1.0/(4*Interval(alpha))).ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band-beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aexp_CH(const Interval& itv){
	Interval res_itv = exp(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);
		alpha = res_itv.diam()/itv.diam();

		//u = log(alpha);
		TEMP1 = res_itv.lb()-alpha*Interval(itv.lb());
		TEMP2 = res_itv.ub()-alpha*Interval(itv.ub());
		if (TEMP1.ub()>TEMP2.ub()) {
			band = Interval((alpha*(1-log(Interval(alpha)))).lb(),TEMP1.ub());
		} else {
			band = Interval((alpha*(1-log(Interval(alpha)))).lb(),TEMP2.ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band-beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Alog_CH(const Interval& itv){
	Interval res_itv = log(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);
		alpha = res_itv.diam()/itv.diam();

		//u = 1/alpha;
		TEMP1 = res_itv.lb()-alpha*Interval(itv.lb());
		TEMP2 = res_itv.ub()-alpha*Interval(itv.ub());
		if (TEMP1.lb()>TEMP2.lb()) {
			band = Interval(TEMP2.lb(),(-log(Interval(alpha))-1).ub());
		}
		else {
			band = Interval(TEMP1.lb(),(-log(Interval(alpha))-1).ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Apow(const Interval &y, const Interval& itvx){
	// return exp(y * log(x));
	this->Alog(itvx);
	*this *= y;
	this->Aexp((y*log(itvx)));
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Acos(const Interval& itv){
	Interval res_itv = cos(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1 ;//, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);
		if (itv.diam()>=Interval::two_pi().lb()) {
			*this = Interval(-1,1);
			return *this;
		}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*x0 +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());


		fxb0 = ::cos(xb0);
		fxb1 = ::cos(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((size()*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;
		Interval u,nb_period;

		// compute the error at _itv.lb() and _itv.ub() and compute the first point such as f'(u) = alpha

		ddelta = (abs(cos(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(cos(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		u = asin(-Interval(alpha));
		nb_period = (itv+Interval::half_pi()) / Interval::pi();

		// evaluate the error at the point such as f'(u) = alpha
		int p1 = ((int) nb_period.lb())-2;
		int p2 = ((int) nb_period.ub())+2;

		int i = p1;

		while (i<=p2) { // looking for a point
			TEMP1 = (itv & (i%2==0? (u + i*Interval::pi()) : (i*Interval::pi() - u)));
			if (!(TEMP1.is_empty())) { // check if maximize the error
				t1 = (abs(cos(TEMP1)-(alpha*TEMP1+beta))).ub();
				if (t1>ddelta)  ddelta= t1;
			}
			i++;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Asin(const Interval& itv){
	Interval res_itv = sin(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1;//, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);
		if (itv.diam()>=Interval::two_pi().lb()) {
			*this = Interval(-1,1);
			return *this;
		}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*x0 +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());


		fxb0 = ::sin(xb0);
		fxb1 = ::sin(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((size()*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;
		Interval u,nb_period;

		// compute the error at _itv.lb() and _itv.ub() and compute the first point such as f'(u) = alpha

		ddelta = (abs(sin(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(sin(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		u = acos(Interval(alpha));
		nb_period = (itv) / Interval::pi();

		// evaluate the error at the point such as f'(u) = alpha
		int p1 = ((int) nb_period.lb())-2;
		int p2 = ((int) nb_period.ub())+2;

		int i = p1;

		while (i<=p2) { // looking for a point
			TEMP1 = (itv & (i%2==0? (u + i*Interval::pi()) : ((i+1)*Interval::pi() - u)));
			if (!(TEMP1.is_empty())) {
				t1 = (abs(sin(TEMP1)-(alpha*TEMP1+beta))).ub();
				if (t1>ddelta)  ddelta= t1;
			}
			i++;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Atan(const Interval& itv){
	Interval res_itv = tan(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1;//, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);
		if (itv.diam()>=Interval::two_pi().lb()) {
			*this = Interval(-1,1);
			return *this;
		}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*x0 +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());


		fxb0 = ::tan(xb0);
		fxb1 = ::tan(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((size()*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;
		Interval u,nb_period;

		// compute the error at _itv.lb() and _itv.ub() and compute the first point such as f'(u) = alpha

		ddelta = (abs(tan(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(tan(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		u = acos(1/sqrt(Interval(alpha)));
		nb_period = (itv) / Interval::pi();

		// evaluate the error at the point such as f'(u) = alpha
		int p1 = ((int) nb_period.lb())-2;
		int p2 = ((int) nb_period.ub())+2;

		int i = p1;

		while (i<=p2) { // looking for a point
			TEMP1 = (itv & ( i*Interval::pi() + u));
			if ((!(TEMP1.is_empty()))) {
				t1 = (abs(tan(TEMP1)-(alpha*TEMP1+beta))).ub();
				if (t1>ddelta)  ddelta= t1;
			}
			TEMP1 = (itv & ( i*Interval::pi() - u ));
			if ((!(TEMP1.is_empty()))) {
				t1 = (abs(tan(TEMP1)-(alpha*TEMP1+beta))).ub();
				if (t1>ddelta)  ddelta= t1;
			}
			i++;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aacos(const Interval& itv){
	Interval res_itv = acos(itv);;
	Interval itv2 = itv & Interval(-1,1);
	// Particular case
	if ( res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv2.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1;//, t2;

		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/::sqrt(2.);
		xb0 = (0.5)*(itv2.diam()*  x0  +itv2.lb()+itv2.ub());
		xb1 = (0.5)*(itv2.diam()*(-x0) +itv2.lb()+itv2.ub());

		fxb0 = ::acos(xb0);
		fxb1 = ::acos(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv2.diam());
		beta   = c0-c1*((itv2.lb()+itv2.ub())/(itv2.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha

		ddelta = (abs(acos(Interval(itv2.lb()))-(alpha*Interval(itv2.lb())+beta))).ub();
		t1     = (abs(acos(Interval(itv2.ub()))-(alpha*Interval(itv2.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// acos'(u)=-1/sqrt(1-u^2) = alpha
		// u = +-sqrt(1-1/(alpha^2))
		Interval TEMP2 = sqrt(1-1/(pow(Interval(alpha),2)));
		if (!((TEMP2 & itv2).is_empty())) {
			t1 = (abs(acos(TEMP2)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv2).is_empty())) {
			t1 = (abs(acos(-TEMP2)-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aasin(const Interval& itv){
	Interval res_itv = asin(itv);;
	Interval itv2 =itv & Interval(-1,1);
	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv2.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1;//, t2;

		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/::sqrt(2.);
		xb0 = (0.5)*(itv2.diam()*  x0  +itv2.lb()+itv2.ub());
		xb1 = (0.5)*(itv2.diam()*(-x0) +itv2.lb()+itv2.ub());

		fxb0 = ::asin(xb0);
		fxb1 = ::asin(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv2.diam());
		beta   = c0-c1*((itv2.lb()+itv2.ub())/(itv2.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha

		ddelta = (abs(asin(Interval(itv2.lb()))-(alpha*Interval(itv2.lb())+beta))).ub();
		t1     = (abs(asin(Interval(itv2.ub()))-(alpha*Interval(itv2.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// asin'(u)=1/sqrt(1-u^2) = alpha
		// u = sqrt(1-1/(alpha^2))
		Interval TEMP2 = sqrt(1/Interval(alpha)-1);
		if (!((TEMP2 & itv2).is_empty())) {
			// tanh(acosh(x)) = sqrt(sqr(x)-1) :x
			t1 = (abs(asin(TEMP2)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv2).is_empty())) {
			// sinh(acosh(-x)) = -sqrt(sqr(x)-1)
			t1 = (abs(asin(-TEMP2)-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aatan(const Interval& itv){
	Interval res_itv = atan(itv);;

	// Particular case
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1;//, t2;

		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*  x0  +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());

		fxb0 = ::atan(xb0);
		fxb1 = ::atan(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha

		ddelta = (abs(atan(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(atan(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// atan'(u)=1/(u^2+1) = alpha
		// u = +-sqrt(1/alpha -1)
		Interval TEMP2 = sqrt(1/Interval(alpha)-1);
		if (!((TEMP2 & itv).is_empty())) {
			t1 = (abs(atan(TEMP2)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv).is_empty())) {
			t1 = (abs(atan(-TEMP2)-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Acosh(const Interval& itv){
	Interval res_itv = cosh(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);
		alpha = ((cosh(Interval(itv.ub()))-cosh(Interval(itv.lb())))/itv.diam()).lb();

		//u = asinh(alpha);
		TEMP1 = cosh(Interval(itv.lb()))-alpha*Interval(itv.lb());
		TEMP2 = cosh(Interval(itv.ub()))-alpha*Interval(itv.ub());
		if (TEMP1.ub()>TEMP2.ub()) {
			// cosh(asinh(alpha)) = sqrt(sqr(alpha)+1)
			band = Interval((sqrt(pow(Interval(alpha),2)+1)-alpha*asinh(Interval(alpha))).lb(),TEMP1.ub());
		}
		else {
			band = Interval((sqrt(pow(Interval(alpha),2)+1)-alpha*asinh(Interval(alpha))).lb(),TEMP2.ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Asinh(const Interval& itv){
	Interval res_itv = sinh(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1;//, t2;
		Interval  TEMP2(0.0);
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*  x0  +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());

		fxb0 = ::sinh(xb0);
		fxb1 = ::sinh(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha
		ddelta = (abs(sinh(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(sinh(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// u = acosh(alpha)
		TEMP2 = acosh(Interval(alpha));
		if (!((TEMP2 & itv).is_empty())) {
			// sinh(acosh(x)) = sqrt(sqr(x)-1)
			t1 = (abs(sqrt(pow(Interval(alpha),2)-1)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv).is_empty())) {
			// sinh(acosh(-x)) = -sqrt(sqr(x)-1)
			t1 = (abs((-sqrt(pow(Interval(alpha),2)-1))-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Atanh(const Interval& itv){
	Interval res_itv = tanh(itv);;

	// Particular case
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1;//, t2;
		Interval  TEMP2(0.0);
		// additional particular case}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*  x0  +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());

		fxb0 = ::tanh(xb0);
		fxb1 = ::tanh(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha
		ddelta = (abs(tanh(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(tanh(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// tanh'(u)=alpha
		// cosh(u)= -2/alpha -1
		// u = +-acosh(-2/alpha -1)
		TEMP2 = acosh(-2/Interval(alpha) -1);
		if (!((TEMP2 & itv).is_empty())) {
			t1 = (abs(tanh(TEMP2)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv).is_empty())) {
			t1 = (abs(tanh(-TEMP2)-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aabs(const Interval& itv){
	Interval res_itv = abs(itv);;

	// Particular case
	if (0<=itv.lb()) return *this;
	if (itv.ub()<=0) {
		this->Aneg();
		return *this;
	}
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
	// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);

		alpha = ((abs(Interval(itv.ub()))-abs(Interval(itv.lb())))/itv.diam()).ub();

		TEMP1 = res_itv.lb()-alpha*Interval(itv.lb());
		TEMP2 = res_itv.ub()-alpha*Interval(itv.ub());
		if (TEMP1.ub()>TEMP2.ub()) {
			// u = 0
			band = Interval(0.0,TEMP1.ub());
		}
		else {
			band = Interval(0.0,TEMP2.ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}




template<class T>
inline AffineMain<T>& AffineMain<T>::Apow(int n, const Interval& itv) {
	//	std::cout << "in power "<<std::endl;

	// Particular case
	if (itv.is_empty()||itv.is_unbounded()) {
		*this = pow(itv,n);
	} else if (!is_actif()) {
		*this = pow(itv,n);
	} else if (itv.diam()< AF_EC) {
		*this = pow(itv,n);
	} else {
	// General Case
		if (n == 0) {
			*this = Interval::one();
			return *this;
		} else if (n == 1)
			return *this;
		else if (n == 2)
			return this->Asqr(itv);
		else if (n<0) {
			this->Apow(-n,itv);
			return this->Ainv(pow(itv,-n));

		} else if (n % 2 == 0) {
			// alpha = (F(sup(x)) - F(inf(x)))/diam(X)
			// u = (f')^{-1}(alpha)
			// d_a = f(inf(x)) -alpha*inf(X)
			// d_b = f(sup(x)) -alpha*sup(x)
			// d_min = min(d_a,d_b)
			// d_max = f(u) - alpha*u
			// beta = Interval(d_min,d_max).mid()
			// zeta = Interval(d_min,d_max).rad()
			double alpha, beta, ddelta, t1, t2;
			Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);

			dmm = pow(itv, n);
			alpha = ((__builtin_powi(itv.ub(),n)-__builtin_powi(itv.lb(),n))/itv.diam());

			TEMP1 = (dmm.lb()) - alpha * Interval(itv.lb());
			TEMP2 = (dmm.ub()) - alpha * Interval(itv.ub());
			// u = (alpha/n)^(1/(n-1))
			if (TEMP1.ub() > TEMP2.ub()) {
				TEMP2 = Interval(alpha) / n;
				band = Interval(
						((1 - n) * TEMP2 * (root(TEMP2, n - 1))).lb(),
						TEMP1.ub());
			} else {
				TEMP1 = Interval(alpha) / n;
				band = Interval(
						((1 - n) * TEMP1 * (root(TEMP1, n - 1))).lb(),
						TEMP2.ub());
			}

			beta = band.mid();
			t1 = (beta - band).ub();
			t2 = (band - beta).ub();
			ddelta = (t1 > t2) ? t1 : t2;

			*this *= alpha;
			*this += beta;
			this->inflate(ddelta);
			//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

		} else {
			// for _itv = [a,b]
			// x0 = 1/sqrt(2)
			// x1= - x0
			// xb0 = 0.5*((b-a)*x0 +(a+b))
			// xb1 = 0.5*((b-a)*x1 +(a+b))
			// c0 = 0.5 (f(xb0)+f(xb1))
			// c1 = x0*f(xb0)+x1*f(xb1)
			// alpha = 2*c1/(b-a)
			// beta = c0-c1*(a+b)/(b-a)
			//  old : ddelta = (b-a)^2 * f''(_itv)/16
			//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

			Interval  TEMP1, TEMP2;
			double alpha, beta, ddelta, t1, x0, xb0, xb1, fxb0, fxb1, c0, c1;

			x0  = 1.0 / ::sqrt(2.);
			xb0 = (0.5) * (itv.diam() * ( x0) + itv.lb() + itv.ub());
			xb1 = (0.5) * (itv.diam() * (-x0) + itv.lb() + itv.ub());
			fxb0 = __builtin_powi(xb0, n);
			fxb1 = __builtin_powi(xb1, n);
			c0 = (0.5) * (fxb0 + fxb1);
			c1 = x0 * fxb0 - x0 * fxb1;

			alpha = 2 * c1 / (itv.diam());
			beta  = c0 - c1 * ((itv.lb() + itv.ub()) / (itv.diam()));
			//ddelta = ((_n * Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub(); //

			// compute the maximal error

			// compute the error at _itv.lb() and _itv.ub()
			ddelta = (abs(
					pow(Interval(itv.lb()), n)
					- (alpha * Interval(itv.lb()) + beta))).ub();
			t1 = (abs(
					pow(Interval(itv.ub()), n)
					- (alpha * Interval(itv.ub()) + beta))).ub();
			if (t1 > ddelta) ddelta= t1 ;
			// u = (alpha/n)^(1/(n-1))
			TEMP2 = pow(Interval(alpha) / n, 1.0 / Interval(n - 1));
			if (!((TEMP2 & itv).is_empty())) {
				t1 = (abs(pow(TEMP2, n) - (alpha * TEMP2 + beta))).ub();
				if (t1 > ddelta) ddelta= t1 ;
			}
			if (!(((-TEMP2) & itv).is_empty())) {
				t1 = (abs(pow(-TEMP2, n) - (alpha * (-TEMP2) + beta))).ub();
				if (t1 > ddelta) ddelta= t1 ;
			}

			*this *= alpha;
			*this += beta;
			this->inflate(ddelta);
			//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);
		}

	}
	//	std::cout << "out power "<<std::endl;
	return *this;
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Apow(double d, const Interval& itv) {
	if ( ((int) (d)) == d) {
		this->Apow((int) (d),itv);
		return *this;
	} else if (d<0) {
		this->Apow(Interval(-d), itv);
		this->Ainv(pow(itv,-d));
		return *this;
	} else {
		this->Apow(Interval(d),itv);
		return *this;
	}
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Aroot(int n, const Interval& itv) {

	if (is_empty()) return *this;
	else if (n==0)  return *this = Interval::one();
	else if (n==1)  return *this;
	else if (is_degenerated()) {
		return *this = pow(Interval(mid()),1.0/n);
	}
	else if (n<0) {
		this->Aroot(-n,itv);
		this->Ainv(root(itv,-n));
		return *this;
	}
	else if (n % 2 == 0) {
		this->Apow(Interval::one()/n,itv);
		return *this; // the negative part of x should be removed
	}
	else if (0 <= itv.lb()) {
		this->Apow(Interval::one()/n,itv);
		return  *this;
	}
	else if (itv.ub() <= 0) {
		this->Aneg();
		this->Apow(Interval::one()/n,-itv);
		this->Aneg();
		return  *this;
	}
	else {
		// TODO do the root when x contains ZERO more properly
		//ibex_error("warning: Affine2Main ROOT non completely well implemented yet");
		//		y=pow(x,e) |  // the negative part of x should be removed
		//	    (-pow(-x,e)); // the positive part of x should be removed
		// BE CAREFULL the result of this union is an INTERVAL, so y lost all its affine form
		return *this = ((pow(itv & Interval::pos_reals(), Interval::one()/n)) | (-pow(-(itv & Interval::neg_reals()),Interval::one()/n)));  // BE CAREFULL the result of this union is an INTERVAL, so y lost all its affine form
		// BE CAREFULL the result of this union is an INTERVAL, so y lost all its affine form
	}

}

//////////////////////////////////////////////////////////////////////////////////
// MinRange Formula


/// for _itv = [a,b]
//if f increase,
// alpha = f'(a)
// band = Interval(a- alpha*a , b-alpha*b)
//if f decrease
// alpha = f'(b)
// band = Interval(a- alpha*b , b-alpha*a)
//
// beta = band.mid()
//  ddelta = band.rad()
//

template<class T>
inline AffineMain<T>& AffineMain<T>::Asqrt_MR(const Interval& itv) {

	Interval res_itv = sqrt(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0), itv2;
		//if ((itv.ub() == POS_INFINITY) || (itv.ub()<0) ) {
		//	_err = sqrt(itv);
		//	_actif = false;
		//	delete[] _val;
		//	_val = NULL;
		//}
		//else {
		if (itv.lb()<0) {
			itv2 = Interval(0.0,itv.ub());
		} else {
			itv2 =itv;
		}

		dmm = sqrt(itv2);
		if (itv2.diam()< AF_EC) {
			alpha = 0.0;
			band =dmm;
		}
		else {
			alpha = (1.0/(2.0*dmm)).lb();  // compute the derivative
			if (alpha<=0) {
				alpha = 0.0;
				band = dmm;
			}
			else {
				TEMP1 = Interval(dmm.lb()) -alpha*itv2.lb();
				TEMP2 = Interval(dmm.ub()) -alpha*itv2.ub();
				if (TEMP1.lb()>TEMP2.ub()) {
					band = Interval(TEMP2.lb(),TEMP1.ub());
					// normally this case never happen
				}
				else {
					band = Interval(TEMP1.lb(),TEMP2.ub());
				}
			}
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band-beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aexp_MR(const Interval& itv) {

	Interval res_itv = exp(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0), itv2;


			dmm = res_itv;
			if (itv.diam()< AF_EC) {
				alpha = 0.0;
				band =dmm;
			}
			else {
				alpha = dmm.lb();  // compute the derivative
				if (alpha<=0) {
					alpha = 0.0;
					band = dmm;
				}
				else {
					TEMP1 = Interval(dmm.lb()) -alpha*itv.lb();
					TEMP2 = Interval(dmm.ub()) -alpha*itv.ub();
					if (TEMP1.lb()>TEMP2.ub()) {
						band = Interval(TEMP2.lb(),TEMP1.ub());
						// normally this case never happen
					}
					else {
						band = Interval(TEMP1.lb(),TEMP2.ub());
					}
				}
			}

			beta = band.mid();
			t1 = (beta -band).ub();
			t2 = (band-beta).ub();
			ddelta = (t1>t2)? t1 : t2;

			*this *= alpha;
			*this += beta;
			this->inflate(ddelta);
			//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);
		}


	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Alog_MR(const Interval& itv) {

	Interval res_itv = log(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0), itv2;

		//if ((itv.lb()<= 0) ||(itv.ub()==POS_INFINITY) ) {
		//	_err = log(itv);
		//	_actif = false;
		//	delete[] _val;
		//	_val = NULL;
		//}
		//else {

			dmm = res_itv;
			if (itv.diam()< AF_EC) {
				alpha = 0.0;
				band =dmm;
			}
			else {
				alpha = (1.0/itv).lb();  // compute the derivative
				if (alpha<=0) {
					alpha = 0.0;
					band = dmm;
				}
				else {
					TEMP1 = Interval(dmm.lb()) -alpha*itv.lb();
					TEMP2 = Interval(dmm.ub()) -alpha*itv.ub();
					if (TEMP1.lb()>TEMP2.ub()) {
						band = Interval(TEMP2.lb(),TEMP1.ub());
						// normally this case never happen
					}
					else {
						band = Interval(TEMP1.lb(),TEMP2.ub());
					}
				}
			}

			beta = band.mid();
			t1 = (beta -band).ub();
			t2 = (band-beta).ub();
			ddelta = (t1>t2)? t1 : t2;

			*this *= alpha;
			*this += beta;
			this->inflate(ddelta);
			//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

		//}

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Ainv_MR(const Interval& itv) {

	Interval res_itv = (1.0/itv);

	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0), itv2;

		dmm = (1.0/abs(itv));
		if (itv.diam()< AF_EC) {
			alpha = 0.0;
			band = dmm;
		}
		else {
			alpha = ((-1.0)/pow(itv,2)).ub(); // compute the derivative
			if (alpha<=0) {
				alpha = 0.0;
				band = dmm;
			}
			else {
				TEMP1 = Interval(dmm.lb()) -alpha*itv.ub();
				TEMP2 = Interval(dmm.ub()) -alpha*itv.lb();
				if (TEMP1.lb()>TEMP2.ub()) {
					band = Interval(TEMP2.lb(),TEMP1.ub());
					// normally this case never happens
				}
				else {
					band = Interval(TEMP1.lb(),TEMP2.ub());
				}
			}
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band-beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		if (itv.lb()<0.0) beta = -beta;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}




/////////////////////////////////////////////////////////
// Specialization



template<>
int AffineMain<AF_fAFFullI>::size() const;

template<>
AffineMain<AF_fAFFullI>& AffineMain<AF_fAFFullI>::operator+=(const Interval& x);

template<>
std::ostream& operator<<(std::ostream& os, const AffineMain<AF_fAFFullI>& x);



#ifdef _IBEX_WITH_AFFINE_EXTENDED_

template<>
AffineMain<AF_No>& AffineMain<AF_No>::Apow(int n, const Interval& itv);

#endif

} // end namespace ibex


#endif /* IBEX_Affine_H_ */


/** \brief atan2(AF[y],AF[x]). */
//Affine2 atan2(const Affine2& y, const Affine2& x);
/** \brief atan2([y],AF[x]). */
//Affine2 atan2(const Interval& y, const Affine2& x);
/** \brief atan2(AF[y],[x]). */
//Affine2 atan2(const Affine2& y, const Interval& x);
/** \brief cosh(AF[x]). */
