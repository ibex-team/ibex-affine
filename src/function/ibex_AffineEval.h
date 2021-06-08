/* ============================================================================
 * I B E X - Evaluation of Affine forms
 * ============================================================================
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Jordan Ninin
 * Created     : April 08, 2013
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_AFFINE_EVAL_H__
#define __IBEX_AFFINE_EVAL_H__

#include "ibex_Function.h"
#include "ibex_AffineVar.h"
#include "ibex_AffineMatrix.h"
#include "ibex_AffineDomain.h"
#include "ibex_FwdAlgorithm.h"
#include "ibex_Domain.h"
#include "ibex_NodeMap.h"


namespace ibex {

typedef AffineEval<AF_Default> Affine2Eval;
typedef AffineEval<AF_Other>  Affine3Eval;

/**
 * \ingroup symbolic
 *
 * \brief Evaluate a function with affine form.
 */
template<class T>
class AffineEval : public FwdAlgorithm {

public:
	/**
	 * \brief Build an Affine evaluator for f.
	 */
	AffineEval(const Function& f);

	/**
	 * \brief Delete this.
	 */
	~AffineEval();

//	/**
//	 * \brief Run the forward algorithm with input domains.
//	 */
//	Domain& eval(const Array<const Domain>& d);
//
//	/**
//	 * \brief Run the forward algorithm with input domains.
//	 */
//	Domain& eval(const Array<Domain>& d);

	/**
	 * \brief Run the forward algorithm with an input box.
	 */
	Domain& eval(const IntervalVector& box);

	/**
	 * \brief Evaluate a subset of components.
	 *
	 * (Specific for vector-valued functions).
	 *
	 * \pre components must be non empty and contain indices in
	 *      [0,f.image_dim()-1] in case of a vector-valued function
	 *      [0,f.expr().dim.nb_rows()] in case of a matrix-valued function.
	 */
	Domain eval(const IntervalVector& box, const BitSet& components);

	/**
	 * \brief Evaluate a submatrix.
	 *
	 * (Specific for matrix_valued functions).
	 *
	 */
	Domain eval(const IntervalVector& box, const BitSet& rows, const BitSet& cols);

	/**
	 * \brief Run the forward algorithm with input domains, and return the result as an Affine domain.
	 */
	TemplateDomain<AffineMain<T> >& eval(const Array<const Domain>& dom, const Array<const TemplateDomain<AffineMain<T> > >& aff);

	/**
	 * \brief Run the forward algorithm with an input \box and return the result as an Affine domain.
	 */
	TemplateDomain<AffineMain<T> >& eval(const AffineMainVector<T>& box);
	TemplateDomain<AffineMain<T> >& eval(const AffineVarMainVector<T>& box);
	std::pair<Domain*,TemplateDomain<AffineMain<T> >* > eval(const IntervalVector& box,const AffineMainVector<T>& aff);
	std::pair<Domain*,TemplateDomain<AffineMain<T> >* > eval(const IntervalVector& box,const AffineVarMainVector<T>& aff);

	/**
	 * \brief Evaluate a subset of components, and return the result as an Affine domain.
	 *
	 * (Specific for vector-valued functions).
	 *
	 * \pre components must be non empty and contain indices in
	 *      [0,f.image_dim()-1] in case of a vector-valued function
	 *      [0,f.expr().dim.nb_rows()] in case of a matrix-valued function.
	 */
	TemplateDomain<AffineMain<T> > eval(const AffineMainVector<T>& box, const BitSet& components);
	TemplateDomain<AffineMain<T> > eval(const AffineVarMainVector<T>& box, const BitSet& components);

	/**
	 * \brief Evaluate a submatrix, and return the result as an Affine domain.
	 *
	 * (Specific for matrix_valued functions).
	 *
	 */
	TemplateDomain<AffineMain<T> > eval(const AffineMainVector<T>& box, const BitSet& rows, const BitSet& cols);
	TemplateDomain<AffineMain<T> > eval(const AffineVarMainVector<T>& box, const BitSet& rows, const BitSet& cols);


	void idx_cp_fwd (int x, int y);
	void idx_fwd (int  x, int y);
	void vector_fwd(int* x, int y);
	void cst_fwd   (int y);
	void symbol_fwd(int y);
	void apply_fwd (int* x, int y);
	void chi_fwd   (int x1, int x2, int x3, int y);
	void gen2_fwd  (int x1, int x2, int y);
	void add_fwd   (int x1, int x2, int y);
	void mul_fwd   (int x1, int x2, int y);
	void sub_fwd   (int x1, int x2, int y);
	void div_fwd   (int x1, int x2, int y);
	void max_fwd   (int x1, int x2, int y);
	void min_fwd   (int x1, int x2, int y);
	void atan2_fwd (int x1, int x2, int y);
	void gen1_fwd  (int x, int y);
	void minus_fwd (int x, int y);
	void minus_V_fwd(int x, int y);
	void minus_M_fwd(int x, int y);
	void trans_V_fwd(int x, int y);
	void trans_M_fwd(int x, int y);
	void sign_fwd   (int x, int y);
	void abs_fwd    (int x, int y);
	void power_fwd  (int x, int y, int p);
	void sqr_fwd    (int x, int y);
	void sqrt_fwd   (int x, int y);
	void exp_fwd    (int x, int y);
	void log_fwd    (int x, int y);
	void cos_fwd    (int x, int y);
	void sin_fwd    (int x, int y);
	void tan_fwd    (int x, int y);
	void cosh_fwd   (int x, int y);
	void sinh_fwd   (int x, int y);
	void tanh_fwd   (int x, int y);
	void acos_fwd   (int x, int y);
	void asin_fwd   (int x, int y);
	void atan_fwd   (int x, int y);
	void acosh_fwd  (int x, int y);
	void asinh_fwd  (int x, int y);
	void atanh_fwd  (int x, int y);
	void floor_fwd  (int x, int y);
	void ceil_fwd   (int x, int y);
	void saw_fwd    (int x, int y);

	void add_V_fwd  (int x1, int x2, int y);
	void add_M_fwd  (int x1, int x2, int y);
	void mul_SV_fwd (int x1, int x2, int y);
	void mul_SM_fwd (int x1, int x2, int y);
	void mul_VV_fwd (int x1, int x2, int y);
	void mul_MV_fwd (int x1, int x2, int y);
	void mul_VM_fwd (int x1, int x2, int y);
	void mul_MM_fwd (int x1, int x2, int y);
	void sub_V_fwd  (int x1, int x2, int y);
	void sub_M_fwd  (int x1, int x2, int y);

	const Function& f;
	ExprDomain d;
	ExprTemplateDomain<AffineMain<T> > af2;
	Agenda** fwd_agenda;         // one agenda for each vector component/matrix row
	Agenda** bwd_agenda;         // one agenda for each vector component/matrix row
	Agenda*** matrix_fwd_agenda; // one agenda for each matrix element
	Agenda*** matrix_bwd_agenda; // one agenda for each matrix element

protected:

	/**
	 * Since there is no affine evaluator in the Function class,
	 * we store here all the evaluators of the called functions
	 * (this avoids to create this object at each evaluation).
	 */
	NodeMap<AffineEval<T>* > apply_eval;


protected:
	/**
	 * Class used internally to interrupt the forward procedure
	 * when an empty domain occurs (<=> the input box is outside
	 * the definition domain of the function).
	 */
	class EmptyBoxException { };

};


/* ============================================================================
 	 	 	 	 	 	 	 implementation
  ============================================================================*/

template<class T>
AffineEval<T>::AffineEval(const Function& f) :
	f(f),
	d(f),
	af2(f),
	fwd_agenda(NULL), bwd_agenda(NULL), matrix_fwd_agenda(NULL), matrix_bwd_agenda(NULL) {

	Dim dim=f.expr().dim;
	int m=dim.vec_size();

	if (m>1) {
		const ExprVector* vec=dynamic_cast<const ExprVector*>(&f.expr());
		if (vec && (vec->orient==ExprVector::COL || dim.type()==Dim::ROW_VECTOR) && m==vec->nb_args) {
			fwd_agenda = new Agenda*[m];
			bwd_agenda = new Agenda*[m];
			for (int i=0; i<m; i++) {
				bwd_agenda[i] = f.cf.agenda(f.nodes.rank(vec->arg(i)));
				fwd_agenda[i] = new Agenda(*bwd_agenda[i],true); // true<=>swap
			}

			if (dim.is_matrix()) {
				int n=dim.nb_cols();

				// check that the matrix is homogeneous (a matrix of scalar expressions)
				bool homogeneous=true; // by default
				for (int i=0; i<m; i++) {
					const ExprVector* fi=dynamic_cast<const ExprVector*>(&vec->arg(i));
					if (!fi || fi->nb_args<n) {
						homogeneous=false;
						break;
					}
				}

				if (homogeneous) {
					matrix_fwd_agenda = new Agenda**[m];
					matrix_bwd_agenda = new Agenda**[m];
					for (int i=0; i<m; i++) {
						const ExprVector& fi=(const ExprVector&) vec->arg(i);
						matrix_fwd_agenda[i] = new Agenda*[n];
						matrix_bwd_agenda[i] = new Agenda*[n];
						for (int j=0; j<n; j++) {
							matrix_bwd_agenda[i][j] = f.cf.agenda(f.nodes.rank(fi.arg(j)));
							matrix_fwd_agenda[i][j] = new Agenda(*matrix_bwd_agenda[i][j],true); // true<=>swap
						}
					}
				}
			}
		}
	}
}

template<class T>
AffineEval<T>::~AffineEval() {
	if (fwd_agenda!=NULL) {
		for (int i=0; i<f.expr().dim.vec_size(); i++) {
			delete fwd_agenda[i];
			delete bwd_agenda[i];
		}
		delete[] fwd_agenda;
		delete[] bwd_agenda;

		if (matrix_fwd_agenda!=NULL) {
			for (int i=0; i<f.expr().dim.nb_rows(); i++) {
				for (int j=0; j<f.expr().dim.nb_cols(); j++) {
					delete matrix_fwd_agenda[i][j];
					delete matrix_bwd_agenda[i][j];
				}
				delete[] matrix_fwd_agenda[i];
				delete[] matrix_bwd_agenda[i];
			}
			delete[] matrix_fwd_agenda;
			delete[] matrix_bwd_agenda;
		}
	}

	for (typename IBEX_NODE_MAP(AffineEval<T>* )::iterator it=apply_eval.begin(); it!=apply_eval.end(); it++) {
		delete it->second;
	}
}

template<class T>
Domain& AffineEval<T>::eval(const IntervalVector& box) {
	d.write_arg_domains(box);
	af2.write_arg_domains(AffineMainVector<T>(AffineVarMainVector<T>(box)));

	try {
		f.forward<AffineEval<T> >(*this);
	} catch(EmptyBoxException&) {
		d.top->set_empty();
		af2.top->set_empty();
	}
	return *d.top;
}
template<class T>
TemplateDomain<AffineMain<T> >& AffineEval<T>::eval(const AffineMainVector<T>& box) {
	d.write_arg_domains(box.itv());
	af2.write_arg_domains(box);

	try {
		f.forward<AffineEval<T> >(*this);
	} catch(EmptyBoxException&) {
		d.top->set_empty();
		af2.top->set_empty();
	}
	return *af2.top;
}
template<class T>
inline TemplateDomain<AffineMain<T> >& AffineEval<T>::eval(const AffineVarMainVector<T>& box) {
	return eval(AffineMainVector<T>(box));
}


template<class T>
std::pair<Domain*,TemplateDomain<AffineMain<T> >* > AffineEval<T>::eval(const IntervalVector& box,const AffineMainVector<T>& aff) {
	d.write_arg_domains(box);
	af2.write_arg_domains(aff);

	try {
		f.forward<AffineEval<T> >(*this);
	} catch(EmptyBoxException&) {
		d.top->set_empty();
		af2.top->set_empty();
	}
	return std::pair<Domain*,TemplateDomain<AffineMain<T> >* >(d.top,af2.top);
}
template<class T>
std::pair<Domain*,TemplateDomain<AffineMain<T> >* > AffineEval<T>::eval(const IntervalVector& box,const AffineVarMainVector<T>& aff) {
	return eval(box,AffineMainVector<T>(aff));
}


template<class T>
TemplateDomain<AffineMain<T> >& AffineEval<T>::eval(const Array<const Domain>& box, const Array<const TemplateDomain<AffineMain<T> > >& aff){
	d.write_arg_domains(box);
	af2.write_arg_domains(aff);

	try {
		f.forward<AffineEval<T> >(*this);
	} catch(EmptyBoxException&) {
		d.top->set_empty();
		af2.top->set_empty();
	}
	return *af2.top;
}




template<class T>
Domain AffineEval<T>::eval(const IntervalVector& box, const BitSet& components) {

	Dim dim=d.top->dim;

	if (dim.type()==Dim::SCALAR) return eval(box);

	d.write_arg_domains(box);
	af2.write_arg_domains(AffineMainVector<T>(AffineVarMainVector<T>(box)));

	assert(!components.empty());

	int m=components.size();

	Domain res(dim.type()==Dim::ROW_VECTOR ?
			Dim(1,m) : // in the case of a row vector, we select columns
			Dim(m,dim.nb_cols())); // in the other cases we select rows

	if (fwd_agenda==NULL) {

		// The vector of expression is heterogeneous (or the expression is scalar).
		//
		// We might be able in this case to use the DAG but
		// - the algorithm is more complex
		// - we might not benefit from possible symbolic simplification due
		//   to the fact that only specific components are required (there is
		//   no simple "on the fly" simplification as in the case of a vector
		//   of homogeneous expressions)
		// so we resort to the components functions f[i] --> symbolic copy+no DAG :(
		int i=0;
		for (BitSet::const_iterator c=components.begin(); c!=components.end(); ++c) {


			AffineEval<T> *func_eval;
			if (!apply_eval.found(f.nodes[c])) {
				func_eval=new AffineEval<T>(f[c]);
				apply_eval.insert(f.nodes[c],func_eval);
			} else {
				func_eval=apply_eval[f.nodes[c]];
			}

			res[i++] = func_eval->eval(box);
		}
		assert(i==m);

		return res;
	}

	// merge all the agendas
	Agenda a(f.nodes.size()); // the global agenda initialized with the maximal possible value
	for (BitSet::const_iterator c=components.begin(); c!=components.end(); ++c) {
		a.push(*(fwd_agenda[c]));
	}

	try {
		f.cf.forward<AffineEval<T> >(*this,a);
		int i=0;
		for (BitSet::const_iterator c=components.begin(); c!=components.end(); ++c) {
			res[i++] = d[bwd_agenda[c]->first()];
		}
	} catch(EmptyBoxException&) {
		d.top->set_empty();
		af2.top->set_empty();
		res.set_empty();
	}

	return res;
}


template<class T>
Domain AffineEval<T>::eval(const IntervalVector& box, const BitSet& rows, const BitSet& cols) {

	Dim dim=d.top->dim;

	switch (dim.type()) {
	case Dim::SCALAR:     return eval(box);
	case Dim::ROW_VECTOR: return eval(box,cols);
	case Dim::COL_VECTOR: return eval(box,rows);
	default : ;// ok continue
	}

	d.write_arg_domains(box);
	af2.write_arg_domains(AffineMainVector<T>(AffineVarMainVector<T>(box)));

	assert(!rows.empty());
	assert(!cols.empty());
	assert(rows.max()<dim.nb_rows());
	assert(cols.max()<dim.nb_cols());

	int m=rows.size();
	int n=cols.size();

	Domain res(Dim(m,n));

	if (matrix_fwd_agenda==NULL) {

		// The vector of expression is heterogeneous (or the expression is scalar).
		//
		// We might be able in this case to use the DAG but
		// - the algorithm is more complex
		// - we might not benefit from possible symbolic simplification due
		//   to the fact that only specific components are required (there is
		//   no simple "on the fly" simplification as in the case of a vector
		//   of homogeneous expressions)
		// so we resort to the components functions f[i] --> symbolic copy+no DAG :(
		int i=0;
		for (BitSet::const_iterator r=rows.begin(); r!=rows.end(); ++r, i++) {
				int j=0;
				for (BitSet::const_iterator c=cols.begin(); c!=cols.end(); ++c, j++) {

					AffineEval<T> *func_eval;
					if (!apply_eval.found(f.nodes[i][j])) {
						func_eval=new AffineEval<T>(f[i][j]);
						apply_eval.insert(f.nodes[i][j],func_eval);
					} else {
						func_eval=apply_eval[f.nodes[i][j]];
					}

					res[i][j] = func_eval->eval(box);
					//res[i][j] = f[r][c].eval_domain(box);
				}
				assert(j==n);
		}
		assert(i==m);

		return res;
	}

	// merge all the agendas
	Agenda a(f.nodes.size()); // the global agenda initialized with the maximal possible value
	for (BitSet::const_iterator r=rows.begin(); r!=rows.end(); ++r) {
		for (BitSet::const_iterator c=cols.begin(); c!=cols.end(); ++c) {
			a.push(*(matrix_fwd_agenda[r][c]));
		}
	}

	try {
		f.cf.forward<AffineEval<T> >(*this,a);
		int i=0;
		for (BitSet::const_iterator r=rows.begin(); r!=rows.end(); ++r, i++) {
			int j=0;
			for (BitSet::const_iterator c=cols.begin(); c!=cols.end(); ++c, j++) {
				res[i][j] = d[matrix_bwd_agenda[r][c]->first()];
			}
		}
	} catch(EmptyBoxException&) {
		d.top->set_empty();
		af2.top->set_empty();
		res.set_empty();
	}

	return res;
}


template<class T>
TemplateDomain<AffineMain<T> > AffineEval<T>::eval(const AffineMainVector<T>& aff, const BitSet& components) {

	Dim dim=d.top->dim;

	if (dim.type()==Dim::SCALAR) return eval(aff);

	d.write_arg_domains(aff.itv());
	af2.write_arg_domains(aff);

	assert(!components.empty());

	int m=components.size();

	TemplateDomain<AffineMain<T> > res(dim.type()==Dim::ROW_VECTOR ?
			Dim(1,m) : // in the case of a row vector, we select columns
			Dim(m,dim.nb_cols())); // in the other cases we select rows

	if (fwd_agenda==NULL) {

		// The vector of expression is heterogeneous (or the expression is scalar).
		//
		// We might be able in this case to use the DAG but
		// - the algorithm is more complex
		// - we might not benefit from possible symbolic simplification due
		//   to the fact that only specific components are required (there is
		//   no simple "on the fly" simplification as in the case of a vector
		//   of homogeneous expressions)
		// so we resort to the components functions f[i] --> symbolic copy+no DAG :(
		int i=0;
		for (BitSet::const_iterator c=components.begin(); c!=components.end(); ++c) {


			AffineEval<T> *func_eval;
			if (!apply_eval.found(f.nodes[c])) {
				func_eval=new AffineEval<T>(f[c]);
				apply_eval.insert(f.nodes[c],func_eval);
			} else {
				func_eval=apply_eval[f.nodes[c]];
			}

			res[i++] = func_eval->eval(aff);
		}
		assert(i==m);

		return res;
	}

	// merge all the agendas
	Agenda a(f.nodes.size()); // the global agenda initialized with the maximal possible value
	for (BitSet::const_iterator c=components.begin(); c!=components.end(); ++c) {
		a.push(*(fwd_agenda[c]));
	}

	try {
		f.cf.forward<AffineEval<T> >(*this,a);
		int i=0;
		for (BitSet::const_iterator c=components.begin(); c!=components.end(); ++c) {
			res[i++] = af2[bwd_agenda[c]->first()];
		}
	} catch(EmptyBoxException&) {
		d.top->set_empty();
		af2.top->set_empty();
		res.set_empty();
	}

	return res;

}
template<class T>
TemplateDomain<AffineMain<T> > AffineEval<T>::eval(const AffineVarMainVector<T>& aff, const BitSet& components) {
	return eval(AffineMainVector<T>(aff), components);
}


template<class T>
TemplateDomain<AffineMain<T> > AffineEval<T>::eval(const AffineMainVector<T>& aff, const BitSet& rows, const BitSet& cols) {

	Dim dim=d.top->dim;

	switch (dim.type()) {
	case Dim::SCALAR:     return eval(aff);
	case Dim::ROW_VECTOR: return eval(aff,cols);
	case Dim::COL_VECTOR: return eval(aff,rows);
	default : ;// ok continue
	}

	d.write_arg_domains(aff.itv());
	af2.write_arg_domains(aff);

	assert(!rows.empty());
	assert(!cols.empty());
	assert(rows.max()<dim.nb_rows());
	assert(cols.max()<dim.nb_cols());

	int m=rows.size();
	int n=cols.size();

	Domain res(Dim(m,n));

	if (matrix_fwd_agenda==NULL) {

		// The vector of expression is heterogeneous (or the expression is scalar).
		//
		// We might be able in this case to use the DAG but
		// - the algorithm is more complex
		// - we might not benefit from possible symbolic simplification due
		//   to the fact that only specific components are required (there is
		//   no simple "on the fly" simplification as in the case of a vector
		//   of homogeneous expressions)
		// so we resort to the components functions f[i] --> symbolic copy+no DAG :(
		int i=0;
		for (BitSet::const_iterator r=rows.begin(); r!=rows.end(); ++r, i++) {
				int j=0;
				for (BitSet::const_iterator c=cols.begin(); c!=cols.end(); ++c, j++) {

					AffineEval<T> *func_eval;
					if (!apply_eval.found(f.nodes[i][j])) {
						func_eval=new AffineEval<T>(f[i][j]);
						apply_eval.insert(f.nodes[i][j],func_eval);
					} else {
						func_eval=apply_eval[f.nodes[i][j]];
					}

					res[i][j] = func_eval->eval(aff);
					//res[i][j] = f[r][c].eval_domain(box);
				}
				assert(j==n);
		}
		assert(i==m);

		return res;
	}

	// merge all the agendas
	Agenda a(f.nodes.size()); // the global agenda initialized with the maximal possible value
	for (BitSet::const_iterator r=rows.begin(); r!=rows.end(); ++r) {
		for (BitSet::const_iterator c=cols.begin(); c!=cols.end(); ++c) {
			a.push(*(matrix_fwd_agenda[r][c]));
		}
	}

	try {
		f.cf.forward<AffineEval<T> >(*this,a);
		int i=0;
		for (BitSet::const_iterator r=rows.begin(); r!=rows.end(); ++r, i++) {
			int j=0;
			for (BitSet::const_iterator c=cols.begin(); c!=cols.end(); ++c, j++) {
				res[i][j] = af2[matrix_bwd_agenda[r][c]->first()];
			}
		}
	} catch(EmptyBoxException&) {
		d.top->set_empty();
		af2.top->set_empty();
		res.set_empty();
	}

	return res;
}


template<class T>
TemplateDomain<AffineMain<T> > AffineEval<T>::eval(const AffineVarMainVector<T>& aff, const BitSet& rows, const BitSet& cols) {
	return eval(AffineMainVector<T>(aff), rows, cols);
}


template<class T>
inline void AffineEval<T>::idx_fwd(int, int ) { /* nothing to do */ }


template<class T>
void AffineEval<T>::idx_cp_fwd(int x, int y) {
	assert(dynamic_cast<const ExprIndex*> (&f.node(y)));

	const ExprIndex& e = (const ExprIndex&) f.node(y);

	d[y] = d[x][e.index];
	af2[y] = af2[x][e.index];
}

template<class T>
inline void AffineEval<T>::symbol_fwd(int) { /* nothing to do */ }

template<class T>
inline void AffineEval<T>::cst_fwd(int y) {
	const ExprConstant& c = (const ExprConstant&) f.node(y);
	switch (c.type()) {
	case Dim::SCALAR:      {
		af2[y].i() = c.get_value();
		d[y].i() = c.get_value();
		break;
	}
	case Dim::ROW_VECTOR:
	case Dim::COL_VECTOR: {
		af2[y].v() = c.get_vector_value();
		d[y].v() = c.get_vector_value();
		break;
	}
	case Dim::MATRIX: {
		af2[y].m() = c.get_matrix_value();
		d[y].m() = c.get_matrix_value();
		break;
	}
	}
}
template<class T>
inline void AffineEval<T>::gen2_fwd(int x1, int x2, int y) {
	const ExprGenericBinaryOp& e = (const ExprGenericBinaryOp&) f.node(y);
	d[y]=e.eval(d[x1],d[x2]);
	switch(e.dim.type()) {
	case Dim::SCALAR :     af2[y].i()= (d[y].i()); break;
	case Dim::ROW_VECTOR :
	case Dim::COL_VECTOR : af2[y].v()=(d[y].v()); break;
	case Dim::MATRIX :     af2[y].m()=(d[y].m()); break;
	}
}

template<class T>
inline void AffineEval<T>::chi_fwd(int x1, int x2, int x3, int y) {
	af2[y].i()=chi(d[x1].i(),af2[x2].i(),af2[x3].i());
	d[y].i()  =chi(d[x1].i(),d[x2].i(),d[x3].i());
}

template<class T>
inline void AffineEval<T>::add_fwd(int x1, int x2, int y) {
	af2[y].i()=af2[x1].i()+af2[x2].i();
	d[y].i()=(af2[y].i().itv() & (d[x1].i()+d[x2].i()));
}

template<class T>
inline void AffineEval<T>::mul_fwd(int x1, int x2, int y) {
	af2[y].i()=af2[x1].i()*af2[x2].i();
	d[y].i()=(af2[y].i().itv() & (d[x1].i()*d[x2].i()));
}

template<class T>
inline void AffineEval<T>::sub_fwd(int x1, int x2, int y) {
	af2[y].i()=af2[x1].i()-af2[x2].i();
	d[y].i()=(af2[y].i().itv() & (d[x1].i()-d[x2].i()));
}

template<class T>
inline void AffineEval<T>::div_fwd(int x1, int x2, int y) {
	af2[y].i()=af2[x2].i();
	af2[y].i().Ainv(d[x2].i());
	af2[y].i() *= af2[x1].i();
	d[y].i()=(af2[y].i().itv() & (d[x1].i()/d[x2].i()));
}

template<class T>
inline void AffineEval<T>::max_fwd(int x1, int x2, int y) {
	d[y].i()= max(d[x1].i(),d[x2].i());
	af2[y].i() = (d[y].i()); // Not implemented in Affine Arithmetic
}

template<class T>
inline void AffineEval<T>::min_fwd(int x1, int x2, int y) {
	d[y].i() = min(d[x1].i(),d[x2].i());
	af2[y].i()= ( d[y].i()); // Not implemented in Affine Arithmetic
}

template<class T>
inline void AffineEval<T>::atan2_fwd(int x1, int x2, int y) {
	d[y].i() = atan2(d[x1].i(),d[x2].i());
	af2[y].i()= (d[y].i()); // Not implemented in Affine Arithmetic
}

template<class T>
inline void AffineEval<T>::gen1_fwd(int x, int y) {
	const ExprGenericUnaryOp& e = (const ExprGenericUnaryOp&) f.node(y);
	d[y]=e.eval(d[x]);
	switch(e.dim.type()) {
	case Dim::SCALAR :     af2[y].i()= (d[y].i()); break;
	case Dim::ROW_VECTOR :
	case Dim::COL_VECTOR : af2[y].v()=(d[y].v()); break;
	case Dim::MATRIX :     af2[y].m()=(d[y].m()); break;
	}
}

template<class T>
inline void AffineEval<T>::minus_fwd(int x, int y) {
	af2[y].i()=-af2[x].i();
	d[y].i()=(af2[y].i().itv() & (-d[x].i()));
}

template<class T>
inline void AffineEval<T>::minus_V_fwd(int x, int y) {
	af2[y].v()=-af2[x].v();
	d[y].v()=(af2[y].v().itv() & (-d[x].v()));
}

template<class T>
inline void AffineEval<T>::minus_M_fwd(int x, int y) {
	af2[y].m()=-af2[x].m();
	d[y].m() = af2[y].m().itv();
	d[y].m()&= (-d[x].m());
}

template<class T>
inline void AffineEval<T>::sign_fwd(int x, int y) {
	d[y].i()=sign(d[x].i());
	af2[y].i()=d[y].i();
}

template<class T>
inline void AffineEval<T>::abs_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Aabs(d[x].i());
	d[y].i()=(af2[y].i().itv() & abs(d[x].i()));
}

template<class T>
inline void AffineEval<T>::power_fwd(int x, int y, int p) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Apow(p,d[x].i());
	d[y].i()=(af2[y].i().itv() & pow(d[x].i(),p));
}

template<class T>
inline void AffineEval<T>::sqr_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Asqr(d[x].i());
	d[y].i()=(af2[y].i().itv() & sqr(d[x].i()));
}

template<class T>
inline void AffineEval<T>::sqrt_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Asqrt(d[x].i());
	d[y].i()=(af2[y].i().itv() & sqrt(d[x].i()));
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::exp_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Aexp(d[x].i());
	d[y].i()=(af2[y].i().itv() & exp(d[x].i()));
}

template<class T>
inline void AffineEval<T>::log_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Alog(d[x].i());
	d[y].i()=(af2[y].i().itv() & log(d[x].i()));
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::cos_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Acos(d[x].i());
	d[y].i()=(af2[y].i().itv() & cos(d[x].i()));
}

template<class T>
inline void AffineEval<T>::sin_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Asin(d[x].i());
	d[y].i()=(af2[y].i().itv() & sin(d[x].i()));
}

template<class T>
inline void AffineEval<T>::tan_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Atan(d[x].i());
	d[y].i()=(af2[y].i().itv() & tan(d[x].i()));
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::cosh_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Acosh(d[x].i());
	d[y].i()=(af2[y].i().itv() & cosh(d[x].i()));
}

template<class T>
inline void AffineEval<T>::sinh_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Asinh(d[x].i());
	d[y].i()=(af2[y].i().itv() & sinh(d[x].i()));
}

template<class T>
inline void AffineEval<T>::tanh_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Atanh(d[x].i());
	d[y].i()=(af2[y].i().itv() & tanh(d[x].i()));
}

template<class T>
inline void AffineEval<T>::acos_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Aacos(d[x].i());
	d[y].i()=(af2[y].i().itv() & acos(d[x].i()));
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::asin_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Aasin(d[x].i());
	d[y].i()=(af2[y].i().itv() & asin(d[x].i()));
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::atan_fwd(int x, int y) {
	af2[y].i()=AffineMain<T>(af2[x].i()).Aatan(d[x].i());
	d[y].i()=(af2[y].i().itv() & atan(d[x].i()));
}

template<class T>
inline void AffineEval<T>::acosh_fwd(int x, int y) {
	d[y].i()=acosh(d[x].i());
	af2[y].i()= (d[y].i()); // Not implemented in Affine Arithmetic
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::asinh_fwd(int x, int y) {
	d[y].i()=asinh(d[x].i());
	af2[y].i()= (d[y].i()); // Not implemented in Affine Arithmetic
}

template<class T>
inline void AffineEval<T>::atanh_fwd(int x, int y) {
	d[y].i()=atanh(d[x].i());
	af2[y].i()= (d[y].i()); // Not implemented in Affine Arithmetic
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::floor_fwd(int x, int y) {
	d[y].i()=floor(d[x].i());
	af2[y].i()= (d[y].i()); // Not implemented in Affine Arithmetic
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::ceil_fwd(int x, int y) {
	d[y].i()=ceil(d[x].i());
	af2[y].i()= (d[y].i()); // Not implemented in Affine Arithmetic
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::saw_fwd(int x, int y) {
	d[y].i()=saw(d[x].i());
	af2[y].i()= (d[y].i()); // Not implemented in Affine Arithmetic
	if ((d[y].i()).is_empty()) throw EmptyBoxException();
}

template<class T>
inline void AffineEval<T>::trans_V_fwd(int x, int y) {
	af2[y].v()=af2[x].v();
	d[y].v()=d[x].v();
}

template<class T>
inline void AffineEval<T>::trans_M_fwd(int x, int y) {
	af2[y].m()=af2[x].m().transpose();
	d[y].m()=d[x].m().transpose();
}

template<class T>
inline void AffineEval<T>::add_V_fwd(int x1, int x2, int y) {
	af2[y].v()=af2[x1].v()+af2[x2].v();
	d[y].v()=(af2[y].v().itv() & (d[x1].v()+d[x2].v()));
}

template<class T>
inline void AffineEval<T>::add_M_fwd(int x1, int x2, int y) {
	af2[y].m()=af2[x1].m()+af2[x2].m();
	d[y].m()=(af2[y].m().itv());
	d[y].m() &= (d[x1].m()+d[x2].m());
}

template<class T>
inline void AffineEval<T>::mul_SV_fwd(int x1, int x2, int y) {
	af2[y].v()=af2[x1].i()*af2[x2].v();
	d[y].v()=(af2[y].v().itv() & (d[x1].i()*d[x2].v()));
}

template<class T>
inline void AffineEval<T>::mul_SM_fwd(int x1, int x2, int y) {
	af2[y].m()=af2[x1].i()*af2[x2].m();
	d[y].m()=(af2[y].m().itv());
	d[y].m() &= (d[x1].i()*d[x2].m());
}

template<class T>
inline void AffineEval<T>::mul_VV_fwd(int x1, int x2, int y) {
	af2[y].i()=af2[x1].v()*af2[x2].v();
	d[y].i()=(af2[y].i().itv() & (d[x1].v()*d[x2].v()));
}

template<class T>
inline void AffineEval<T>::mul_MV_fwd(int x1, int x2, int y) {
	af2[y].v()=af2[x1].m()*af2[x2].v();
	d[y].v()=(af2[y].v().itv() & (d[x1].m()*d[x2].v()));
}

template<class T>
inline void AffineEval<T>::mul_VM_fwd(int x1, int x2, int y) {
	af2[y].v()=af2[x1].v()*af2[x2].m();
	d[y].v()=(af2[y].v().itv() & (d[x1].v()*d[x2].m()));
}

template<class T>
inline void AffineEval<T>::mul_MM_fwd(int x1, int x2, int y) {
	af2[y].m()=af2[x1].m()*af2[x2].m();
	d[y].m()=(af2[y].m().itv()) ;
	d[y].m() &=  (d[x1].m()*d[x2].m());
}

template<class T>
inline void AffineEval<T>::sub_V_fwd(int x1, int x2, int y) {
	af2[y].v()=af2[x1].v()-af2[x2].v();
	d[y].v()=(af2[y].v().itv() & (d[x1].v()-d[x2].v()));
}

template<class T>
inline void AffineEval<T>::sub_M_fwd(int x1, int x2, int y) {
	af2[y].m()=af2[x1].m()-af2[x2].m();
	d[y].m() = af2[y].m().itv();
	d[y].m() &= (d[x1].m()-d[x2].m());
}

template<class T>
inline void AffineEval<T>::apply_fwd(int* x, int y) {
	assert(dynamic_cast<const ExprApply*> (&f.node(y)));

	const ExprApply& a = (const ExprApply&) f.node(y);

	assert(&a.func!=&f); // recursive calls not allowed

	Array<const Domain> d2(a.func.nb_arg());
	Array<const TemplateDomain<AffineMain<T> > > af22(a.func.nb_arg());

	for (int i=0; i<a.func.nb_arg(); i++) {
		d2.set_ref(i,d[x[i]]);
		af22.set_ref(i,af2[x[i]]);
	}

	AffineEval<T> *func_eval;
	if (!apply_eval.found(a)) {
		func_eval=new AffineEval<T>(a.func);
		apply_eval.insert(a,func_eval);
	} else {
		func_eval=apply_eval[a];
	}

	func_eval->eval(d2,af22);

	d[y] = *func_eval->d.top;
	af2[y] = *func_eval->af2.top;
}

template<class T>
inline void AffineEval<T>::vector_fwd(int *x, int y) {
	assert(dynamic_cast<const ExprVector*>(&(f.node(y))));

	const ExprVector& v = (const ExprVector&) f.node(y);

	assert(v.type()!=Dim::SCALAR);

	if (v.dim.is_vector()) {
		for (int i=0; i<v.length(); i++)  {
			af2[y].v()[i]=af2[x[i]].i();
			d[y].v()[i]=d[x[i]].i();
		}
	}
	else {
		if (v.row_vector())
			for (int i=0; i<v.length(); i++) {
				af2[y].m().set_col(i,af2[x[i]].v());
				d[y].m().set_col(i,d[x[i]].v());
			}
		else
			for (int i=0; i<v.length(); i++) {
				af2[y].m().set_row(i,af2[x[i]].v());
				d[y].m().set_row(i,d[x[i]].v());
			}
	}
}
} // namespace ibex


#endif /* __IBEX_AFFINE_EVAL_H__ */
