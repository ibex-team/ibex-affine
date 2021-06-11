/* ============================================================================
 * I B E X - Example
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENCE.
 *
 * Author(s)   : Jordan Ninin
 * Created     : June 6, 2020
 * ---------------------------------------------------------------------------- */

#include "ibex.h"
#include "ibex_Affine.h"
#include <time.h>

using namespace std;
using namespace ibex;

int main() {

	{
		cout << "==========================================" << endl;
		cout << "==========================================" << endl;
		cout << "TEST 1: Definitions " << endl;


		IntervalVector I(3, Interval(1, 2));
		I[1] = Interval(1, 3);

		Affine2Variables AF(I);
		cout << "==  Initialization of the Affine2 Form with 3 epsilons variables: " << endl << AF << endl;

		cout << "==  Affine2Variables is the only way to construct an correct Affine2 Form from an Interval: " << endl ;
		Affine2Variables xx(1 , I[1]);
		cout<< xx[0] << endl;
		Affine2Variables yy(1);
		yy[0] = I[1];
		cout<< yy[0] << endl;

		cout << "==  Constructor Affine2(Interval& itv) and Constructor Affine2Variables(Affine2& af) is not allowed to avoid misunderstandings or bad use. " << endl ;
		cout << "==  Affine2 is used for intermediate result or constant: " << endl ;
		Affine2 faa;
		faa = I[1];
		cout << "== Constant faa="<<I[1]<< " : " << faa << endl;
		faa = AF[0]+AF[1];
		cout << "== faa=AF[0]+AF[1] : " << faa << endl;
		faa = faa+2*faa ;
		cout << "== faa=faa+2*faa : " << faa << endl;
		faa = faa*AF[0];
		cout << "== faa = faa*AF[0] : " << faa << endl;


		Variable x(2, "x");
		Variable y("y");
		Function ff(x, y, x[0] * pow(x[1], 2) - exp(x[0] * x[1])+y, "myf");
		cout << "== Using Eval function : "<<ff<< endl;
		Interval fi = ff.eval(I);
		cout << "== Interval results using interval computation of the function: "<< endl << fi << endl;
		Affine2Eval eval_af(ff);
		Affine2Domain dom_faa = eval_af.eval(AF);
		faa = dom_faa.i();
		cout << "== Affine2 results using Affine2 computation of the function: "<< endl << faa << endl;
		Domain dom_itv = eval_af.eval(I);
		fi = dom_itv.i();
		cout << "== Interval results using Affine2 computation of the function: "<< endl << fi << endl;

		cout << "== You can get in the same time Interval and Affine result with one evaluation : " << endl;
		pair<Domain*, Affine2Domain*> res = eval_af.eval(I,AF);
		fi = res.first->i();
		faa = res.second->i();
		cout << fi << endl;
		cout << faa << endl;
		cout <<  endl;


		Function lininf(x, y,
				faa.mid() - faa.err()
				+ faa.val(0) * (2 * x[0] - (I[0].lb() + I[0].ub()))
				/ (I[0].diam())
				+ faa.val(1) * (2 * x[1] - (I[1].lb() + I[1].ub()))
				/ (I[1].diam())
				+ faa.val(2) * (2 * y - (I[2].lb() + I[2].ub()))
				/ (I[2].diam()));
		Function linsup(x, y,
				faa.mid() + faa.err()
				+ faa.val(0) * (2 * x[0] - (I[0].lb() + I[0].ub()))
				/ (I[0].diam())
				+ faa.val(1) * (2 * x[1] - (I[1].lb() + I[1].ub()))
				/ (I[1].diam())
				+ faa.val(2) * (2 * y - (I[2].lb() + I[2].ub()))
				/ (I[2].diam()));


		Function f_inf(x, y, lininf(x, y) - ff(x, y));
		Function f_sup(x, y, ff(x, y) - linsup(x, y));
		NumConstraint c_inf(f_inf, GT);
		NumConstraint c_sup(f_sup, GT);


		SystemFactory sysfac1, sysfac2;
		sysfac1.add_var(x);
		sysfac1.add_var(y);
		sysfac1.add_ctr(c_inf);
		System sys1(sysfac1);

		sysfac2.add_var(x);
		sysfac2.add_var(y);
		sysfac2.add_ctr(c_sup);
		System sys2(sysfac2);


		CtcFwdBwd ct1(f_inf, GT);   // HC4revise Algorithm in the forward step.
		CtcFixPoint ft1(ct1, 0.001);

		CtcFwdBwd ct2(f_sup, GT);  // HC4revise Algorithm in the forward step.
		CtcFixPoint ft2(ct2, 0.001);

		LargestFirst bbb(0.001);
		CellStack ccc;

		// Build a solver
		DefaultSolver sol1(sys1);
		DefaultSolver sol2(sys2);
		Solver::Status stat1 = sol1.solve(I);
		Solver::Status stat2 = sol2.solve(I);

		cout << "==  Validation of the enclosure the Affine2 form : ok?  inf " << stat1 << "  /  sup "<< stat2 << endl;
		if (stat1!=1 || stat2!=1) {
			cout << " size?  :  inf " << sol1.get_nb_cells() << " / sup " << sol2.get_nb_cells() << endl;
			cout << sol1.get_data() << endl;
			cout << "==========================================" << endl;
			cout << sol2.get_data() << endl;
		}

		cout << "==========================================" << endl;
		cout << "==========================================" << endl;
		cout << "== In the same way, you can use Affine3Variables, Affine3 and Affine3Eval. " <<endl;
		cout << "== Affine2 has a static internal structure to improve computation performance. " <<endl;
		cout << "== Affine3 has a dynamic internal structure to use it with the DynIbex project. " <<endl;
		cout << "==========================================" << endl;
		cout << "==========================================" << endl;

	}
	{

		IntervalVector I(3, Interval(1, 2));
		I[1] = Interval(1, 3);

		Affine3Variables AF(I);
		cout << "==  Initialization of the Affine3 Form with 3 epsilons variables: " << endl << AF << endl;

		cout << "==  Affine3Variables is the only way to construct an correct Affine3 Form from an Interval: " << endl ;
		Affine3Variables xx(1 , I[1]);
		cout<< xx[0] << endl;
		Affine3Variables yy(1);
		yy[0] = I[1];
		cout<< yy[0] << endl;

		cout << "==  Constructor Affine3(Interval& itv) and Constructor Affine3Variables(Affine2& af) is not allowed to avoid misunderstandings or bad use. " << endl ;
		cout << "==  Affine3 is used for intermediate result or constant: " << endl ;
		Affine3 faa;
		faa = I[1];
		cout << "== Constant faa="<<I[1]<< " : " << faa << endl;
		faa = AF[0]+AF[1];
		cout << "== faa=AF[0]+AF[1] : " << faa << endl;
		faa = faa+2*faa ;
		cout << "== faa=faa+2*faa : " << faa << endl;
		faa = faa*AF[0];
		cout << "== faa = faa*AF[0] : " << faa << endl;

		Variable x(2, "x");
		Variable y("y");
		Function ff(x, y, x[0] * pow(x[1], 2) - exp(x[0] * x[1])+y, "myf");
		cout << "== Using Eval function : "<<ff<< endl;
		Interval fi = ff.eval(I);
		cout << "== Interval results using interval computation of the function: "<< endl << fi << endl;
		Affine3Eval eval_af(ff);
		Affine3Domain dom_faa = eval_af.eval(AF);
		faa = dom_faa.i();
		cout << "== Affine3 results using Affine3 computation of the function: "<< endl << faa << endl;
		Domain dom_itv = eval_af.eval(I);
		fi = dom_itv.i();
		cout << "== Interval results using Affine3 computation of the function: "<< endl << fi << endl;

		cout << "== You can get in the same time Interval and Affine3 result with one evaluation : " << endl;
		pair<Domain*, Affine3Domain*> res = eval_af.eval(I,AF);
		fi = res.first->i();
		faa = res.second->i();
		cout << fi << endl;
		cout << faa << endl;
		cout <<  endl;

	}

	{
		cout << "==========================================" << endl;
		cout << "==========================================" << endl;
		cout << "TEST 2: " << endl;
		Interval fia, i1, i2;

		i1 = Interval(1, 2);
		i2 = Interval(1, 3);
		fia = i1 * pow(i2, 2) - exp(i1 * i2);

		Affine2Variables aa(2);
		aa[0] = i1;
		aa[1] = i2;
		Affine2 faa = aa[0] * pow(aa[1], 2) - exp(aa[0] * aa[1]);

		cout << aa[0] << endl;
		cout << aa[1] << endl;
		cout << fia << endl;
		cout << faa << endl;

		Affine2Variables ff(1, Interval(0, 1));

		cout << ff[0] << endl;
		cout << pow(ff[0], 2) << endl;

		Affine2 g;
		cout << g << endl;
		g = faa;
		cout << g << endl;
		g = ff[0];
		cout << g << endl;

		/*		cout << " test de l'erreur non modifiable " << endl;
		cout << faa.err() << endl;
		faa.err() += Interval(10000);
		cout << faa.err() << endl;
		 */
		Affine2Variables fff(2);
		fff[0] = Interval(0.5, 1);
		fff[1] = Interval(2, 3);
		cout << fff << endl;

		cout << "test add" << endl;
		cout << fff[0] + fff[1] << endl;
		cout << "test minus" << endl;
		cout << fff[0] - fff[1] << endl;
		cout << "test mul" << endl;
		cout << fff[0] * fff[1] << endl;
		cout << "test div" << endl;
		cout << fff[0] / fff[1] << endl;

		cout << "==========================================" << endl;
		cout << "test log" << endl;
		cout << log(fff[0]) << endl;
		cout << "test inv" << endl;
		cout << 1.0 / (fff[0]) << endl;
		cout << "test exp" << endl;
		cout << exp(fff[0]) << endl;
		cout << "test sqrt" << endl;
		cout << sqrt(fff[0]) << endl;
		cout << "test pow 2" << endl;
		cout << pow(fff[0], 2) << endl;

		cout << "test pow 3" << endl;
		cout << pow(fff[0], 3) << endl;

		cout << "test pow 0" << endl;
		cout << pow(fff[0], 0) << endl;

		cout << "test pow 1" << endl;
		cout << pow(fff[0], 1) << endl;

		cout << "==========================================" << endl;
		cout << "test log" << endl;
		cout << log(fff[1]) << endl;
		cout << "test inv" << endl;
		cout << 1.0 / (fff[1]) << endl;
		cout << "test exp" << endl;
		cout << exp(fff[1]) << endl;
		cout << "test sqrt" << endl;
		cout << sqrt(fff[1]) << endl;
		cout << "test pow 2" << endl;
		cout << pow(fff[1], 2) << endl;

		cout << "test pow 3" << endl;
		cout << pow(fff[1], 3) << endl;

		cout << "test pow 0" << endl;
		cout << pow(fff[1], 0) << endl;

		cout << "test pow 1" << endl;
		cout << pow(fff[1], 1) << endl;
		cout << "==========================================" << endl;

	}


	{
		cout << "==========================================" << endl;
		cout << "==========================================" << endl;
		int n = 1.e4;
		cout << "TEST 3 Performance : " << n << " Evaluations of the Sheckel-5 Function "<< endl;

		double A[5][4] = { { 4, 4, 4, 4 }, { 1, 1, 1, 1 }, { 8, 8, 8, 8 }, { 6,
				6, 6, 6 }, { 3, 7, 3, 7 } };

		double c[5] = { 0.1, 0.2, 0.2, 0.4, 0.4 };
		clock_t start, endtime;
		double cpuTime;
		{
			double f, z;
			double x[4] = { 4, 4, 4, 4 };
			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					for (int j = 0; j < 4; j++) {
						z = z + ::pow((x[j] - A[i][j]), 2);
					}
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / ((double) CLOCKS_PER_SEC);
			cout << "double : CPU-time = " << cpuTime << " seconds" << endl;
			cout << "Output : " << f << endl;
		}
		{
			Interval f, z;
			IntervalVector x(4, Interval(3.9, 4.1));
			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					for (int j = 0; j < 4; j++) {
						z = z + pow((x[j] - A[i][j]), 2);
					}
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Interval : CPU-time = " << cpuTime << " seconds" << endl;
			cout << "Output : " << f << endl;
		}
		{
			Affine2 f, z;
			Affine2Variables x(4, Interval(3.9, 4.1));  // Initialization with x[i] = Affine2(4,i+1,Interval(3.9, 4.1));
			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					for (int j = 0; j < 4; j++) {
						z = z + pow((x[j] - A[i][j]), 2);
					}
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Affine2 : CPU-time = " << cpuTime << " seconds" << endl;
			cout << "Output : " << f << endl;
		}
		{
			Affine3 f, z;
			Affine3Variables x(4, Interval(3.9, 4.1));  // Initialization with x[i] = Affine2(4,i+1,Interval(3.9, 4.1));
			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					for (int j = 0; j < 4; j++) {
						z = z + pow((x[j] - A[i][j]), 2);
					}
					f = f - 1.0 / (z + c[i]);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Affine3 : CPU-time = " << cpuTime << " seconds" << endl;
			cout << "Output : " << f << endl;
		}
		{
			Affine3 f, z;
			Affine3Variables x(4, Interval(3.9, 4.1));  // Initialization with x[i] = Affine2(4,i+1,Interval(3.9, 4.1));
			start = clock();
			for (int k = 0; k < n; k++) {
				f = 0;
				for (int i = 1; i < 5; i++) {
					z = 0;
					for (int j = 0; j < 4; j++) {
						z = z + pow((x[j] - A[i][j]), 2);
					}
					f = f - 1.0 / (z + c[i]);
					f.compact(1.e-7);
				}
			}
			endtime = clock();
			cpuTime = difftime(endtime, start) / CLOCKS_PER_SEC;
			cout << "Affine3 with compact: CPU-time = " << cpuTime << " seconds" << endl;
			cout << "Output : " << f << endl;
		}
//-------------------------------------
// with ibex-lib with GAOL and CLP
//-------------------------------------
// WITH ibex-affine with DEBUG mode
//		TEST 3 Performance : 1.e5 evaluations of the Sheckel-5 Function
//		double : CPU-time = 3.00001e-06 seconds
//		Output : -0.153195
//		Interval : CPU-time = 0.163087 seconds
//		Output : [-0.1663975282993113, -0.1415123357615466]
//		Affine2 : CPU-time = 4.75134 seconds
//		Output : [-0.1611323182907459, -0.1460443472713157] : -0.153588 + -0.000750549 eps_0 + -0.00268127 eps_1 + -0.000750549 eps_2 + -0.00268127 eps_3 + 0.000680337 [-1,1]
//		Affine3 : CPU-time = 105.288 seconds
//		Output : [-0.1614737468029863, -0.1457536702041349] : -0.153613 + -0.000751076 eps_0 + -0.00268241 eps_1 + -0.000751076 eps_2 + -0.00268241 eps_3 + 3.83015e-06 eps_5399996 + 5.78023e-18 eps_5399997 + 7.65888e-06 eps_5399998 + 3.83015e-06 eps_5399999 + 5.78023e-18 eps_5400000 + 7.65888e-06 eps_5400001 + 3.83015e-06 eps_5400002 + 5.78023e-18 eps_5400003 + 7.65888e-06 eps_5400004 + 3.83015e-06 eps_5400005 + 5.78023e-18 eps_5400006 + 7.65888e-06 eps_5400007 + 5.44129e-18 eps_5400008 + -6.40484e-05 eps_5400009 + 1.2157e-06 eps_5400010 + 3.45426e-18 eps_5400011 + 2.43115e-06 eps_5400012 + 1.2157e-06 eps_5400013 + 3.45426e-18 eps_5400014 + 2.43115e-06 eps_5400015 + 1.2157e-06 eps_5400016 + 3.45426e-18 eps_5400017 + 2.43115e-06 eps_5400018 + 1.72731e-18 eps_5400019 + 1.2157e-06 eps_5400020 + 3.45426e-18 eps_5400021 + 2.43115e-06 eps_5400022 + 3.45462e-18 eps_5400023 + 3.45462e-18 eps_5400024 + -2.01278e-05 eps_5400025 + 1.87494e-05 eps_5400026 + 3.74832e-05 eps_5400027 + 1.87494e-05 eps_5400028 + 3.74832e-05 eps_5400029 + 1.87494e-05 eps_5400030 + 3.74832e-05 eps_5400031 + 9.98474e-18 eps_5400032 + 1.87494e-05 eps_5400033 + 3.74832e-05 eps_5400034 + 1.66413e-17 eps_5400035 + -0.00031525 eps_5400036 + 1.21046e-05 eps_5400037 + 2.41687e-05 eps_5400038 + 1.20746e-05 eps_5400039 + 1.82223e-17 eps_5400040 + 2.41447e-05 eps_5400041 + 1.21046e-05 eps_5400042 + 2.41687e-05 eps_5400043 + 9.65316e-18 eps_5400044 + 1.20746e-05 eps_5400045 + 1.82223e-17 eps_5400046 + 2.41447e-05 eps_5400047 + 8.57686e-18 eps_5400048 + -0.000163172 eps_5400049 + 6.07154e-17[-1,1]
//
//-------------------------------------
// WITH ibex-affine with RELEASE mode
//		TEST 3 Performance : 1.e5 evaluations of the Sheckel-5 Function
//		double : CPU-time = 4.00001e-06 seconds
//		Output : -0.153195
//		Interval : CPU-time = 0.164941 seconds
//		Output : [-0.1663975282993113, -0.1415123357615466]
//		Affine2 : CPU-time = 0.892618 seconds
//		Output : [-0.1611323182907459, -0.1460443472713157] : -0.153588 + -0.000750549 eps_0 + -0.00268127 eps_1 + -0.000750549 eps_2 + -0.00268127 eps_3 + 0.000680337 [-1,1]
//		Affine3 : CPU-time = 9.95956 seconds
//		Output : [-0.1614737468029863, -0.1457536702041349] : -0.153613 + -0.000751076 eps_0 + -0.00268241 eps_1 + -0.000751076 eps_2 + -0.00268241 eps_3 + 3.83015e-06 eps_5399996 + 5.78023e-18 eps_5399997 + 7.65888e-06 eps_5399998 + 3.83015e-06 eps_5399999 + 5.78023e-18 eps_5400000 + 7.65888e-06 eps_5400001 + 3.83015e-06 eps_5400002 + 5.78023e-18 eps_5400003 + 7.65888e-06 eps_5400004 + 3.83015e-06 eps_5400005 + 5.78023e-18 eps_5400006 + 7.65888e-06 eps_5400007 + 5.44129e-18 eps_5400008 + -6.40484e-05 eps_5400009 + 1.2157e-06 eps_5400010 + 3.45426e-18 eps_5400011 + 2.43115e-06 eps_5400012 + 1.2157e-06 eps_5400013 + 3.45426e-18 eps_5400014 + 2.43115e-06 eps_5400015 + 1.2157e-06 eps_5400016 + 3.45426e-18 eps_5400017 + 2.43115e-06 eps_5400018 + 1.72731e-18 eps_5400019 + 1.2157e-06 eps_5400020 + 3.45426e-18 eps_5400021 + 2.43115e-06 eps_5400022 + 3.45462e-18 eps_5400023 + 3.45462e-18 eps_5400024 + -2.01278e-05 eps_5400025 + 1.87494e-05 eps_5400026 + 3.74832e-05 eps_5400027 + 1.87494e-05 eps_5400028 + 3.74832e-05 eps_5400029 + 1.87494e-05 eps_5400030 + 3.74832e-05 eps_5400031 + 9.98474e-18 eps_5400032 + 1.87494e-05 eps_5400033 + 3.74832e-05 eps_5400034 + 1.66413e-17 eps_5400035 + -0.00031525 eps_5400036 + 1.21046e-05 eps_5400037 + 2.41687e-05 eps_5400038 + 1.20746e-05 eps_5400039 + 1.82223e-17 eps_5400040 + 2.41447e-05 eps_5400041 + 1.21046e-05 eps_5400042 + 2.41687e-05 eps_5400043 + 9.65316e-18 eps_5400044 + 1.20746e-05 eps_5400045 + 1.82223e-17 eps_5400046 + 2.41447e-05 eps_5400047 + 8.57686e-18 eps_5400048 + -0.000163172 eps_5400049 + 6.07154e-17[-1,1]
//		Affine3 avec compact: CPU-time = 9.41401 seconds
//		Output : [-0.1614737468029855, -0.1457536702041357] : -0.153613 + -0.000751076 eps_0 + -0.00268241 eps_1 + -0.000751076 eps_2 + -0.00268241 eps_3 + -0.00031525 eps_11200034 + 0.00036965 eps_11200035 + 0.000308158 eps_11200049 + 0[-1,1]

//-------------------------------------
// WITH ibex-affine with RELEASE mode
//		TEST 3 Performance : 1.e6 evaluations of the Sheckel-5 Function
//		double : CPU-time = 7.00001e-06 seconds
//		Output : -0.153195
//		Interval : CPU-time = 0.585353 seconds
//		Output : [-0.1663975282993113, -0.1415123357615466]
//		Affine2 : CPU-time = 8.00415 seconds
//		Output : [-0.1611323182907459, -0.1460443472713157] : -0.153588 + -0.000750549 eps_0 + -0.00268127 eps_1 + -0.000750549 eps_2 + -0.00268127 eps_3 + 0.000680337 [-1,1]
//		Affine3 : CPU-time = 99.003 seconds
//		Output : [-0.1614737468029863, -0.1457536702041349] : -0.153613 + -0.000751076 eps_0 + -0.00268241 eps_1 + -0.000751076 eps_2 + -0.00268241 eps_3 + 3.83015e-06 eps_53999996 + 5.78023e-18 eps_53999997 + 7.65888e-06 eps_53999998 + 3.83015e-06 eps_53999999 + 5.78023e-18 eps_54000000 + 7.65888e-06 eps_54000001 + 3.83015e-06 eps_54000002 + 5.78023e-18 eps_54000003 + 7.65888e-06 eps_54000004 + 3.83015e-06 eps_54000005 + 5.78023e-18 eps_54000006 + 7.65888e-06 eps_54000007 + 5.44129e-18 eps_54000008 + -6.40484e-05 eps_54000009 + 1.2157e-06 eps_54000010 + 3.45426e-18 eps_54000011 + 2.43115e-06 eps_54000012 + 1.2157e-06 eps_54000013 + 3.45426e-18 eps_54000014 + 2.43115e-06 eps_54000015 + 1.2157e-06 eps_54000016 + 3.45426e-18 eps_54000017 + 2.43115e-06 eps_54000018 + 1.72731e-18 eps_54000019 + 1.2157e-06 eps_54000020 + 3.45426e-18 eps_54000021 + 2.43115e-06 eps_54000022 + 3.45462e-18 eps_54000023 + 3.45462e-18 eps_54000024 + -2.01278e-05 eps_54000025 + 1.87494e-05 eps_54000026 + 3.74832e-05 eps_54000027 + 1.87494e-05 eps_54000028 + 3.74832e-05 eps_54000029 + 1.87494e-05 eps_54000030 + 3.74832e-05 eps_54000031 + 9.98474e-18 eps_54000032 + 1.87494e-05 eps_54000033 + 3.74832e-05 eps_54000034 + 1.66413e-17 eps_54000035 + -0.00031525 eps_54000036 + 1.21046e-05 eps_54000037 + 2.41687e-05 eps_54000038 + 1.20746e-05 eps_54000039 + 1.82223e-17 eps_54000040 + 2.41447e-05 eps_54000041 + 1.21046e-05 eps_54000042 + 2.41687e-05 eps_54000043 + 9.65316e-18 eps_54000044 + 1.20746e-05 eps_54000045 + 1.82223e-17 eps_54000046 + 2.41447e-05 eps_54000047 + 8.57686e-18 eps_54000048 + -0.000163172 eps_54000049 + 6.07154e-17[-1,1]
//		Affine3 avec compact: CPU-time = 93.3908 seconds
//		Output : [-0.1614737468029855, -0.1457536702041357] : -0.153613 + -0.000751076 eps_0 + -0.00268241 eps_1 + -0.000751076 eps_2 + -0.00268241 eps_3 + -0.00031525 eps_112000034 + 0.00036965 eps_112000035 + 0.000308158 eps_112000049 + 0[-1,1]

	}

	return 0;
}



/* OUTPUT of ./ex_affineform
==========================================
==========================================
TEST 1: Definitions
==  Initialization of the Affine2 Form with 3 epsilons variables:
([1, 2] : 1.5 + 0.5 eps_0 + 0 eps_1 + 0 eps_2 + 0 [-1,1]  ; [1, 3] : 2 + 0 eps_0 + 1 eps_1 + 0 eps_2 + 0 [-1,1]  ; [1, 2] : 1.5 + 0 eps_0 + 0 eps_1 + 0.5 eps_2 + 0 [-1,1] )
==  Affine2Variables is the only way to construct an correct Affine2 Form from an Interval:
[1, 3] : 2 + 1 eps_0 + 0 [-1,1]
[1, 3] : 2 + 1 eps_0 + 0 [-1,1]
==  Constructor Affine2(Interval& itv) and Constructor Affine2Variables(Affine2& af) is not allowed to avoid misunderstandings or bad use.
==  Affine2 is used for intermediate result or constant:
== Constant faa=[1, 3] : [1, 3] : 2 + 1 [-1,1]
== faa=AF[0]+AF[1] : [2, 5] : 3.5 + 0.5 eps_0 + 1 eps_1 + 0 eps_2 + 0 [-1,1]
== faa=faa+2*faa : [6, 15] : 10.5 + 1.5 eps_0 + 3 eps_1 + 0 eps_2 + 0 [-1,1]
== faa = faa*AF[0] : [2.249999999999995, 30.00000000000001] : 16.125 + 7.5 eps_0 + 4.5 eps_1 + 0 eps_2 + 1.87501 [-1,1]
== Using Eval function : myf:(x[2],y)->(((x(1)*x(2)^2)-exp((x(1)*x(2))))+y)
== Interval results using interval computation of the function:
[-401.4287934927352, 17.28171817154096]
== Affine2 results using Affine2 computation of the function:
[-390.4287934927361, 274.6849531133789] : -57.8719 + -77.8921 eps_0 + -114.213 eps_1 + 0.5 eps_2 + 139.952 [-1,1]
== Interval results using Affine2 computation of the function:
[-390.428793492736, 17.28171817154096]
== You can get in the same time Interval and Affine result with one evaluation :
[-390.428793492736, 17.28171817154096]
[-390.4287934927361, 274.6849531133789] : -57.8719 + -77.8921 eps_0 + -114.213 eps_1 + 0.5 eps_2 + 139.952 [-1,1]

==  Validation of the enclosure the Affine2 form : ok?  inf 1  /  sup 1
==========================================
==========================================
== In the same way, you can use Affine3Variables, Affine3 and Affine3Eval.
== Affine2 has a static internal structure to improve computation performance.
== Affine3 has a dynamic internal structure to use it with the DynIbex project.
==========================================
==========================================
==  Initialization of the Affine3 Form with 3 epsilons variables:
([1, 2] : 1.5 + 0.5 eps_0 + 0[-1,1] ; [1, 3] : 2 + 1 eps_1 + 0[-1,1] ; [1, 2] : 1.5 + 0.5 eps_2 + 0[-1,1])
==  Affine3Variables is the only way to construct an correct Affine3 Form from an Interval:
[1, 3] : 2 + 1 eps_0 + 0[-1,1]
[1, 3] : 2 + 1 eps_0 + 0[-1,1]
==  Constructor Affine3(Interval& itv) and Constructor Affine3Variables(Affine2& af) is not allowed to avoid misunderstandings or bad use.
==  Affine3 is used for intermediate result or constant:
== Constant faa=[1, 3] : [1, 3] : 2 + 1 eps_50 + 0[-1,1]
== faa=AF[0]+AF[1] : [2, 5] : 3.5 + 0.5 eps_0 + 1 eps_1 + 0[-1,1]
== faa=faa+2*faa : [6, 15] : 10.5 + 1.5 eps_0 + 3 eps_1 + 0[-1,1]
== faa = faa*AF[0] : [1.5, 30] : 15.75 + 7.5 eps_0 + 4.5 eps_1 + 2.25 eps_51 + 0[-1,1]
== Using Eval function : myf:(x[2],y)->(((x(1)*x(2)^2)-exp((x(1)*x(2))))+y)
== Interval results using interval computation of the function:
[-401.4287934927352, 17.28171817154096]
== Affine3 results using Affine3 computation of the function:
[-392.5569075618581, 276.8379956914853] : -57.8594 + -77.8879 eps_0 + -114.213 eps_1 + 0.5 eps_2 + -40.071 eps_52 + 0.810508 eps_53 + 3.87722e-15 eps_54 + 1.54805 eps_55 + -96.8805 eps_56 + 2.78619 eps_57 + 7.10543e-15 eps_58 + 0[-1,1]
== Interval results using Affine3 computation of the function:
[-392.5569075618581, 17.28171817154096]
== You can get in the same time Interval and Affine3 result with one evaluation :
[-392.5569075618581, 17.28171817154096]
[-392.5569075618581, 276.8379956914853] : -57.8594 + -77.8879 eps_0 + -114.213 eps_1 + 0.5 eps_2 + -40.071 eps_66 + 0.810508 eps_67 + 3.87722e-15 eps_68 + 1.54805 eps_69 + -96.8805 eps_70 + 2.78619 eps_71 + 7.10543e-15 eps_72 + 0[-1,1]

==========================================
==========================================
TEST 2:
[1, 2] : 1.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1]
[1, 3] : 2 + 0 eps_0 + 1 eps_1 + 0 [-1,1]
[-402.4287934927352, 15.28171817154096]
[-391.4287934927366, 216.5149344057413] : -87.4569 + -64.8214 eps_0 + -94.6071 eps_1 + 144.544 [-1,1]
[-0, 1] : 0.5 + 0.5 eps_0 + 0 [-1,1]
[-0.2500000000000005, 1.000000000000001] : 0.375 + 0.5 eps_0 + 0.125001 [-1,1]
[-inf, inf] : Affine Form is not enabled.
[-391.4287934927366, 216.5149344057413] : -87.4569 + -64.8214 eps_0 + -94.6071 eps_1 + 144.544 [-1,1]
[-0, 1] : 0.5 + 0.5 eps_0 + 0 [-1,1]
([0.5, 1] : 0.75 + 0.25 eps_0 + 0 eps_1 + 0 [-1,1]  ; [2, 3] : 2.5 + 0 eps_0 + 0.5 eps_1 + 0 [-1,1] )
test add
[2.5, 4] : 3.25 + 0.25 eps_0 + 0.5 eps_1 + 0 [-1,1]
test minus
[-2.5, -1] : -1.75 + 0.25 eps_0 + -0.5 eps_1 + 0 [-1,1]
test mul
[0.7499999999999996, 3.000000000000001] : 1.875 + 0.625 eps_0 + 0.375 eps_1 + 0.125001 [-1,1]
test div
[0.1123724356957941, 0.5000000000000005] : 0.306187 + 0.102063 eps_0 + -0.0625 eps_1 + 0.0292518 [-1,1]
==========================================
test log
[-0.6931471805599458, 0.05966010114161001] : -0.316743 + 0.346574 eps_0 + 0 eps_1 + 0.0298301 [-1,1]
test inv
[0.8284271247461895, 2.000000000000001] : 1.41422 + -0.5 eps_0 + 0 eps_1 + 0.0857865 [-1,1]
test exp
[1.582104563562883, 2.718281828459048] : 2.1502 + 0.534781 eps_0 + 0 eps_1 + 0.0333084 [-1,1]
test sqrt
[0.7071067811865473, 1.012563132923543] : 0.859835 + 0.146447 eps_0 + 0 eps_1 + 0.00628157 [-1,1]
test pow 2
[0.1874999999999998, 1.000000000000001] : 0.59375 + 0.375 eps_0 + 0 eps_1 + 0.0312501 [-1,1]
test pow 3
[-0.01562500000000024, 1.000000000000001] : 0.492188 + 0.429688 eps_0 + 0 eps_1 + 0.0781251 [-1,1]
test pow 0
<1, 1> : 1 + 0 eps_0 + 0 eps_1 + 0 [-1,1]
test pow 1
[0.5, 1] : 0.75 + 0.25 eps_0 + 0 eps_1 + 0 [-1,1]
==========================================
test log
[0.6931471805599449, 1.119115780042374] : 0.906132 + 0 eps_0 + 0.202733 eps_1 + 0.0102518 [-1,1]
test inv
[0.3164965809277256, 0.5000000000000004] : 0.408249 + 0 eps_0 + -0.0833333 eps_1 + 0.00841838 [-1,1]
test exp
[5.823560187970371, 20.08553692318768] : 12.9546 + 0 eps_0 + 6.34825 eps_1 + 0.782748 [-1,1]
test sqrt
[1.414213562373094, 1.740077828072841] : 1.57715 + 0 eps_0 + 0.158919 eps_1 + 0.00401352 [-1,1]
test pow 2
[3.749999999999999, 9.000000000000002] : 6.375 + 0 eps_0 + 2.5 eps_1 + 0.125001 [-1,1]
test pow 3
[6.124999999999992, 27.00000000000003] : 16.5626 + 0 eps_0 + 9.43751 eps_1 + 1.00001 [-1,1]
test pow 0
<1, 1> : 1 + 0 eps_0 + 0 eps_1 + 0 [-1,1]
test pow 1
[2, 3] : 2.5 + 0 eps_0 + 0.5 eps_1 + 0 [-1,1]
==========================================
==========================================
==========================================
TEST 3 Performance : 10000 Evaluations of the Sheckel-5 Function
double : CPU-time = 3.00001e-06 seconds
Output : -0.153195
Interval : CPU-time = 0.0219041 seconds
Output : [-0.1663975282993113, -0.1415123357615466]
Affine2 : CPU-time = 0.28383 seconds
Output : [-0.1611323182907459, -0.1460443472713157] : -0.153588 + -0.000750549 eps_0 + -0.00268127 eps_1 + -0.000750549 eps_2 + -0.00268127 eps_3 + 0.000680337 [-1,1]
Affine3 : CPU-time = 1.01594 seconds
Output : [-0.1614737468029863, -0.1457536702041349] : -0.153613 + -0.000751076 eps_0 + -0.00268241 eps_1 + -0.000751076 eps_2 + -0.00268241 eps_3 + 3.83015e-06 eps_540019 + 5.78023e-18 eps_540020 + 7.65888e-06 eps_540021 + 3.83015e-06 eps_540022 + 5.78023e-18 eps_540023 + 7.65888e-06 eps_540024 + 3.83015e-06 eps_540025 + 5.78023e-18 eps_540026 + 7.65888e-06 eps_540027 + 3.83015e-06 eps_540028 + 5.78023e-18 eps_540029 + 7.65888e-06 eps_540030 + 5.44129e-18 eps_540031 + -6.40484e-05 eps_540032 + 1.2157e-06 eps_540033 + 3.45426e-18 eps_540034 + 2.43115e-06 eps_540035 + 1.2157e-06 eps_540036 + 3.45426e-18 eps_540037 + 2.43115e-06 eps_540038 + 1.2157e-06 eps_540039 + 3.45426e-18 eps_540040 + 2.43115e-06 eps_540041 + 1.72731e-18 eps_540042 + 1.2157e-06 eps_540043 + 3.45426e-18 eps_540044 + 2.43115e-06 eps_540045 + 3.45462e-18 eps_540046 + 3.45462e-18 eps_540047 + -2.01278e-05 eps_540048 + 1.87494e-05 eps_540049 + 3.74832e-05 eps_540050 + 1.87494e-05 eps_540051 + 3.74832e-05 eps_540052 + 1.87494e-05 eps_540053 + 3.74832e-05 eps_540054 + 9.98474e-18 eps_540055 + 1.87494e-05 eps_540056 + 3.74832e-05 eps_540057 + 1.66413e-17 eps_540058 + -0.00031525 eps_540059 + 1.21046e-05 eps_540060 + 2.41687e-05 eps_540061 + 1.20746e-05 eps_540062 + 1.82223e-17 eps_540063 + 2.41447e-05 eps_540064 + 1.21046e-05 eps_540065 + 2.41687e-05 eps_540066 + 9.65316e-18 eps_540067 + 1.20746e-05 eps_540068 + 1.82223e-17 eps_540069 + 2.41447e-05 eps_540070 + 8.57686e-18 eps_540071 + -0.000163172 eps_540072 + 6.07154e-17[-1,1]
Affine3 with compact: CPU-time = 0.941321 seconds
Output : [-0.1614737468029855, -0.1457536702041357] : -0.153613 + -0.000751076 eps_0 + -0.00268241 eps_1 + -0.000751076 eps_2 + -0.00268241 eps_3 + -0.00031525 eps_1120057 + 0.00036965 eps_1120058 + 0.000308158 eps_1120072 + 0[-1,1]


*/









