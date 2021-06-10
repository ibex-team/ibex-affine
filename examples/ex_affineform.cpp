//============================================================================
//                                  I B E X
// File        : arith03.cpp
// Author      : Jordan Ninin
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Fev 28, 2013
// Last Update : Fev 28, 2013
//============================================================================

#include "ibex.h"
#include "ibex_Affine.h"
#include <time.h>

using namespace std;
using namespace ibex;

int main() {

	{
		cout << "==========================================" << endl;
		cout << "==========================================" << endl;
		cout << "TEST 1: " << endl;
		Variable x(2);
		Variable x1, x2;

		Function ff(x, x2, x[0] * pow(x[1], 2) - exp(x[0] * x[1]));

		IntervalVector I(3, Interval(1, 2));
		I[1] = Interval(1, 3);
		Affine2Variables AF(I);
		Interval fi = ff.eval(I);
		cout << fi << endl;
		Affine2Eval eval_af(ff);
		Affine2Domain dom_faa = eval_af.eval(AF);
		Affine2 faa = dom_faa.i();
		/* OR this way, it is the same
		std::pair<Domain*, Affine2Domain*> res = eval_af.eval(I,AF);
		Affine2 faa = res.second->i();
		 */
		cout << faa << endl;

		//		Function lininf(x, faa.val(0)-faa.err().ub() + faa.val(1)*(2*x[0]-(I[0].lb()+I[0].ub()))/(I[0].diam()) + faa.val(2)*(2*x[1]-(I[1].lb()+I[1].ub()))/(I[1].diam())) ;
		//		Function linsup(x, faa.val(0)+faa.err().ub() + faa.val(1)*(2*x[0]-(I[0].lb()+I[0].ub()))/(I[0].diam()) + faa.val(2)*(2*x[1]-(I[1].lb()+I[1].ub()))/(I[1].diam())) ;

		Function lininf(x, x2,
				faa.mid() - faa.err()
				+ faa.val(0) * (2 * x[0] - (I[0].lb() + I[0].ub()))
				/ (I[0].diam())
				+ faa.val(1) * (2 * x[1] - (I[1].lb() + I[1].ub()))
				/ (I[1].diam())
				+ faa.val(2) * (2 * x2 - (I[2].lb() + I[2].ub()))
				/ (I[2].diam()));
		Function linsup(x, x2,
				faa.mid() + faa.err()
				+ faa.val(0) * (2 * x[0] - (I[0].lb() + I[0].ub()))
				/ (I[0].diam())
				+ faa.val(1) * (2 * x[1] - (I[1].lb() + I[1].ub()))
				/ (I[1].diam())
				+ faa.val(2) * (2 * x2 - (I[2].lb() + I[2].ub()))
				/ (I[2].diam()));

		cout << lininf << endl;

		Function f_inf(x, x2, lininf(x, x2) - ff(x, x2));
		Function f_sup(x, x2, ff(x, x2) - linsup(x, x2));
		NumConstraint c_inf(f_inf, GT);
		NumConstraint c_sup(f_sup, GT);


		SystemFactory sysfac1, sysfac2;
		sysfac1.add_var(x);
		sysfac1.add_var(x2);
		sysfac1.add_ctr(c_inf);
		System sys1(sysfac1);

		sysfac2.add_var(x);
		sysfac2.add_var(x2);
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

		cout << " ok?  : inf " << stat1 << "  /  sup "<< stat2 << endl;
		if (stat1!=1 || stat2!=1) {
			cout << " size?  :  inf " << sol1.get_nb_cells() << " / sup " << sol2.get_nb_cells() << endl;
			cout << sol1.get_data() << endl;
			cout << "==========================================" << endl;
			cout << sol2.get_data() << endl;
		}
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
		int n = 1.e6;
		cout << "TEST 3 Performance : " << n << " evaluations of the Sheckel-5 Function "<< endl;
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
			cout << "Affine3 avec compact: CPU-time = " << cpuTime << " seconds" << endl;
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

