//============================================================================
//                                  I B E X
// File        : optimizer04.cpp
// Author      : Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 11, 2014
// Last Update : Oct 15, 2019
//============================================================================

#include "ibex_Optimizer04Config_Affine.h"
#include "ibex_Optimizer.h"

using namespace ibex;

int main(int argc, char** argv) {

	Optimizer04Config_Affine config(argc,argv);

	Optimizer o(config);

	o.optimize(config.sys->box);

	o.report();

}
