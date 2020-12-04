//============================================================================
//                                  I B E X                                   
// File        : ibex_Optimizer04Config_Affine.h
// Author      : Bertrand Neveu, Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Dec 11, 2014
// Last Update : Dec 04, 2020
//============================================================================

#ifndef __IBEX_OPTIMIZER_04_CONFIG_AFFINE__
#define __IBEX_OPTIMIZER_04_CONFIG_AFFINE__

#include "ibex_Optimizer04Config.h"

namespace ibex {

class Optimizer04Config_Affine : public Optimizer04Config {
public:

	Optimizer04Config_Affine(int argc, char** argv);

protected:

	virtual Linearizer* get_linear_relax();
};

} // namespace ibex

#endif // __IBEX_OPTIMIZER_04_CONFIG_AFFINE__
