//============================================================================
//                                  I B E X                                   
// File        : Ponts30.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Apr 10, 2012
// Last Update : Apr 10, 2012
//============================================================================

#include "Ponts30.h"

namespace ibex {


void Ponts30::build_box() {

	double _init_box[30][2]= {
			{ -46.516999999999996, 53.483000000000004 } ,
			{ -53.574999999999996, 46.425000000000004 } ,
			{ -46.759999999999998, 53.240000000000002 } ,
			{ -53.517999999999994, 46.482000000000006 } ,
			{ -46.589999999999996, 53.410000000000004 } ,
			{ -53.335999999999999, 46.664000000000001 } ,
			{ -46.602999999999994, 53.397000000000006 } ,
			{ -53.085999999999999, 46.914000000000001 } ,
			{ -46.832999999999998, 53.167000000000002 } ,
			{ -53.278999999999996, 46.721000000000004 } ,
			{ -47.003999999999998, 52.996000000000002 } ,
			{ -53.461999999999996, 46.538000000000004 } ,
			{ -49.599999999999994, 50.400000000000006 } ,
			{ -51.958999999999996, 48.041000000000004 } ,
			{ -47.596999999999994, 52.403000000000006 } ,
			{ -48.497999999999998, 51.502000000000002 } ,
			{ -46.552999999999997, 53.447000000000003 } ,
			{ -52.840999999999994, 47.159000000000006 } ,
			{ -46.789999999999999, 53.210000000000001 } ,
			{ -52.919999999999995, 47.080000000000005 } ,
			{ -46.976999999999997, 53.023000000000003 } ,
			{ -52.753999999999998, 47.246000000000002 } ,
			{ -46.739999999999995, 53.260000000000005 } ,
			{ -52.674999999999997, 47.325000000000003 } ,
			{ -46.503, 53.497 } ,
			{ -52.595999999999997, 47.404000000000003 } ,
			{ 0, 55 } ,
			{ -50, 50 } ,
			{ -50, 50 } ,
			{ -50, 50 } };

	double bounds[30][30][2] = {
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-49.6,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-50,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-47.597,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-46.503,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-7.243,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-3,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-7,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-3,53.497} , {-52.596,47.404} , {0,55} , {-50,50} , {-5,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-3,53.497} , {-52.596,47.404} , {0,55} , {0,0} , {-5,50} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-3,53.497} , {-52.596,47.404} , {0,55} , {0,0} , {0,0} , {-50,50}},
	{{-46.517,53.483} , {-53.575,46.425} , {-46.76,53.24} , {-53.518,46.482} , {-46.59,53.41} , {-53.336,46.664} , {-46.603,53.397} , {-53.086,46.914} , {-46.833,53.167} , {-53.279,46.721} , {-47.004,52.996} , {-53.462,46.538} , {-5,50.4} , {-51.959,48.041} , {-3,52.403} , {-48.498,51.502} , {-46.553,53.447} , {-52.841,47.159} , {-46.79,53.21} , {-52.92,47.08} , {-46.977,53.023} , {-52.754,47.246} , {-46.74,53.26} , {-52.675,47.325} , {-3,53.497} , {-52.596,47.404} , {0,55} , {0,0} , {0,0} , {0,0}} };


	double _hc4_box[30][2]= {
			{1,3.91742} ,
			{-4,4} ,
			{1,3.66742} ,
			{-4,4} ,
			{1.25,3.91742} ,
			{-3.75,3.75} ,
			{1.5,3.96742} ,
			{-3.5,3.5} ,
			{1.2,3.66742} ,
			{-3.8,3.8} ,
			{0.95,3.41742} ,
			{-4.05,4.05} ,
			{-4.44089e-15,0.417424} ,
			{-2,2} ,
			{2,4.41742} ,
			{-2.94289,2.94289} ,
			{1.75,4.21742} ,
			{-3.25,3.25} ,
			{1.5,4.21742} ,
			{-3.5,3.5} ,
			{1.5,4.46742} ,
			{-3.5,3.5} ,
			{1.75,4.46742} ,
			{-3.25,3.25} ,
			{2,4.46742} ,
			{-2.95235,2.95235} ,
			{5,5} ,
			{-0,0} ,
			{0,0} ,
			{-0,0} };

	init_box=IntervalVector(30,_init_box);

	for (int i=0; i<30; i++) {
		hc4r_box[i]=IntervalVector(30,bounds[i]);
	}

	hc4_box=IntervalVector(30,_hc4_box);
}

void Ponts30::build_equ() {

	const ExprSymbol& O_y=ExprSymbol::new_("O_y");
	const ExprSymbol& O_x=ExprSymbol::new_("O_x");
	const ExprSymbol& N_y=ExprSymbol::new_("N_y");
	const ExprSymbol& N_x=ExprSymbol::new_("N_x");
	const ExprSymbol& M_y=ExprSymbol::new_("M_y");
	const ExprSymbol& M_x=ExprSymbol::new_("M_x");
	const ExprSymbol& L_y=ExprSymbol::new_("L_y");
	const ExprSymbol& L_x=ExprSymbol::new_("L_x");
	const ExprSymbol& K_y=ExprSymbol::new_("K_y");
	const ExprSymbol& K_x=ExprSymbol::new_("K_x");
	const ExprSymbol& J_y=ExprSymbol::new_("J_y");
	const ExprSymbol& J_x=ExprSymbol::new_("J_x");
	const ExprSymbol& I_y=ExprSymbol::new_("I_y");
	const ExprSymbol& I_x=ExprSymbol::new_("I_x");
	const ExprSymbol& H_y=ExprSymbol::new_("H_y");
	const ExprSymbol& H_x=ExprSymbol::new_("H_x");
	const ExprSymbol& G_y=ExprSymbol::new_("G_y");
	const ExprSymbol& G_x=ExprSymbol::new_("G_x");
	const ExprSymbol& F_y=ExprSymbol::new_("F_y");
	const ExprSymbol& F_x=ExprSymbol::new_("F_x");
	const ExprSymbol& E_y=ExprSymbol::new_("E_y");
	const ExprSymbol& E_x=ExprSymbol::new_("E_x");
	const ExprSymbol& D_y=ExprSymbol::new_("D_y");
	const ExprSymbol& D_x=ExprSymbol::new_("D_x");
	const ExprSymbol& C_y=ExprSymbol::new_("C_y");
	const ExprSymbol& C_x=ExprSymbol::new_("C_x");
	const ExprSymbol& B_y=ExprSymbol::new_("B_y");
	const ExprSymbol& B_x=ExprSymbol::new_("B_x");
	const ExprSymbol& A_y=ExprSymbol::new_("A_y");
	const ExprSymbol& A_x=ExprSymbol::new_("A_x");

	Array<const ExprSymbol> symbols(30);
	int i=0;
	symbols.set_ref(i++,O_y);
	symbols.set_ref(i++,O_x);
	symbols.set_ref(i++,N_y);
	symbols.set_ref(i++,N_x);
	symbols.set_ref(i++,M_y);
	symbols.set_ref(i++,M_x);
	symbols.set_ref(i++,L_y);
	symbols.set_ref(i++,L_x);
	symbols.set_ref(i++,K_y);
	symbols.set_ref(i++,K_x);
	symbols.set_ref(i++,J_y);
	symbols.set_ref(i++,J_x);
	symbols.set_ref(i++,I_y);
	symbols.set_ref(i++,I_x);
	symbols.set_ref(i++,H_y);
	symbols.set_ref(i++,H_x);
	symbols.set_ref(i++,G_y);
	symbols.set_ref(i++,G_x);
	symbols.set_ref(i++,F_y);
	symbols.set_ref(i++,F_x);
	symbols.set_ref(i++,E_y);
	symbols.set_ref(i++,E_x);
	symbols.set_ref(i++,D_y);
	symbols.set_ref(i++,D_x);
	symbols.set_ref(i++,C_y);
	symbols.set_ref(i++,C_x);
	symbols.set_ref(i++,B_y);
	symbols.set_ref(i++,B_x);
	symbols.set_ref(i++,A_y);
	symbols.set_ref(i++,A_x);

	// (reverse order)
//	const ExprSymbol& A_y=ExprSymbol::new_("A_y");
//	const ExprSymbol& B_x=ExprSymbol::new_("B_x");
//	const ExprSymbol& B_y=ExprSymbol::new_("B_y");
//	const ExprSymbol& C_x=ExprSymbol::new_("C_x");
//	const ExprSymbol& C_y=ExprSymbol::new_("C_y");
//	const ExprSymbol& D_x=ExprSymbol::new_("D_x");
//	const ExprSymbol& D_y=ExprSymbol::new_("D_y");
//	const ExprSymbol& E_x=ExprSymbol::new_("E_x");
//	const ExprSymbol& E_y=ExprSymbol::new_("E_y");
//	const ExprSymbol& F_x=ExprSymbol::new_("F_x");
//	const ExprSymbol& F_y=ExprSymbol::new_("F_y");
//	const ExprSymbol& G_x=ExprSymbol::new_("G_x");
//	const ExprSymbol& G_y=ExprSymbol::new_("G_y");
//	const ExprSymbol& H_x=ExprSymbol::new_("H_x");
//	const ExprSymbol& H_y=ExprSymbol::new_("H_y");
//	const ExprSymbol& I_x=ExprSymbol::new_("I_x");
//	const ExprSymbol& I_y=ExprSymbol::new_("I_y");
//	const ExprSymbol& J_x=ExprSymbol::new_("J_x");
//	const ExprSymbol& J_y=ExprSymbol::new_("J_y");
//	const ExprSymbol& K_x=ExprSymbol::new_("K_x");
//	const ExprSymbol& K_y=ExprSymbol::new_("K_y");
//	const ExprSymbol& L_x=ExprSymbol::new_("L_x");
//	const ExprSymbol& L_y=ExprSymbol::new_("L_y");
//	const ExprSymbol& M_x=ExprSymbol::new_("M_x");
//	const ExprSymbol& M_y=ExprSymbol::new_("M_y");
//	const ExprSymbol& N_x=ExprSymbol::new_("N_x");
//	const ExprSymbol& N_y=ExprSymbol::new_("N_y");
//	const ExprSymbol& O_x=ExprSymbol::new_("O_x");
//	const ExprSymbol& O_y=ExprSymbol::new_("O_y");

	Array<const ExprNode> equ(30);
	i=0;

	equ.set_ref(i++,(	sqr(N_x - O_x) + sqr(N_y - O_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(M_x - O_x) + sqr(M_y - O_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(M_x - N_x) + sqr(M_y - N_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(J_x - N_x) + sqr(J_y - N_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(L_x - M_x) + sqr(L_y - M_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(K_x - M_x) + sqr(K_y - M_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(G_x - L_x) + sqr(G_y - L_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(K_x - L_x) + sqr(K_y - L_y) - 0.089999999999997 ));
	equ.set_ref(i++,(	sqr(J_x - K_x) + sqr(J_y - K_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(K_x - N_x) + sqr(K_y - N_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(I_x - J_x) + sqr(I_y - J_y) - 9 ));
	equ.set_ref(i++,(	sqr(H_x - J_x) + sqr(H_y - J_y) - 25 ));
	equ.set_ref(i++,(	sqr(B_x - I_x) + sqr(B_y - I_y) - 25 ));
	equ.set_ref(i++,(	sqr(A_x - I_x) + sqr(A_y - I_y) - 4 ));
	equ.set_ref(i++,(	sqr(B_x - H_x) + sqr(B_y - H_y) - 9 ));
	equ.set_ref(i++,(	sqr(H_x - I_x) + sqr(H_y - I_y) - 16 ));
	equ.set_ref(i++,(	sqr(F_x - G_x) + sqr(F_y - G_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(C_x - G_x) + sqr(C_y - G_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(F_x - L_x) + sqr(F_y - L_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(D_x - F_x) + sqr(D_y - F_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(D_x - E_x) + sqr(D_y - E_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(E_x - F_x) + sqr(E_y - F_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(D_x - G_x) + sqr(D_y - G_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(C_x - D_x) + sqr(C_y - D_y) - 0.0625 ));
	equ.set_ref(i++,(	sqr(C_x - H_x) + sqr(C_y - H_y) - 18.003049000000004 ));
	equ.set_ref(i++,(	sqr(B_x - C_x) + sqr(B_y - C_y) - 9 ));
	equ.set_ref(i++,(	sqr(A_x - B_x) + sqr(A_y - B_y) - 25 ));
	equ.set_ref(i++,(	B_x ));
	equ.set_ref(i++,(	A_y ));
	equ.set_ref(i++,(	A_x ));

	f =  new Function(symbols,ExprVector::new_col(equ).simplify(ExprNode::default_simpl_level),"ponts30");
}

Ponts30::Ponts30() : init_box(30), hc4r_box(30,30), hc4_box(30) {
	build_equ();
	build_box();
}


} // end namespace ibex
