//General FEFLOW tools module
#pragma once

/* Infos on the module:
Purpose: include in a FEFLOW IFM plugin code written in C++, to simplify the call to some API functions
Author and copyright: Marc Laurencelle
Last updates: Oct. 2017
*/

#include "stdifm.h"

//BETA!
bool genFFtools_extract_nodal_param(IfmDocument pDoc, IfmParamID nId, double* destpValarray)
{
	bool good;
	int totNofNodes = IfmGetNumberOfNodes(pDoc); //number of nodes (global)
	int parsiz = IfmGetParamSize(pDoc, nId);
	if (parsiz != totNofNodes) IfmWarning(pDoc, "[genFFtools_extract_nodal_param] STRANGELY, param.size != nb.of.Nodes !!?");
	int Npvread = IfmGetParamValues(pDoc, nId, destpValarray, 0, totNofNodes); //number of parameter values read successfully
	good = Npvread == totNofNodes;
	if (!good) IfmWarning(pDoc, "[genFFtools_extract_nodal_param] STRANGELY, nb.param.values.read != nb.of.Nodes !!?");
	//IF DONE with no warning, should be OK then, and ready for use by other analysis/post-processing functions.
	return good;
}
