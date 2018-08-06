/* Infos on the module:
Purpose: include in a FEFLOW IFM plugin code written in C++, to export results of a transient run
Author and copyright: Marc Laurencelle
Last updates: Oct. 2017
*/

/* TODO list:
- Replace the NC_SHARE mode (which may prevent perf. optim. due to frequent write ops to disk) by autom. calling of
  nc_sync(...) at regular time intervals. Store time vars in the nc object. With t=now when initializing nc.
  [not a serious optimization, except in cases with high-freq. export...; therefore, would be nice for the General uses]

- Consider using CZ_dynA_Eindeces for the usual CZsums computations, to avoid looping over all elements. I may even
  remove from memory the CZid arrays once CZ_dynA_Eindeces is built, then!...

- Consider using the vector-specialized library includes for simpler calculations maybe: <functional> <numeric>. However, it may be harder to read the code, since it seems it requires a transform call in order to work... (e.g. look http://www.cplusplus.com/reference/functional/plus/)

DONE | ABANDONED list:
x ABANDONED IDEA: Change the order of dimensions for CZ_quantiles so that it's easier to use the data in R: currently it's «conc_quantiles : num [1:101, 1:4, 1:199]». [CZ#, Q#, T#] would be better I think.
*/


#include "stdifm.h"
#include "Pathcch.h"
#include "MarcGeneralFEFLOWtools.h"

#include <stdio.h>
#include <string.h>
#include <netcdf.h>

#include <vector>
#include <algorithm>
#include <numeric>

#define NDIMS2 2
#define VECTCNDIMS3 3 //for vectors' components (Darcy fluxes...); here, 3 is for the export-array dimensions[time,vect.comp.,node]
#define CZNDIMS2 2
#define CZQNDIMS3 3
#define MZNDIMS2 2
#define NBNEGCNDIMS2 2 //related to NbCuppBnds...

#define VECTC_FIXED_LEN 2 //2D models ONLY, for now (~NUMBER OF SPATIAL DIMENSIONS IN THE MODEL)

#define NOD_DIM_NAME "nodei"
#define REC_DIM_NAME "time"
#define VECTC_DIM_NAME "vectcompi"
#define CZ_DIM_NAME "contzonei"
#define QPROB_DIM_NAME "qprob"
#define MZ_DIM_NAME "monitzonei"
#define NEGCUBND_DIM_NAME "negconc_ubnd"

/* Names of things. */
#define UNITS "units"
#define NODEFF_NAME "nodenumber"
#define NODIDV_NAME "nodeidvalue"
#define XCOORD_NAME "x_local"
#define YCOORD_NAME "y_local"
#define YOFTOPROCK_NAME "yoftoprock_local"
#define HEAD_NAME "hydraulic_head"
#define CONC_NAME "concentration"
#define DFLUX_NAME "Darcy_flux_nodal"
#define MRBUDG_NAME "mass_rate_budget"
#define CQUANT_NAME "conc_quantiles"
#define VCOMD_NAME "vCOM_depth"
#define DPFD_NAME "DPF_depth"
#define GLOBMINCONC_NAME "glob_min_concentration"
#define GLOBMAXCONC_NAME "glob_max_concentration"
#define VCONCGRAD_NAME "vertical_conc_gradient"
#define NBNEGC_NAME "nb_negative_concs"

/* For the units attributes. */
#define UNITS "units"
//a. for Dimensions:
#define NODE_UNITS "node_as_ordered"
#define TIME_UNITS "days_since_sim_start"
#define VECTC_UNITS "vect_component_index"
#define CZ_UNITS "content_zone_raw_index"
#define MZ_UNITS "monit_zone_raw_index"
#define QPROB_UNITS "probability_at_quantile"
#define NEGCUBND_UNITS "negative_conc_upper_bound"
//b. for Variables:
#define NODEFF_UNITS "index_in_FF_mesh"
#define NODIDV_UNITS "id_value_by_user"
#define XCOORD_UNITS "m_rel_to_local_orig"
#define YCOORD_UNITS "m_rel_to_local_orig"
#define YOFTOPROCK_UNITS "m_rel_to_local_orig"
#define HEAD_UNITS "meter"
#define CONC_UNITS "mg_per_liter"
#define DFLUX_UNITS "meter_per_day"
#define MRBUDG_UNITS "g_per_day"
#define CQUANT_UNITS "quantile_at_Prob"
#define VCOMD_UNITS "meters below ref. top (vertical center of solute mass)"
#define DPFD_UNITS "meters below ref. top (depth of deepest plume front)"
#define VCONCGRAD_UNITS "mg_per_liter_per_meter"
#define NBNEGC_UNITS "number of nodal conc. values below the upper bounds"

//#define MAX_ATT_LEN 80

/* Handle errors by printing an error message and exiting with a non-zero status. */
#define ERR(e) {sprintf_s(nh.errtxtbuffer, 180, "[NetCDF] Error: %s", nc_strerror(e)); IfmWarning(pDoc, nh.errtxtbuffer); return 2;}
#define ERRpropag(e) {return e;}

const int MAXCTYPE = IfmENERGY_TOTAL; //DO NOT use this constant directly; instead, please use NBCTYPES!
const int NBCTYPES = MAXCTYPE + 1;

const int HWnbProbDivs = 100; //HARD-WIRED, for now... (TODO minor: Make it customizable from the FEM?)

const double nodata_double_ncdf = std::numeric_limits<double>::quiet_NaN();

#define NbNegCuppBnds 6
const double negCubounds[NbNegCuppBnds] = { nodata_double_ncdf, 0.0, -1000.0, -5000.0, -10000.0, -50000.0 };

//based on http://stackoverflow.com/questions/2185443/does-c-support-constant-arrays-of-type-string
const char* CTYPE_namesA[NBCTYPES] = { "total_volume", "void_volume", "fluid_content", "diluted_mass", "sorbed_mass", "energy_fluid", "energy_solid", "total_energy" };
const char* CTYPE_unitsA[NBCTYPES] = { "m3", "m3", "m3", "g", "g", "J", "J", "J" }; //based on help documentation on 'IfmGetElementalContent'

struct ncdfheader {
	//OPTIONS FOR CONTROLLING WHAT TO EXPORT:
	bool putNodeFF; //node index in the FEFLOW model mesh (beware: it changes very often when mesh geometry is modified)
	bool putX, putY;
	bool putYoftoprock; //NEW; when true, requires a fully set-up Yoftoprock nodal user data
	bool puthead, putconc;
	bool putDfluxcomponents;
	bool putmratebudg;
	bool putconcquantiles; //flag for activation of quantiles computing & export, throughout the simulation; ignored IF nbcontentzones == 0
	bool putmassconvectmonit; //flag for potential activation of MZ recording, during the _init_ stage ONLY; later, we rather check if nbmonitzones > 0 to know if it's effectively active!
	bool putnbnegconcs; //NEW flag... BETA
	bool putglobextrconcs; //NEW flag... BETA
	bool putvconcgrads; //NEW flag... BETA
	
	//Some parameters which can be read from outside this module (but which are set internally ONLY!)
	int nbcontentzones = 0;
	int nbmonitzones = 0;
	bool arefixednodes = false; //NEW!

	//Special arrays pointers to customize some variables
	double *PtocompmrbudgetA = nullptr;

	//Internal copies of nodal data for all nodes in their original FEFLOW ordering
	// (all have a fixed length = the numb 
	// ('fixed' ones are filled with steady values only once, during the netCDF initialization)
	double *copyoffixednodalXcoordsA = nullptr;
	double *copyoffixednodalYcoordsA = nullptr;
	// ('transient' ones are updated with new values during each export-time-step write call)
	double *copyoftransientnodalconcsA = nullptr;

	//Arrays for computation of some types of export data
	int compnbnegconcsA[NbNegCuppBnds]; //array for computation of the number of negative concentrations below each of the conc. upper bounds hard-wired in this code; (fixed-length static array, hence: no need to deallocate it at the end)
	double *compvconcgradsA = nullptr; //array for computation of vertical concentration gradients at the export nodes

	//Data arrays for the export nodes (expN)
	int *Pnodeff_expN_outA; //array of the 'effective' i.e. current node FF-indeces for the export nodes (values change occasionally if not .arefixednodes)
	double *Pxcoord_expN_outA; //array of the current X coordinates for the export nodes (values change occasionally if not .arefixednodes)
	double *Pycoord_expN_outA; //array of the current Y coordinates for the export nodes (values change occasionally if not .arefixednodes)
	double *Pyoftoprock_expN_outA; //array of the current Y-of-top-rock user-data-specified coordinates for the export nodes (values are updated occasionally if not .arefixednodes... though this normally should not change the values)

	long opt_Yoftoprock_rID = -1; //Optional nodal reference distribution ID for the 'Yoftoprock' r.d. (required only if putYoftoprock == true)

	//Special dyn.-sized array of vectors storing elem. indeces for faster CZ quantile, or MZ indicators, computations
	//(NOTE: They are sized only once, at netcdf initialization, and keep a constant size afterwards.)
	std::vector<int> *CZ_dynA_Eindeces = nullptr; //one vector per CZ id
	std::vector<int> *MZ_dynA_Eindeces = nullptr; //one vector per MZ id

	std::vector<double> CZ_fixed_qprobsV; //fixed probabilities P(q) for quantile computations
	double *CZ_2DquantsA; //row i = CZ index; column j = quantile index
	double *DFLUX_2DvaluesA; //row i = vector-component index; column j = node-as-ordered index

	double *MZ_vCOMdA; //internal storage of most recent computed vCOM_depth values: array[nbMZ]
	double *MZ_DPFdA; //internal storage of most recent computed DPF_depth values: array[nbMZ]

	double comp_globminc; //computed global minimum concentration at current export time step; from scan of all active nodes (i.e. non-NaN conc. values)
	double comp_globmaxc; //computed global maximum concentration at current export time step; from scan of all active nodes (i.e. non-NaN conc. values)

	//FILE ID
	char filepath[_MAX_PATH];
	int ncid; //ID for the netCDF file

	//For handling error messages
	char errtxtbuffer[180];

	//DIMENSIONS IDs
	int nod_dimid; //ID for the nodal dimension (not the same as FEFLOW node index)
	int rectime_dimid; //ID for the recorded time step (not the same as FEFLOW time step index)
	int vectc_dimid; //ID for the vectors' components dimension (Darcy flux export...)
	int dimids[NDIMS2];
	int Dflux_dimids[VECTCNDIMS3];
	//variants for CZ:
	int CZ_dimid;
	int qprob_dimid;
	int CZdimids[CZNDIMS2];
	int CZQdimids[CZQNDIMS3];
	//variants for MZ:
	int MZ_dimid;
	int MZdimids[MZNDIMS2];
	//variants for NEGC:
	int negcubnd_dimid;
	int nbnegcdimids[NBNEGCNDIMS2];

	//VARIABLE IDs
	//a. for Dimensions
	int rectime_varid;
	int nodindex_varid;
	int vectc_varid;
	int CZindex_varid;
	int qprob_varid;
	int MZindex_varid;
	int negcubnd_varid;
	//b. for Variables (dep. on the put___ selection)
	int nodeff_varid;
	int nodidv_varid; //NEW: id value by user
	int xcoord_varid;
	int ycoord_varid;
	int yoftoprock_varid;
	int head_varid;
	int conc_varid;
	int Dflux_varid;
	int mrbudg_varid;
	int cquant_varid;
	int *CZ_varidA; //NEW
	int vCOMd_varid;
	int DPFd_varid;
	int vconcgrad_varid; //NEW
	int globminconc_varid; //NEW
	int globmaxconc_varid; //NEW
	int nbnegc_varid; //NEW

	/* Special external-assigned parameters to make it possible to compute DPFd (if putmassconvectmonit...) */
	// (these parameters MUST be set prior to the simulation run, if .putmassconvectmonit & then nbmonitzones > 0)
	double DPF_normconcmin = 0.0; //minimum normalized concentration to look for, to locate the DPF...
	double ref_C0 = -9999.0; //reference conc. for minimum density (in mg/L): rho(C0) = rho0 + 0.0*densityratio
	double ref_CS = -9999.0; //reference conc. for maximum density (in mg/L): rho(CS) = rho0 + 1.0*densityratio
	
	int nbelems_fixed = 0; //used by several tools of this module; must be defined prior to calling the netcdftools_init_... functions; TODO-minor: Could add a protection which verifies that nbelems_fixed>0 before using it
	int nbnodes_fixed = 0; //used by several tools of this module; must be defined prior to calling the netcdftools_init_... functions; TODO-minor: Could add a protection which verifies that nbnodes_fixed>0 before using it
	int *contentzoneidsA = nullptr; //optional; set from within the main IFM plugin code, if the related elemental User Data is found; otherwise, it remains null while nbcontentzones = 0 still.
	int *uniqsortedczidsA = nullptr;
	double **CZcomputedsums; //BETA: first dimension is based on the IfmContentType enum type, which current last value in FEFLOW IFM document.h is IfmENERGY_TOTAL (--> MAXCTYPE) on which NBCTYPES is defined
	bool dothiscontenttype[NBCTYPES]; //booleans[0 to MAXCTYPE]
	bool isporosityloaded = false;
	bool istotvolumloaded = false;
	double *intern_steadyPorosA; //array of steady (constant-in-time) elemental EFFECTIVE porosities (though also playing the role of a 'total porosity' in most contexts, due to the presence of n_eff in the transient storage term of the ADE!
	double *intern_steadyTotVolumA; //array of steady (constant-in-time) elemental 'total volumes' (solid + saturated voids)
	bool steadyActiveElems = false; //tells if ALL elements remain active throughout the entire simulation (true), or else that the active-state of the elements should be updated-stored-verified during the simulation due to changes in element states (activated/deactivated) by the WCM mode

	bool *extern_transientACTIVelemA = nullptr; //array of transient elemental active/inactive states (which may be changed at some time steps by the host plugin
	double *intern_steadyXelemA; //array of steady (constant-in-time) elemental centroids.X position (local coord. syst.)
	double *intern_steadyYelemA; //array of steady (constant-in-time) elemental centroids.Y position (local coord. syst.)
	double *intern_steadyZelemA; //array of steady (constant-in-time) elemental centroids.Z position (local coord. syst.)
	double *intern_steadyDelemA; //array of steady (constant-in-time) elemental centroids.Depth below the refered 'top': calculated depth below top as specified in netcdftools_init_monitzones(...)
	int *monitzoneidsA = nullptr; //optional; set from within the main IFM plugin code, if the related elemental User Data is found; otherwise, it remains null while nbmonitzones = 0 still.
	int *uniqsortedmzidsA = nullptr;

	int *nodeindexA; //FEFLOW-node-index array; (values can change but not the length of the array); array of export nodes only (this is NOT a full array!)
	int *nodeifixedA; //ordered for-export node-index ('nodei') array; (fixed integer values starting with 0)
	int *nodeidvaluesA = nullptr; //NEW optional array of IDs for the export nodes, to be written once, as time-indep. data; this is indeed time-indep. because ONLY internal special idvalues actually correspond to non-fixed nodes... (e.g. the adaptative glob.top.nodes with Id = expNautotopId = -999 from PaleoSea2D.cpp)
	size_t nodeindex_len = 0; //length of nodeindex and nodeifixed arrays (this length must be kept constant through time!)
	//REMOVED: size_t Len_DimNodes = 0; //Length in the 'node' dimension (with an initial value to know when it is initialized); essentially for type conversion of nodeindex_len from int to size_t


	/* The start and count arrays will tell the netCDF library where to write our data. */
	size_t start[NDIMS2], count[NDIMS2];
	size_t CZstart[CZNDIMS2], CZcount[CZNDIMS2];
	size_t CZQstart[CZQNDIMS3], CZQcount[CZQNDIMS3];
	size_t MZstart[MZNDIMS2], MZcount[MZNDIMS2];
	size_t DFLUXstart[VECTCNDIMS3], DFLUXcount[VECTCNDIMS3];
	size_t NBNEGCstart[NBNEGCNDIMS2], NBNEGCcount[NBNEGCNDIMS2];

	//time record index (starts at 0; unlimited)
	size_t rec_index;
};

//Note: If the model contains elements which are initially inactive, please use an external array as alreadyV, which should contain the extracted IfmTOTAL_VOLUME elem. content for ALL elements (by tricking the IFM function by temporarily activating the i^th element in the extraction loop)
bool compute_all_element_centroids(IfmDocument pDoc, double *outX, double *outY, double *outZ = nullptr, double *alreadyV = nullptr)
{
	bool is3D = IfmGetNumberOfDimensions(pDoc) == IfmNDM_3D;
	bool doX = outX != nullptr;
	bool doY = outY != nullptr;
	bool doZ = is3D && outZ != nullptr;
	bool extV = alreadyV != nullptr; //use the external alreadyV array of total elem. volume data, if specified
	//TODO: Add warning msgs

	int nnpe = IfmGetNumberOfNodesPerElement(pDoc); //number of nodes per element (constant; at least assumed & required)
	int ne = IfmGetNumberOfElements(pDoc); //number of elements

	double Evolum; //current element's volume
	double *EnodcoordA;
	EnodcoordA = new double[nnpe];

	int ei; //current element index
	int j, cni; //current loop iterator inside current element, and actual node index
	for (ei = 0; ei < ne; ei++) {
		Evolum = extV ? alreadyV[ei] : IfmGetElementalContent(pDoc, IfmTOTAL_VOLUME, ei);
		if (doX) {
			for (j = 0; j < nnpe;j++) {
				cni = IfmGetNode(pDoc, ei, j);
				EnodcoordA[j] = IfmGetX(pDoc, cni);
			}
			outX[ei] = IfmIntegrateNodalQuantitiesOfElement(pDoc, ei, EnodcoordA) / Evolum;
		}
		if (doY) {
			for (j = 0; j < nnpe;j++) {
				cni = IfmGetNode(pDoc, ei, j);
				EnodcoordA[j] = IfmGetY(pDoc, cni);
			}
			outY[ei] = IfmIntegrateNodalQuantitiesOfElement(pDoc, ei, EnodcoordA) / Evolum;
		}
		if (doZ) {
			for (j = 0; j < nnpe;j++) {
				cni = IfmGetNode(pDoc, ei, j);
				EnodcoordA[j] = IfmGetZ(pDoc, cni);
			}
			outZ[ei] = IfmIntegrateNodalQuantitiesOfElement(pDoc, ei, EnodcoordA) / Evolum;
		}
	}
	if(is3D) IfmInfo(pDoc, "BETA: compute_all_element_centroids: done... but not yet VERIFIED for 3D models!");
	return true;
}

bool compute_quantiles(std::vector<double> &values, std::vector<double> &probs, std::vector<double> &quants_out)
{
	std::vector<double> sv; //sorted values
	sv = values; //assigns the values (~copy); works: http://www.cplusplus.com/reference/vector/vector/operator=/

	size_t nv = sv.size(); //number of values
	size_t np = probs.size(); //number of quantiles requested (= nb of probs)

	if (nv > 0)  std::sort(sv.begin(), sv.end()); //do sort the values so that sv indeed becomes sorted

    //vector<size_t> qiv;
	//qiv.resize(np);
	quants_out.resize(np);

	size_t i;
	size_t qvi;
	for (i = 0; i < np; i++) {
		qvi = (size_t)floor((nv - 1.0)*probs[i]);
		//qiv[i] = qvi;
		quants_out[i] = nv > 0 ? sv[qvi] : nodata_double_ncdf; //sv[qiv[i]];
	}

	return true;
}

template<class T>
size_t count_duplicates_in_vector(std::vector<T> &values, std::vector<int> &counts_out, bool alreadysorted = false)
{
	std::vector<T> sv; //sorted values
	sv = values; //assigns the values (~copy); works: http://www.cplusplus.com/reference/vector/vector/operator=/
	if (!alreadysorted) std::sort(sv.begin(), sv.end()); //do sort the values so that sv indeed becomes sorted

	size_t nv = sv.size(); //number of values

	counts_out.clear();

	if (nv == 0) return 0; //and counts_out is empty...

	T currval, prevval;
	int currcnt = 0;
	bool newval;

	prevval = sv[0];

	size_t i;
	for (i = 0; i < nv; i++) {
		currval = sv[i];
		newval = currval != prevval;
		if (newval) {
			counts_out.push_back(currcnt); //stores the count for the previous distinct value & its duplicates
			currcnt = 1; //initializes the count for the current (~new) distinct value & its duplicates
			prevval = currval;
		}
		else {
			currcnt++;
		}
	}
	counts_out.push_back(currcnt);

	return counts_out.size();
}

////BETA!
//bool netcdftools_extract_nodal_param(IfmDocument pDoc, ncdfheader &nh, IfmParamID nId, double* destpValarray)
//{
//	bool good;
//	int parsiz = IfmGetParamSize(pDoc, nId);
//	if (parsiz != nh.nbnodes_fixed) IfmWarning(pDoc, "[netcdftools_extract_nodal_param] STRANGELY, param.size != nb.Nodes !!?");
//	int nbread = IfmGetParamValues(pDoc, nId, destpValarray, 0, nh.nbnodes_fixed);
//	good = nbread == nh.nbnodes_fixed;
//	if (!good) IfmWarning(pDoc, "[netcdftools_extract_nodal_param] STRANGELY, nb.par.values.read != nb.Nodes !!?");
//	//IF DONE with no warning, should be OK then, and ready for use by other 'compute' functions.
//	return good;
//}

//Note: netcdftools_init_contentzones should be called BERFORE calling this function (as nh.dothiscontenttype is needed)
void netcdftools_init_steadydata(IfmDocument pDoc, ncdfheader &nh, bool forceporo = false, bool forcetotvolum = false)
{
	nh.isporosityloaded = false; //initial default
	bool isporosityneeded = nh.dothiscontenttype[IfmVOID_VOLUME] || nh.dothiscontenttype[IfmFLUID_CONTENT] || nh.putmassconvectmonit || forceporo;

	nh.istotvolumloaded = false;
	bool istotvolumneeded = ((nh.dothiscontenttype[IfmTOTAL_VOLUME] || isporosityneeded) && true) || nh.putmassconvectmonit || forcetotvolum; //RETRYING SOMETHING (TODO: If that works nicely, remove old comments, and deBETAize it): nh.steadyActiveElems; //The last condition is important, since inactive elements will return a volume of 0.0. Therefore, as we don't want to interfere with the model setup by temporarily activating all elements to extract their volume, we'll rather get (up-to-date, current) elem. volumes each time this info is needed. (reviewed on Dec. 16th, 2016)
	const bool BETAtmpactivElem = true; //NEW BETA necessary operation: seems to work properly!
	int tmpEactiv; //temporary storage of the activ. state in the initial model

	if (isporosityneeded) nh.intern_steadyPorosA = new double[nh.nbelems_fixed];
	if (istotvolumneeded) nh.intern_steadyTotVolumA = new double[nh.nbelems_fixed];

	int i;
	for (i = 0; i < nh.nbelems_fixed; i++) {
		/* Porosity (if needed) */
		if (isporosityneeded) {
			nh.intern_steadyPorosA[i] = IfmGetMatMassPorosity(pDoc, i);
		}
		/* Elemental volume (if needed) */
		/* 1. First, temporarily activating the element so that V != 0.0 (which is the usual returned value for inactive elements, since we use IfmGetElementalContent(...) */
		if (istotvolumneeded) {
			if (BETAtmpactivElem) {
				tmpEactiv = IfmGetMatElementActive(pDoc, i);
				if (tmpEactiv == 0) IfmSetMatElementActive(pDoc, i, 1);
			}
			/* 2. Extracting V */
			nh.intern_steadyTotVolumA[i] = IfmGetElementalContent(pDoc, IfmTOTAL_VOLUME, i);
			/* 3. Finally, setting the element back to its previous active/inactive state */
			if (BETAtmpactivElem && tmpEactiv == 0) {
				IfmSetMatElementActive(pDoc, i, 0);
			}
			//TODO++: Verify if that perturbates the initial values in the model parts which started as inactive...
		}
	}
	nh.isporosityloaded = isporosityneeded;
	nh.istotvolumloaded = istotvolumneeded;

	if (nh.putvconcgrads) {
		//X coordinates (nodal, full domain, active & inactive):
		nh.copyoffixednodalXcoordsA = new double[nh.nbnodes_fixed];
		if (!genFFtools_extract_nodal_param(pDoc, IfmP_MSH_X, nh.copyoffixednodalXcoordsA)) {
			delete[] nh.copyoffixednodalXcoordsA;
			nh.putvconcgrads = false;
			IfmWarning(pDoc, "[netcdftools_init_steadydata.putvconcgrads] Unsuccessful read of IfmP_MSH_X FF param. array, hence: disabling .putvconcgrads (= false).");
			return;
		}

		//Y coordinates (nodal, full domain, active & inactive):
		nh.copyoffixednodalYcoordsA = new double[nh.nbnodes_fixed];
		if (!genFFtools_extract_nodal_param(pDoc, IfmP_MSH_Y, nh.copyoffixednodalYcoordsA)) {
			delete[] nh.copyoffixednodalYcoordsA;
			nh.putvconcgrads = false;
			IfmWarning(pDoc, "[netcdftools_init_steadydata.putvconcgrads] Unsuccessful read of IfmP_MSH_X FF param. array, hence: disabling .putvconcgrads (= false).");
			return;
		}

		//temp. buffer storing computed vertical concentration gradient values... (nodal, only at export nodes)
		nh.compvconcgradsA = new double[nh.nodeindex_len];

		//SEEMS TO WORK PROPERLY: IfmInfo(pDoc, "BETA: [netcdftools_init_steadydata].putvconcgrads done.");
	}
}

void netcdftools_init_dynamicdatabuffers(IfmDocument pDoc, ncdfheader &nh)
{
	if (nh.putnbnegconcs || nh.putglobextrconcs || nh.putvconcgrads) {
		nh.copyoftransientnodalconcsA = new double[nh.nbnodes_fixed];
	}
}

//TODO: Use the functional bool output to decide if plugin should continue or not (by completing both this code and the plugin's main code)
//Note: netcdftools_init_steadydata should be called AFTER having called this function
bool netcdftools_init_contentzones(IfmDocument pDoc, ncdfheader &nh, const IfmContentType* types2compute, const int ntypes, const char* elCZrdname = nullptr, long elCZrdid = -1)
{
	int i;
	bool validCZrd = (elCZrdname != nullptr) || (elCZrdid >= 0);
	if (validCZrd && elCZrdid < 0) {
		elCZrdid = IfmGetElementalRefDistrIdByName(pDoc, elCZrdname);
		validCZrd = elCZrdid >= 0;
	}
	nh.contentzoneidsA = nullptr; //pre-initialization null pointer, to know if it must be deallocated or not if !validCZrd

	if (validCZrd) {
		/* Determine which contents should be computed */
		/* 1. Initialization of the array = all false */
		for (i = 0; i < NBCTYPES; i++) {
			nh.dothiscontenttype[i] = false;
		}
		/* 2. Assigning true to desired c. types */
		for (i = 0; i < ntypes; i++) {
			nh.dothiscontenttype[types2compute[i]] = true;
		}
		/* 3. Counting desired c. types */
		int nbuniqctypes = 0;
		for (i = 0; i < NBCTYPES; i++) {
			if (nh.dothiscontenttype[i]) nbuniqctypes++;
		}
		/* 4. Finally, disabling CZ monitoring if no c. type is true */
		if (nbuniqctypes < 1) validCZrd = false;
	}

	if (validCZrd) {
		std::vector<int> czidV(nh.nbelems_fixed); //vector of CZ ids
		double readzv; //CZ Id value for current i^th element

		/* importing zone id data */
		for (i = 0; i < nh.nbelems_fixed; i++) {
			/* ContentZone Id */
			readzv = IfmGetElementalRefDistrValue(pDoc, elCZrdid, i);
			czidV[i] = isnan(readzv) ? -1 : (int)floor(readzv + 0.5);
			if (czidV[i] < 0) czidV[i] = -1;
			//Therefore, values will be either valid integers >= 0, else -1... but never lower.
		}
		nh.contentzoneidsA = new int[nh.nbelems_fixed];
		std::copy(czidV.begin(), czidV.end(), nh.contentzoneidsA); //based on http://stackoverflow.com/questions/2923272/how-to-convert-vector-to-array-c or http://stackoverflow.com/questions/12866912/fastest-way-to-copy-the-contents-of-a-vector-into-an-array
		//Compiler build warning msg explanations for std copy: http://stackoverflow.com/questions/903064/compiler-error-function-call-with-parameters-that-may-be-unsafe
		//Anyway, we can ignore this warning with no stress, as both source vector and dest. array are sized with the same length = nh.nbelems_fixed!

		/* extracting sorted unique zone id values */
		//Based on http://en.cppreference.com/w/cpp/algorithm/unique
		// and http://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector
		std::sort(czidV.begin(), czidV.end());

		std::vector<int> nbelperCZ; //element count for each (sorted) CZ Id (aligned with nh.uniqsortedczidsA, thus)
		size_t nbCZconfirm = count_duplicates_in_vector(czidV, nbelperCZ, true);

		auto tmpV = std::unique(czidV.begin(), czidV.end());
		czidV.erase(tmpV, czidV.end());

		/* ...then we must remove the first value if it is -1 */
		if (czidV.size() > 0 && czidV[0] == -1) {
			czidV.erase(czidV.begin() + 0);
			nbelperCZ.erase(nbelperCZ.begin() + 0);
			nbCZconfirm--;
		}

		validCZrd = czidV.size() > 0;

		if (validCZrd) {
			nh.nbcontentzones = (int)czidV.size();
			nh.uniqsortedczidsA = new int[nh.nbcontentzones];
			std::copy(czidV.begin(), czidV.end(), nh.uniqsortedczidsA);
			//Compiler build warning msg related to std copy can be ignored here again (see explanations above). 

			char txtbuff[180];
			sprintf_s(txtbuff, 180, "[NetCDF] netcdftools_init_contentzones: %d different CZ IDs were detected.", nh.nbcontentzones);
			IfmInfo(pDoc, txtbuff);

			/* Allocate the summation arrays */
			// based on http://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new and http://stackoverflow.com/questions/14829105/2d-dynamic-memory-allocation-array-in-c
			nh.CZcomputedsums = new double*[NBCTYPES];
			for (i = 0; i < NBCTYPES; i++) {
				nh.CZcomputedsums[i] = nh.dothiscontenttype[i] ? new double[nh.nbcontentzones] : nullptr;
			}

			/* memorize the indeces of the elements of each CZ */
			if (nh.putconcquantiles) {
				/* Mostly copied from the code of function netcdftools_compute_contentzones(...) */
				int *foundczptr; //pointer to the found value
				int fndczintpos; //position in the internal array
				int currCZ;
				int *findfirst = nh.uniqsortedczidsA + 0;
				int *findlastplus1 = nh.uniqsortedczidsA + nh.nbcontentzones + 1;
				int cnterrs = 0;

				nh.CZ_dynA_Eindeces = new std::vector<int>[nh.nbcontentzones]; //based on: http://stackoverflow.com/a/29358808/3433903
				for (i = 0; i < nh.nbcontentzones; i++) {
					nh.CZ_dynA_Eindeces[i].reserve(nbelperCZ[i]); //should match exactly with final 'actual' size of the vectors
				}

				for (i = 0; i < nh.nbelems_fixed; i++) {
					currCZ = nh.contentzoneidsA[i];
					if (currCZ < 0) continue;
					//find, based on http://www.cplusplus.com/reference/algorithm/find/
					foundczptr = std::find(findfirst, findlastplus1, currCZ);
					if (foundczptr != findlastplus1) {
						fndczintpos = (int)(foundczptr - findfirst);
						nh.CZ_dynA_Eindeces[fndczintpos].push_back(i); //TODO: Consider some optimization here, since this repeated insertion at vector end is not cool at all!
					}
					else {
						cnterrs++;
					}
				}
				for (i = 0; i < nh.nbcontentzones; i++) {
					nh.CZ_dynA_Eindeces[i].shrink_to_fit();
				}

				nh.CZ_fixed_qprobsV.resize(HWnbProbDivs + 1);
				for (i = 0; i <= HWnbProbDivs; i++) { //Programmer's note: Note the '<=' here.
					nh.CZ_fixed_qprobsV[i] = (double)i / HWnbProbDivs;
				}

				nh.CZ_2DquantsA = new double[nh.nbcontentzones * (HWnbProbDivs + 1)];

				if (cnterrs > 0) IfmWarning(pDoc, "[NetCDF] There were errors due to invalid Content Zone Ids while preparing CZ elem. indeces storage...");
			}
		}
	}
	if (!validCZrd) {
		if (nh.contentzoneidsA != nullptr) {
			delete[] nh.contentzoneidsA;
			nh.contentzoneidsA = nullptr;
		}
		nh.contentzoneidsA = nullptr;
		nh.uniqsortedczidsA = nullptr;
		nh.nbcontentzones = 0;
		for (i = 0; i < NBCTYPES; i++) {
			nh.dothiscontenttype[i] = false;
		}
		nh.CZ_dynA_Eindeces = nullptr;
		IfmInfo(pDoc, "[NetCDF] netcdftools_init_contentzones: INACTIVE, since NO related user data is present, or it is all nodata.");
	}
	return true;
}

//TODO: Use the functional bool output to decide if plugin should continue or not (by completing both this code and the plugin's main code)
//(TODO-minor-longterm: make it compatible with 3D problems... at least the typical ones with structured mesh)
//Note: assign nullptr to unused char, and -1 to unused long arguments; netcdftools_init_steadydata must be called before
bool netcdftools_init_monitzones(IfmDocument pDoc, ncdfheader &nh, const char* elMZrdname, long elMZrdid, const char* eltoprdname, long eltoprdid)
{
	/* 2D problems only (X,Y) for now */
	bool valid;
	int i;
	char txtbuff[180];

	if (!nh.putmassconvectmonit) {
		nh.nbmonitzones = 0;
		IfmWarning(pDoc, "[NetCDF] netcdftools_init_monitzones was called although ncdf.putmassconvectmonit == false. Please improve the code to avoid this.");
		return false;
	}

	if (nh.ref_C0 < 0.0) {
		IfmWarning(pDoc, "[NetCDF] netcdftools_init_monitzones: ref_C0 was not set to >= 0.0 mg/L, in the plugin code.");
		IfmWarning(pDoc, "[NetCDF] netcdftools_init_monitzones: TODO program an Alert + Stop of the plugin right after this event.");
		return false;
	}
	if (nh.ref_CS <= nh.ref_C0) {
		IfmWarning(pDoc, "[NetCDF] netcdftools_init_monitzones: ref_CS was not set to a value > ref_C0, in the plugin code.");
		IfmWarning(pDoc, "[NetCDF] netcdftools_init_monitzones: TODO program an Alert + Stop of the plugin right after this event.");
		return false;
	}

	/* Prepare centroid (~barycenter) coordinates of ALL elements */
	if (!nh.istotvolumloaded) {
		IfmWarning(pDoc, "[NetCDF] netcdftools_init_monitzones: TotVolum was not loaded as expected. The plugin will likely bug very soon.");
		IfmWarning(pDoc, "[NetCDF] netcdftools_init_monitzones: TODO program an Alert + Stop of the plugin right after this event.");
		return false;
	}
	else {
		nh.intern_steadyXelemA = new double[nh.nbelems_fixed];
		nh.intern_steadyYelemA = new double[nh.nbelems_fixed];
		valid = compute_all_element_centroids(pDoc, nh.intern_steadyXelemA, nh.intern_steadyYelemA, nullptr, nh.intern_steadyTotVolumA);
	}

	/* Detecting the user Data ref. distribution for the local top of the aquifer of interest (~datum) (OPTIONAL) */
	bool validtoprd = (eltoprdname != nullptr) || (eltoprdid >= 0);
	if (validtoprd && eltoprdid < 0) {
		eltoprdid = IfmGetElementalRefDistrIdByName(pDoc, eltoprdname);
		validtoprd = eltoprdid >= 0;
	}
	/* importing top data, and calculating depths */
	double readtv; //top value for current i^th element
	nh.intern_steadyDelemA = new double[nh.nbelems_fixed]; //allocating the complementary array with calculated depth in the rock aquifer
	for (i = 0; i < nh.nbelems_fixed; i++) {
		/* top elevation data */
		readtv = validtoprd ? IfmGetElementalRefDistrValue(pDoc, eltoprdid, i) : 0.0; //Note: This way, it works even if not valid top ref.distrib. was specified... but it won't give usable | inferable depth values, except in very simple models with a flat constant top elevation
		/* calculated depth in aquifer of interest (below specified top) */
		nh.intern_steadyDelemA[i] = readtv - nh.intern_steadyYelemA[i];
	}
	if (validtoprd) {
		char tmptoprdname[60];
		IfmGetElementalRefDistrName(pDoc, eltoprdid, tmptoprdname);
		sprintf_s(txtbuff, 180, "[NetCDF] netcdftools_init_monitzones: elem. depths (at centroid) calculated using '%s' as top elev. ref. distrib.", tmptoprdname);
		IfmInfo(pDoc, txtbuff);
	}
	else {
		IfmWarning(pDoc, "[NetCDF] netcdftools_init_monitzones: elem. depths (at centroid) calculated using top=0.0m since no valid ref.distrib. was specified.");
	}

//COPIED###########################

	/* Detecting the user Data ref. distribution for the monit. zone Ids (REQUIRED!) */
	bool validMZrd = (elMZrdname != nullptr) || (elMZrdid >= 0);
	if (validMZrd && elMZrdid < 0) {
		elMZrdid = IfmGetElementalRefDistrIdByName(pDoc, elMZrdname);
		validMZrd = elMZrdid >= 0;
	}
	nh.monitzoneidsA = nullptr; //pre-initialization null pointer, to know if it must be deallocated or not if !validMZrd

	if (validMZrd) {
		std::vector<int> mzidV(nh.nbelems_fixed); //vector of MZ ids
		double readzv; //MZ Id value for current i^th element

		/* importing monit zone id data */
		for (i = 0; i < nh.nbelems_fixed; i++) {
			/* Monit Zone Id */
			readzv = IfmGetElementalRefDistrValue(pDoc, elMZrdid, i);
			mzidV[i] = isnan(readzv) ? -1 : (int)floor(readzv + 0.5);
			if (mzidV[i] < 0) mzidV[i] = -1;
			//Therefore, values will be either valid integers >= 0, else -1... but never lower.
		}
		nh.monitzoneidsA = new int[nh.nbelems_fixed];
		std::copy(mzidV.begin(), mzidV.end(), nh.monitzoneidsA); //see comments for similar call in netcdftools_init_contentzones procedure

	    /* extracting sorted unique zone id values (first, sorting the data) */
		std::sort(mzidV.begin(), mzidV.end());

		/* (then, as a bonus: counting nb of elements per MZ zone */
		std::vector<int> nbelperMZ; //element count for each (sorted) MZ Id (aligned with nh.uniqsortedmzidsA, thus)
		size_t nbMZconfirm = count_duplicates_in_vector(mzidV, nbelperMZ, true);

		/* (then, removing duplicates so that only the distinct values remain in the (sorted) vector */
		auto tmpV = std::unique(mzidV.begin(), mzidV.end());
		mzidV.erase(tmpV, mzidV.end());

		/* ...then we must remove the first value if it is -1 */
		if (mzidV.size() > 0 && mzidV[0] == -1) {
			mzidV.erase(mzidV.begin() + 0);
			nbelperMZ.erase(nbelperMZ.begin() + 0);
			nbMZconfirm--;
		}

		validMZrd = mzidV.size() > 0;

		/* Next: continue if everything looks valid and there is >0 monit zone */
		if (validMZrd) {
			nh.nbmonitzones = (int)mzidV.size();
			nh.uniqsortedmzidsA = new int[nh.nbmonitzones];
			std::copy(mzidV.begin(), mzidV.end(), nh.uniqsortedmzidsA);
			//Compiler build warning msg related to std copy can be ignored here again (see explanations above). 

			/* Lastly, create the small arrays which will store the computed values of interest */
			nh.MZ_vCOMdA = new double[nh.nbmonitzones];
			nh.MZ_DPFdA = new double[nh.nbmonitzones];

			sprintf_s(txtbuff, 180, "[NetCDF] netcdftools_init_monitzones: %d different MZ IDs were detected.", nh.nbmonitzones);
			IfmInfo(pDoc, txtbuff);

			/* Memorize the indeces of the elements of each CZ */
			/* (mostly copied from the code of function netcdftools_init_contentzones, late section for putquantiles) */
			int *foundMZptr; //pointer to the found value
			int fndMZintpos; //position in the internal array
			int currMZ;
			int *findfirst = nh.uniqsortedmzidsA + 0;
			int *findlastplus1 = nh.uniqsortedmzidsA + nh.nbmonitzones + 1;
			int cnterrs = 0;

			nh.MZ_dynA_Eindeces = new std::vector<int>[nh.nbmonitzones]; //based on: http://stackoverflow.com/a/29358808/3433903
			for (i = 0; i < nh.nbmonitzones; i++) {
				nh.MZ_dynA_Eindeces[i].reserve(nbelperMZ[i]); //should match exactly with final 'actual' size of the vectors
			}

			for (i = 0; i < nh.nbelems_fixed; i++) {
				currMZ = nh.monitzoneidsA[i];
				if (currMZ < 0) continue;
				//find, based on http://www.cplusplus.com/reference/algorithm/find/
				foundMZptr = std::find(findfirst, findlastplus1, currMZ);
				if (foundMZptr != findlastplus1) {
					fndMZintpos = (int)(foundMZptr - findfirst);
					nh.MZ_dynA_Eindeces[fndMZintpos].push_back(i); //TODO: Consider some optimization here, since this repeated insertion at vector end is not cool at all! >> NO MORE a problem, since I've reserved the appropriate size for the vector data.
				}
				else {
					cnterrs++;
				}
			}
			for (i = 0; i < nh.nbmonitzones; i++) {
				nh.MZ_dynA_Eindeces[i].shrink_to_fit();
			}

			if (cnterrs > 0) IfmWarning(pDoc, "[NetCDF] There were errors due to invalid Monit Zone Ids while preparing MZ elem. indeces storage...");

		}
	}
//END COPY
	if (!validMZrd) {
		//TODO: Write code very similar to that of initCZ... but adapted to the dyn. arrays created here.
		nh.putmassconvectmonit = false;
		nh.nbmonitzones = 0;
	}
	IfmWarning(pDoc, "TODO-BETA write code at end of initMZ to delete unused dyn. arrays if NOT validMZrd...");
	IfmWarning(pDoc, "TODO-BETA write code... also so that it informs (returns) false if failed to prepare MZ!");
	return true;
}

void netcdftools_compute_MZ_indicators(IfmDocument pDoc, ncdfheader &nh)
{
	if (nh.nbmonitzones == 0) return; //protection

	int mzi, ej, elemidx;
	double Mval; //elemental diluted-mass content value at the j^th element
	double PVval; //calculated porous volume within the j^th element (used to transform mass <--> concentration)
	double RMval; //relative diluted-mass above background conc. (i.e. M - C0*PV)
	double Dval; //depth of the centroid of the j^th element
	double Cval, NCval; //calculated conc. value; calc. normalized conc. (using C0 to CS bounds)
	std::vector<int> currEiV; //vector storing a temporary copy of the elem. indeces for the current MZ...

	/* accumulators (reset to zero at start of each MZ) */
	double cumRM; //cumulative summing of elem. relative Mass [g] (i.e. after subtracting C0*PV background conc.)
	double cumRMtD; //cum. summing of elem. rel.Mass*Depth [g*m]
	double updDmax; //updating through the inner loop, max. Depth [m] of element with calc. concentration > DPF_concmin.

	/* Shorter names to point to already allocated arrays; for computing Cval with cleaner code */
	double *tmpTotVolA = nh.intern_steadyTotVolumA;
	double *tmpPorosA = nh.intern_steadyPorosA;
	/* (other shorter names) */
	double Cmin = nh.ref_C0;
	double Cmax = nh.ref_CS;

	for (mzi = 0; mzi < nh.nbmonitzones; mzi++) {
		currEiV = nh.MZ_dynA_Eindeces[mzi];
		cumRM = 0.0;
		cumRMtD = 0.0;
		updDmax = -1.0e40; //very large negative number, so that it'll be replaced right away by an in-domain depth value
		for (ej = 0; ej < currEiV.size(); ej++) {
			elemidx = currEiV[ej];
			if (!nh.steadyActiveElems && !nh.extern_transientACTIVelemA[elemidx]) continue; //to skip inactive elements
			Mval = IfmGetElementalContent(pDoc, IfmDILUTED_MASS, elemidx);
			PVval = tmpTotVolA[elemidx] * tmpPorosA[elemidx];
			RMval = Mval - Cmin * PVval;
			Dval = nh.intern_steadyDelemA[elemidx];
			Cval = Mval / PVval;
			NCval = (Cval - Cmin) / (Cmax - Cmin);
			cumRM = cumRM + RMval;
			cumRMtD = cumRMtD + RMval*Dval;
			if (NCval > nh.DPF_normconcmin && Dval > updDmax) updDmax = Dval;
		}
		nh.MZ_vCOMdA[mzi] = cumRMtD / cumRM;
		nh.MZ_DPFdA[mzi] = updDmax;
	}
	/* END BETA DEMO CODE */
}

void netcdftools_compute_CZ_quantiles(IfmDocument pDoc, ncdfheader &nh)
{
	/* STRUCTURE VERY SIMILAR TO THAT OF netcdftools_compute_contentzones(...) */
	if (nh.nbcontentzones == 0) return; //protection

	/* 
	if (ctype != IfmDILUTED_MASS) { //other protection (may still bug before actually stopping the simulation)
		IfmAlert(pDoc, NULL, "  OK  ", "Unmanaged ctype != IfmDILUTED_MASS in the call netcdftools_compute_CZ_quantiles(...). Please rectify the plugin code!");
		IfmSetSimulationControlFlag(pDoc, IfmCTL_ABORT);
		return;
	} */

	//Shorter names to point to already allocated arrays
	//	double *destsumA = nh.CZcomputedsums[ctype];
	double *tmpTotVolA = nh.intern_steadyTotVolumA;
	double *tmpPorosA = nh.intern_steadyPorosA;

	int czi, ej, qj, elemidx;
	double dmassval, concval; //elemental diluted-mass content value at the j^th element; calculated conc. value
	std::vector<int> currEiV;
	std::vector<double> concsV;
	std::vector<double> calcquantsV;
	bool cquantok;

	for (czi = 0; czi < nh.nbcontentzones; czi++) {
		currEiV = nh.CZ_dynA_Eindeces[czi];
		/* First, compute & extract the current elemental concentrations */
		concsV.clear(); //Important so that the vector is cleared and again of zero length
		concsV.reserve(currEiV.size()); //Allocates the vector capacity to the max number of values it may have to stack
		for (ej = 0; ej < currEiV.size(); ej++) {
			elemidx = currEiV[ej];
			if (!nh.steadyActiveElems && !nh.extern_transientACTIVelemA[elemidx]) continue; //to skip inactive elements
			dmassval = IfmGetElementalContent(pDoc, IfmDILUTED_MASS, elemidx);
			concval = dmassval / (tmpTotVolA[elemidx] * tmpPorosA[elemidx]);
			// TODO: Verify if a protection should be added in case of inactive elements?...
			//REPLACED: concsV[ej] = concval;
			concsV.push_back(concval); //adding/appending the concval of the current ACTIVE element to the stack-vector
		}
		/* Then, compute the quantiles */
		cquantok = compute_quantiles(concsV, nh.CZ_fixed_qprobsV, calcquantsV);
		if (cquantok) {
			/* ...and copy the results into the 2D array */
			for (qj = 0; qj <= HWnbProbDivs; qj++) { //Programmer's note: Note the '<=' here.
				// array[i][j] is then rewritten as array[i*sizeY + j] ; source: http://stackoverflow.com/a/936709/3433903
				nh.CZ_2DquantsA[czi*(HWnbProbDivs + 1) + qj] = calcquantsV[qj];
			}
		}
		else {
			IfmWarning(pDoc, "netcdftools_compute_CZ_quantiles: UNMANAGED compute_quantiles == false!... Might bug soon.");
		}
	}
	//[LOOKS GOOD (so: no more Beta?)] IfmInfo(pDoc, "TMP BETA: netcdftools_compute_CZ_quantiles finished.");
}

/* The following function was adapted & verified to give adequate results for the following content types:
IfmTOTAL_VOLUME		Total volume
*IfmVOID_VOLUME		Void volume (but since it's a fully saturated problem, this is the volume of 'effective porosity' rather than 'total porosity'... at least until the workflow evolves towards Double Porosity modeling (with 'solid' phase conc. in the remaining major part of the total porosity)
*IfmFLUID_CONTENT	Fluid content (same as the Void volume, since it's a fully saturated problem with single liquid fluid
IfmDILUTED_MASS		Diluted mass (fluid phase)
 note: Types with an asterisk needed adaptation, since we work with fully saturated models with only an effective porosity. */
void netcdftools_compute_contentzones(IfmDocument pDoc, ncdfheader &nh, IfmContentType ctype)
{
	if (nh.nbcontentzones == 0) return; //protection
	int i;
	int *foundczptr; //pointer to the found value
	int fndczintpos; //position in the internal array
	int currCZ;
	double readcval;
	int *findfirst = nh.uniqsortedczidsA + 0;
	int *findlastplus1 = nh.uniqsortedczidsA + nh.nbcontentzones + 1;
	int cnterrs = 0;
	
	bool poroVmode = (ctype == IfmVOID_VOLUME) || (ctype == IfmFLUID_CONTENT);
	if (poroVmode && !nh.isporosityloaded) {
		poroVmode = false;
		IfmWarning(pDoc, "[NetCDF] netcdftools_compute_contentzones: Porosity was not loaded as expected. Computed VoidV and FluidC contents will not be adequate.");
	}

	bool totVmode = (ctype == IfmTOTAL_VOLUME);

	//Shorter names to point to already allocated arrays
	double *destsumA = nh.CZcomputedsums[ctype];
	double *tmpTotVolA = nh.intern_steadyTotVolumA;
	double *tmpPorosA = nh.intern_steadyPorosA;

	/* if (poroVmode && tmpTotVolA == nullptr) {
		poroVmode = false;
		IfmWarning(pDoc, "netcdftools_compute_contentzones: Strangely, tmpTotVolA is null. Computed VoidV and FluidC contents will not be adequate.");
	} */

	for (i = 0; i < nh.nbcontentzones; i++) {
		destsumA[i] = 0.0;
	}
	
	double poroVval;

	for (i = 0; i < nh.nbelems_fixed; i++) {
		currCZ = nh.contentzoneidsA[i];
		if (currCZ < 0) continue;
		if (!nh.steadyActiveElems && !nh.extern_transientACTIVelemA[i]) continue; //to skip currently inactive elements

		//find, based on http://www.cplusplus.com/reference/algorithm/find/
		foundczptr = std::find(findfirst, findlastplus1, currCZ);
		if (foundczptr != findlastplus1) {
			poroVval = poroVmode ? (nh.istotvolumloaded ? tmpTotVolA[i] : IfmGetElementalContent(pDoc, IfmTOTAL_VOLUME, i)) * tmpPorosA[i] : 0.0;
			readcval = poroVmode ? poroVval : (totVmode && nh.istotvolumloaded ? tmpTotVolA[i] : IfmGetElementalContent(pDoc, ctype, i));
			fndczintpos = (int)(foundczptr - findfirst);
			destsumA[fndczintpos] = destsumA[fndczintpos] + readcval;
		}
		else {
			cnterrs++;
		}
	}
	if (cnterrs > 0) IfmWarning(pDoc, "[NetCDF] There were errors due to invalid Content Zone Ids while computing zone sums...");
}


//This function must NOT be called if nbcontentzones == 0
int netcdftools_write_CZ_headings(IfmDocument pDoc, ncdfheader &nh)
{
	/* Error handling. */
	int retval;

	/* Define the dimensions for CZ. */
	if ((retval = nc_def_dim(nh.ncid, CZ_DIM_NAME, nh.nbcontentzones, &nh.CZ_dimid)))
		ERR(retval);

	/* Define the coordinate variables for CZ. */
	if ((retval = nc_def_var(nh.ncid, CZ_DIM_NAME, NC_INT, 1, &nh.CZ_dimid,
		&nh.CZindex_varid)))
		ERR(retval);

	/* Assign units attributes to coordinate variables. */
	if ((retval = nc_put_att_text(nh.ncid, nh.CZindex_varid, UNITS,
		strlen(CZ_UNITS), CZ_UNITS)))
		ERR(retval);

	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. Both of the netCDF variables we are
	creating share the same four dimensions. In C, the
	unlimited dimension must come first on the list of dimids. */
	nh.CZdimids[0] = nh.rectime_dimid;
	nh.CZdimids[1] = nh.CZ_dimid;

	nh.CZ_varidA = new int[NBCTYPES];

	int cti;
	int nbvarsdefined = 0;
	for (cti = 0; cti < NBCTYPES; cti++) {
		if (!nh.dothiscontenttype[cti]) continue;
		/* Define the netCDF variable for the "CTYPE[cti]" data. */
		if ((retval = nc_def_var(nh.ncid, CTYPE_namesA[cti], NC_DOUBLE, CZNDIMS2,
			nh.CZdimids, &nh.CZ_varidA[cti])))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.CZ_varidA[cti], UNITS,
			strlen(CTYPE_unitsA[cti]), CTYPE_unitsA[cti])))
			ERR(retval);
		nbvarsdefined++;
	}
	char txtbuffer[180];
	sprintf_s(txtbuffer, 180, "[NetCDF] netcdftools_write_CZ_headings: %d content-type output variables were created.", nbvarsdefined);
	IfmInfo(pDoc, txtbuffer);

	/* NOW, QUANTILES */
	if (nh.putconcquantiles) {
		/* Define the dimension qprob. */
		if ((retval = nc_def_dim(nh.ncid, QPROB_DIM_NAME, HWnbProbDivs + 1, &nh.qprob_dimid)))
			ERR(retval);

		/* Define the coordinate variable for qprob. */
		if ((retval = nc_def_var(nh.ncid, QPROB_DIM_NAME, NC_DOUBLE, 1, &nh.qprob_dimid,
			&nh.qprob_varid)))
			ERR(retval);

		/* Assign units attributes to the coordinate variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.qprob_varid, UNITS,
			strlen(QPROB_UNITS), QPROB_UNITS)))
			ERR(retval);

		nh.CZQdimids[0] = nh.rectime_dimid;
		nh.CZQdimids[1] = nh.CZ_dimid;
		nh.CZQdimids[2] = nh.qprob_dimid;

		/* Define the netCDF variable for the quantiles data. */
		if ((retval = nc_def_var(nh.ncid, CQUANT_NAME, NC_DOUBLE, CZQNDIMS3,
			nh.CZQdimids, &nh.cquant_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.cquant_varid, UNITS,
			strlen(CQUANT_UNITS), CQUANT_UNITS)))
			ERR(retval);

	}

	return NC_NOERR;
}

//This function should NOT be called if nbcontentzones == 0
int netcdftools_write_CZsums(IfmDocument pDoc, ncdfheader &nh)
{
	/* Error handling. */
	int retval;

	if (nh.nbcontentzones < 1) ERR(NC_NOERR);

	nh.CZstart[0] = nh.rec_index;

	int cti;
	for (cti = 0; cti < NBCTYPES; cti++) {
		if (!nh.dothiscontenttype[cti]) continue;
		if ((retval = nc_put_vara_double(nh.ncid, nh.CZ_varidA[cti], nh.CZstart, nh.CZcount,
			&nh.CZcomputedsums[cti][0])))
			ERR(retval);
	}
	//IfmInfo(pDoc, "BETA PROGRAMMING: netcdftools_write_CZsums done for this time step.");
	return NC_NOERR;
}

//This function must NOT be called if !nh.putconcquantiles
int netcdftools_write_CZquantiles(IfmDocument pDoc, ncdfheader &nh)
{
	/* Error handling. */
	int retval;

	if (nh.nbcontentzones < 1) ERR(NC_NOERR);

	nh.CZQstart[0] = nh.rec_index;

	if ((retval = nc_put_vara_double(nh.ncid, nh.cquant_varid, nh.CZQstart, nh.CZQcount,
		&nh.CZ_2DquantsA[0])))
		ERR(retval);
	//IfmInfo(pDoc, "BETA PROGRAMMING: netcdftools_write_CZsums done for this time step.");
	return NC_NOERR;
}

//This function must NOT be called if nbmonitzones == 0
int netcdftools_write_MZ_headings(IfmDocument pDoc, ncdfheader &nh)
{
	/* Error handling. */
	int retval;

	/* Define the dimensions for MZ. */
	if ((retval = nc_def_dim(nh.ncid, MZ_DIM_NAME, nh.nbmonitzones, &nh.MZ_dimid)))
		ERR(retval);

	/* Define the coordinate variable for MZ. */
	if ((retval = nc_def_var(nh.ncid, MZ_DIM_NAME, NC_INT, 1, &nh.MZ_dimid,
		&nh.MZindex_varid)))
		ERR(retval);

	/* Assign units attributes to the coordinate variable. */
	if ((retval = nc_put_att_text(nh.ncid, nh.MZindex_varid, UNITS,
		strlen(MZ_UNITS), MZ_UNITS)))
		ERR(retval);

	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. Both of the netCDF variables we are
	creating share the same four dimensions. In C, the
	unlimited dimension must come first on the list of dimids. */
	nh.MZdimids[0] = nh.rectime_dimid;
	nh.MZdimids[1] = nh.MZ_dimid;

	/* Define the netCDF variable for the vCOM_depth data. */
	if ((retval = nc_def_var(nh.ncid, VCOMD_NAME, NC_DOUBLE, MZNDIMS2,
		nh.MZdimids, &nh.vCOMd_varid)))
		ERR(retval);

	/* Assign units attributes to the netCDF variable. */
	if ((retval = nc_put_att_text(nh.ncid, nh.vCOMd_varid, UNITS,
		strlen(VCOMD_UNITS), VCOMD_UNITS)))
		ERR(retval);

	/* Define the netCDF variable for the DPF_depth data. */
	if ((retval = nc_def_var(nh.ncid, DPFD_NAME, NC_DOUBLE, MZNDIMS2,
		nh.MZdimids, &nh.DPFd_varid)))
		ERR(retval);

	/* Assign units attributes to the netCDF variable. */
	if ((retval = nc_put_att_text(nh.ncid, nh.DPFd_varid, UNITS,
		strlen(DPFD_UNITS), DPFD_UNITS)))
		ERR(retval);

	return NC_NOERR;
}

//This function must NOT be called if nbmonitzones == 0
int netcdftools_write_MZindicators(IfmDocument pDoc, ncdfheader &nh)
{
	/* Error handling. */
	int retval;

	if (nh.nbmonitzones == 0) ERR(NC_NOERR);

	nh.MZstart[0] = nh.rec_index;

	if ((retval = nc_put_vara_double(nh.ncid, nh.vCOMd_varid, nh.MZstart, nh.MZcount, &nh.MZ_vCOMdA[0])))
		ERR(retval);
	if ((retval = nc_put_vara_double(nh.ncid, nh.DPFd_varid, nh.MZstart, nh.MZcount, &nh.MZ_DPFdA[0])))
		ERR(retval);

	//[LOOKS GOOD (so: no more Beta?)] IfmInfo(pDoc, "BETA PROGRAMMING: netcdftools_write_MZindicators done for this time step.");
	return NC_NOERR;
}

int netcdftools_write_NBNEGC_headings(IfmDocument pDoc, ncdfheader &nh)
{
	/* Error handling. */
	int retval;

	/* Define the dimensions for NEGCUBND in preparation for NBNEGC. */
	if ((retval = nc_def_dim(nh.ncid, NEGCUBND_DIM_NAME, NbNegCuppBnds, &nh.negcubnd_dimid)))
		ERR(retval);

	/* Define the coordinate variable for NEGCUBND in preparation for NBNEGC. */
	if ((retval = nc_def_var(nh.ncid, NEGCUBND_DIM_NAME, NC_DOUBLE, 1, &nh.negcubnd_dimid,
		&nh.negcubnd_varid)))
		ERR(retval);

	/* Assign units attributes to the coordinate variable. */
	if ((retval = nc_put_att_text(nh.ncid, nh.negcubnd_varid, UNITS,
		strlen(NEGCUBND_UNITS), NEGCUBND_UNITS)))
		ERR(retval);

	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. Both of the netCDF variables we are
	creating share the same four dimensions. In C, the
	unlimited dimension must come first on the list of dimids. */
	nh.nbnegcdimids[0] = nh.rectime_dimid;
	nh.nbnegcdimids[1] = nh.negcubnd_dimid;

	/* Define the netCDF variable for the nb_negative_concs data. */
	if ((retval = nc_def_var(nh.ncid, NBNEGC_NAME, NC_INT, NBNEGCNDIMS2,
		nh.nbnegcdimids, &nh.nbnegc_varid)))
		ERR(retval);

	/* Assign units attributes to the netCDF variable. */
	if ((retval = nc_put_att_text(nh.ncid, nh.nbnegc_varid, UNITS,
		strlen(NBNEGC_UNITS), NBNEGC_UNITS)))
		ERR(retval);

	return NC_NOERR;
}


#define innanconcvalue isnan
/* REPLACED; TESTING
bool innanconcvalue(double cval)
{
	return isnan(cval) ? true : (cval > -1.0e-10 & cval < 1.0e-10);
}
*/ 
//IMPORTANT NOTE: netcdftools_extract_allnodalconcs must be run prior to this function!
void netcdftools_compute_nbnegconcs(IfmDocument pDoc, ncdfheader &nh)
{
	int j; //neg.conc.level index
	int i; //nodal index (nested inside j)
	double cval; //current[i] conc. value
	bool currcNaN; //is current[i] conc. value == NaN?

	double* nodalconcs = nh.copyoftransientnodalconcsA; //array of nodal concentrations for ALL nodes in their natural FEFLOW ordering; just a shorter name for the array pointer assigned...

	/* counting NaN's (only) */
	nh.compnbnegconcsA[0] = 0;
	for (i = 0; i < nh.nbnodes_fixed; i++) {
		currcNaN = innanconcvalue(nodalconcs[i]);
		if (currcNaN) nh.compnbnegconcsA[0]++;
	}

	/* counting negative values (excluding NaN's) */
	for (j = 1; j < NbNegCuppBnds; j++) {
		nh.compnbnegconcsA[j] = 0;
		for (i = 0; i < nh.nbnodes_fixed; i++) {
			currcNaN = innanconcvalue(nodalconcs[i]);
			cval = !currcNaN ? round(nodalconcs[i]) : nodata_double_ncdf;
			if (!currcNaN && cval < negCubounds[j]) nh.compnbnegconcsA[j]++;
		}
	}
	//SEEMS TO WORK PROPERLY: IfmInfo(pDoc, "BETA-temp: netcdftools_compute_nbnegconcs done.");
}

//IMPORTANT NOTE: netcdftools_extract_allnodalconcs must be run prior to this function!
void netcdftools_compute_globextrconcs(IfmDocument pDoc, ncdfheader &nh)
{
	int i; //nodal index (nested inside j)
	double cval; //current[i] conc. value
	bool currcNaN; //is current[i] conc. value == NaN?

	double* nodalconcs = nh.copyoftransientnodalconcsA; //array of nodal concentrations for ALL nodes in their natural FEFLOW ordering; just a shorter name for the array pointer assigned...

	double compgmin = 1e100;
	double compgmax = -1e100;
	for (i = 0; i < nh.nbnodes_fixed; i++) {
		cval = nodalconcs[i];
		currcNaN = innanconcvalue(cval);
		if (!currcNaN) {
			if (cval < compgmin)compgmin = cval;
			if (cval > compgmax)compgmax = cval;
		}
	}
	nh.comp_globminc = compgmin;
	nh.comp_globmaxc = compgmax;
	//SEEMS TO WORK PROPERLY: IfmInfo(pDoc, "BETA-temp: netcdftools_compute_globextrconcs done.");
}

void netcdftools_compute_vconcgrads(IfmDocument pDoc, ncdfheader &nh)
{
	double* nodalconcs = nh.copyoftransientnodalconcsA; //array of nodal concentrations for ALL nodes in their natural FEFLOW ordering; just a shorter name for the array pointer assigned...
	double* nodalXcoords = nh.copyoffixednodalXcoordsA; //array of nodal X coordinates for ALL nodes in their natural FEFLOW ordering; just a shorter name for the array pointer assigned...
	double* nodalYcoords = nh.copyoffixednodalYcoordsA; //array of nodal Y coordinates for ALL nodes in their natural FEFLOW ordering; just a shorter name for the array pointer assigned...
	
	double* compdeltac = nh.compvconcgradsA; //destination array for the computed values (pointer); shorter name for...
	const double deltaY = 1e-3; //delta vertical distance for calculating the v.gradient; in meters (currently equals 1 mm)

	double Cat, Cbelow, Xat, Yat, Ybelow;

	int i;
	int nodi;
	IfmBool ptindom; //is the point within the domain?
	for (i = 0; i < nh.nodeindex_len; i++) {
		nodi = nh.nodeindexA[i];
		Xat = nodalXcoords[nodi];
		Yat = nodalYcoords[nodi];
		Ybelow = Yat - deltaY;
		Cat = nodalconcs[nodi];
		Cbelow = IfmGetResultsTransportMassValueAtXYZ(pDoc, Xat, Ybelow, 0.0, &ptindom);
		compdeltac[i] = ptindom ? (Cat - Cbelow) / (Yat - Ybelow) : nodata_double_ncdf;
	}
	//SEEMS TO WORK PROPERLY: IfmInfo(pDoc, "BETA: netcdftools_compute_vconcgrads done.");
}

int netcdf_specify_export_nodes(IfmDocument pDoc, ncdfheader &nh, int *nodesA, int arraylen, bool fixednodes, int *nodeidvalues = nullptr, bool force_initial = false)
{
	bool initializing = (nh.nodeindex_len == 0) || force_initial;
	if (initializing || (arraylen == nh.nodeindex_len)) {
		nh.nodeindexA = nodesA; //TOCHG!
		nh.nodeidvaluesA = nodeidvalues;
		nh.nodeindex_len = arraylen;
		/* Prepare the nodei fixed array */
		if (initializing) {
			nh.arefixednodes = fixednodes;
			nh.nodeifixedA = new int[nh.nodeindex_len];
			int i;
			for (i = 0; i < nh.nodeindex_len; i++) {
				nh.nodeifixedA[i] = i;
			}
		}
		return NC_NOERR;
	}
	else {
		IfmWarning(pDoc, "[NetCDF] ERROR-to-be: Number of export nodes has just changed! That's not allowed! Stopping...");
		return NC_EINVALCOORDS;
	}
	//AND THEN, DON'T FORGET TO DEFINE ALSO the .puthead and others, right after calling this function...
}

//This function writes X,Y,nodeFF data only (for allowing writing this data as time-indep. if arefixednodes
int netcdf_writecoordsdata(IfmDocument pDoc, ncdfheader &nh)
{
	/* Error handling. */
	int retval;

	int i;
	int nodi;

	/* Extracting the new coordinate data at current export nodes (expN) */
	for (i = 0; i < nh.nodeindex_len; i++) {
		nodi = nh.nodeindexA[i];
		if (nh.putNodeFF) nh.Pnodeff_expN_outA[i] = nodi;
		if (nh.putX) nh.Pxcoord_expN_outA[i] = IfmGetX(pDoc, nodi);
		if (nh.putY) nh.Pycoord_expN_outA[i] = IfmGetY(pDoc, nodi);
		if (nh.putYoftoprock) nh.Pyoftoprock_expN_outA[i] = IfmGetNodalRefDistrValue(pDoc, nh.opt_Yoftoprock_rID, nodi);
	}

	size_t start1D = 0, count1D = nh.nodeindex_len;
	size_t *start_adapt = nh.arefixednodes ? &start1D : nh.start;
	size_t *count_adapt = nh.arefixednodes ? &count1D : nh.count;

	if (nh.putNodeFF) {
		//NodeFF:
		if ((retval = nc_put_vara_int(nh.ncid, nh.nodeff_varid, start_adapt, count_adapt,
			&nh.Pnodeff_expN_outA[0])))
			ERR(retval);
	}
	if (nh.putX) {
		//X:
		if ((retval = nc_put_vara_double(nh.ncid, nh.xcoord_varid, start_adapt, count_adapt,
			&nh.Pxcoord_expN_outA[0])))
			ERR(retval);
	}
	if (nh.putY) {
		//Y:
		if ((retval = nc_put_vara_double(nh.ncid, nh.ycoord_varid, start_adapt, count_adapt,
			&nh.Pycoord_expN_outA[0])))
			ERR(retval);
	}
	if (nh.putYoftoprock) {
		//Yoftoprock:
		if ((retval = nc_put_vara_double(nh.ncid, nh.yoftoprock_varid, start_adapt, count_adapt,
			&nh.Pyoftoprock_expN_outA[0])))
			ERR(retval);
	}

	return NC_NOERR;
}

int netcdf_create_file_headers(IfmDocument pDoc, ncdfheader &nh)
{
	//EXAMPLE code from: http://www.unidata.ucar.edu/software/netcdf/docs/examples1.html

	//HARDWIRE TO REMOVE SOON:
	//int Nnodes = IfmGetNumberOfNodes(pDoc);
	//nh.Len_DimNodes = Nnodes; //TOCHG!
	//nh.nodeindex_len = Nnodes; //TOCHG!

    /* Error handling. */
	int retval;

	/* Create the path (folder) if it does not exist */
	//TODO SOON+++ TO DEAL WITH THE QUITE FREQUENT CASE WHERE THE USER FORGETS TO CREATE A RESULTS FOLDER! +++++++++++
	//char tmpfullfilepath[_MAX_PATH];
	//char tmpfolderpath[_MAX_PATH];
	//PathCchRemoveFileSpec(tmpfolderpath, _MAX_PATH)
	//IfmMakeDirectory(pDoc, )

	/* Create the file. */                          //NEW TESTING of NC_SHARE! (20 nov.): seems to work Okay; brings the advantage of flushing the write buffer regularly, and of updating the (time)record count in the header, making it possible to read data for most if not all time steps even if the simulation is brutally terminated.
	if ((retval = nc_create(nh.filepath, NC_CLOBBER | NC_64BIT_OFFSET | NC_SHARE, &nh.ncid))) //BETA: NC_NETCDF4&NC_64BIT_OFFSET under test! Check compatibility with R!
		ERR(retval);

	/* Define the dimensions. The record dimension is defined to have
	* unlimited length - it can grow as needed. In this example it is
	* the time dimension.*/
	if ((retval = nc_def_dim(nh.ncid, NOD_DIM_NAME, nh.nodeindex_len, &nh.nod_dimid)))
		ERR(retval);
	if ((retval = nc_def_dim(nh.ncid, REC_DIM_NAME, NC_UNLIMITED, &nh.rectime_dimid)))
		ERR(retval);
	if ((retval = nc_def_dim(nh.ncid, VECTC_DIM_NAME, VECTC_FIXED_LEN, &nh.vectc_dimid)))
		ERR(retval);

	/* Define the coordinate variables. We will only define coordinate
	variables for [......].  Ordinarily we would need to provide
	an array of dimension IDs for each variable's dimensions, but
	since coordinate variables only have one dimension, we can
	simply provide the address of that dimension ID (&lat_dimid) and
	similarly for (&lon_dimid). */
	if ((retval = nc_def_var(nh.ncid, NOD_DIM_NAME, NC_INT, 1, &nh.nod_dimid,
		&nh.nodindex_varid))) //NOTE: The nodei (index) array is constant, since it only gives the index (1..N) within the array.
		ERR(retval);
	if ((retval = nc_def_var(nh.ncid, REC_DIM_NAME, NC_DOUBLE, 1, &nh.rectime_dimid,
		&nh.rectime_varid)))
		ERR(retval);
	if ((retval = nc_def_var(nh.ncid, VECTC_DIM_NAME, NC_INT, 1, &nh.vectc_dimid,
		&nh.vectc_varid)))
		ERR(retval);

	/* Assign units attributes to coordinate variables. */
	if ((retval = nc_put_att_text(nh.ncid, nh.nodindex_varid, UNITS,
		strlen(NODE_UNITS), NODE_UNITS)))
		ERR(retval);
	if ((retval = nc_put_att_text(nh.ncid, nh.rectime_varid, UNITS,
		strlen(TIME_UNITS), TIME_UNITS)))
		ERR(retval);
	if ((retval = nc_put_att_text(nh.ncid, nh.vectc_varid, UNITS,
		strlen(VECTC_UNITS), VECTC_UNITS)))
		ERR(retval);

	/* The dimids array is used to pass the dimids of the dimensions of
	the netCDF variables. Both of the netCDF variables we are
	creating share the same four dimensions. In C, the
	unlimited dimension must come first on the list of dimids. */
	nh.dimids[0] = nh.rectime_dimid;
	nh.dimids[1] = nh.nod_dimid;

	nh.Dflux_dimids[0] = nh.rectime_dimid;
	nh.Dflux_dimids[1] = nh.vectc_dimid;
	nh.Dflux_dimids[2] = nh.nod_dimid;

	if (nh.nodeidvaluesA != nullptr) {
		/* Define the netCDF variable for the nodeidvalue data. */
		if ((retval = nc_def_var(nh.ncid, NODIDV_NAME, NC_INT, 1,
			&nh.nod_dimid, &nh.nodidv_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.nodidv_varid, UNITS,
			strlen(NODIDV_UNITS), NODIDV_UNITS)))
			ERR(retval);
	}

	//Allocate memory for fixed-constant |OR| transient coordinate data for the export nodes (expN)
	if (nh.putNodeFF) nh.Pnodeff_expN_outA = new int[nh.nodeindex_len];
	if (nh.putX) nh.Pxcoord_expN_outA = new double[nh.nodeindex_len];
	if (nh.putY) nh.Pycoord_expN_outA = new double[nh.nodeindex_len];
	if (nh.putYoftoprock) {
		nh.opt_Yoftoprock_rID = IfmGetNodalRefDistrIdByName(pDoc, "Yoftoprock");
		if (nh.opt_Yoftoprock_rID < 0) {
			IfmWarning(pDoc, "[NetCDF] netcdf_create_file_headers: Sorry: forced to DISABLE .putYoftoprock, since Yoftoprock nodal user-data distribution is MISSING!");
			nh.putYoftoprock = false;
		}
		else {
			nh.Pyoftoprock_expN_outA = new double[nh.nodeindex_len];
		}
	}

	if (nh.putNodeFF) {
		/* Define the netCDF variable for the NodeFF data. */
		if ((retval = nc_def_var(nh.ncid, NODEFF_NAME, NC_INT, nh.arefixednodes ? 1 : NDIMS2,
			nh.arefixednodes ? &nh.nod_dimid : nh.dimids, &nh.nodeff_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.nodeff_varid, UNITS,
			strlen(NODEFF_UNITS), NODEFF_UNITS)))
			ERR(retval);
	}

	if (nh.putX) {
		/* Define the netCDF variable for the X data. */
		if ((retval = nc_def_var(nh.ncid, XCOORD_NAME, NC_DOUBLE, nh.arefixednodes ? 1 : NDIMS2,
			nh.arefixednodes ? &nh.nod_dimid : nh.dimids, &nh.xcoord_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.xcoord_varid, UNITS,
			strlen(XCOORD_UNITS), XCOORD_UNITS)))
			ERR(retval);
	}

	if (nh.putY) {
		/* Define the netCDF variable for the Yoftoprock data. */
		if ((retval = nc_def_var(nh.ncid, YCOORD_NAME, NC_DOUBLE, nh.arefixednodes ? 1 : NDIMS2,
			nh.arefixednodes ? &nh.nod_dimid : nh.dimids, &nh.ycoord_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.ycoord_varid, UNITS,
			strlen(YCOORD_UNITS), YCOORD_UNITS)))
			ERR(retval);
	}

	if (nh.putYoftoprock) {
		/* Define the netCDF variable for the Y data. */
		if ((retval = nc_def_var(nh.ncid, YOFTOPROCK_NAME, NC_DOUBLE, nh.arefixednodes ? 1 : NDIMS2,
			nh.arefixednodes ? &nh.nod_dimid : nh.dimids, &nh.yoftoprock_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.yoftoprock_varid, UNITS,
			strlen(YOFTOPROCK_UNITS), YOFTOPROCK_UNITS)))
			ERR(retval);
	}

	if (nh.puthead) {
		/* Define the netCDF variable for the Head data. */
		if ((retval = nc_def_var(nh.ncid, HEAD_NAME, NC_DOUBLE, NDIMS2,
			nh.dimids, &nh.head_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.head_varid, UNITS,
			strlen(HEAD_UNITS), HEAD_UNITS)))
			ERR(retval);
	}

	if (nh.putconc) {
		/* Define the netCDF variable for the Concentration data. */
		if ((retval = nc_def_var(nh.ncid, CONC_NAME, NC_DOUBLE, NDIMS2,
			nh.dimids, &nh.conc_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.conc_varid, UNITS,
			strlen(CONC_UNITS), CONC_UNITS)))
			ERR(retval);
	}

	if (nh.putvconcgrads) {
		/* Define the netCDF variable for the Vertical Concentration Gradient data. */
		if ((retval = nc_def_var(nh.ncid, VCONCGRAD_NAME, NC_DOUBLE, NDIMS2,
			nh.dimids, &nh.vconcgrad_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.vconcgrad_varid, UNITS,
			strlen(VCONCGRAD_UNITS), VCONCGRAD_UNITS)))
			ERR(retval);
	}

	if (nh.putglobextrconcs) {
		/* Define the netCDF variable for the Global Minimum & Maximum Concentration data. */
		if ((retval = nc_def_var(nh.ncid, GLOBMINCONC_NAME, NC_DOUBLE, 1, &nh.rectime_dimid,
			&nh.globminconc_varid)))
			ERR(retval);
		if ((retval = nc_def_var(nh.ncid, GLOBMAXCONC_NAME, NC_DOUBLE, 1, &nh.rectime_dimid,
			&nh.globmaxconc_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.globminconc_varid, UNITS,
			strlen(CONC_UNITS), CONC_UNITS)))
			ERR(retval);
		if ((retval = nc_put_att_text(nh.ncid, nh.globmaxconc_varid, UNITS,
			strlen(CONC_UNITS), CONC_UNITS)))
			ERR(retval);
	}

	if (nh.putnbnegconcs) {
		if ((retval = netcdftools_write_NBNEGC_headings(pDoc, nh))) ERR(retval);
	}

	if (nh.putDfluxcomponents) {
		/* Define the netCDF variable for the Darcy-flux-components data. */
		if ((retval = nc_def_var(nh.ncid, DFLUX_NAME, NC_DOUBLE, VECTCNDIMS3,
			nh.Dflux_dimids, &nh.Dflux_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.Dflux_varid, UNITS,
			strlen(DFLUX_UNITS), DFLUX_UNITS)))
			ERR(retval);

		nh.DFLUX_2DvaluesA = new double[VECTC_FIXED_LEN * nh.nodeindex_len];
	}

	if (nh.putmratebudg) {
		/* Define the netCDF variable for the Mass Rate Budget data. */
		if ((retval = nc_def_var(nh.ncid, MRBUDG_NAME, NC_DOUBLE, NDIMS2,
			nh.dimids, &nh.mrbudg_varid)))
			ERR(retval);

		/* Assign units attributes to the netCDF variable. */
		if ((retval = nc_put_att_text(nh.ncid, nh.mrbudg_varid, UNITS,
			strlen(MRBUDG_UNITS), MRBUDG_UNITS)))
			ERR(retval);
	}

	if (nh.nbcontentzones > 0) {
		if ((retval = netcdftools_write_CZ_headings(pDoc, nh))) ERR(retval);
	}

	if (nh.nbmonitzones > 0) {
		if ((retval = netcdftools_write_MZ_headings(pDoc, nh))) ERR(retval);
	}

	/* End define mode. */
	if ((retval = nc_enddef(nh.ncid)))
		ERR(retval);

	//======================================================================//
	/* Write the coordinate variable data (HERE ONLY IF FIXED export nodes) */
	

	if ((retval = nc_put_var_int(nh.ncid, nh.nodindex_varid, &nh.nodeifixedA[0])))
		ERR(retval);

	/* If nodes are fixed, coordinate variables can be written now, once for all. */
	if (nh.arefixednodes) {
		if ((retval = netcdf_writecoordsdata(pDoc, nh))) ERR(retval);
	}

	/* Write the "id value by user" variable data (only if available). */
	//NOTE: This array can indeed be written once for all, because ONLY special internal nodeidvalues correspond to non-fixed nodes... (see also the array's comment)
	if (nh.nodeidvaluesA != nullptr) {
		if ((retval = nc_put_var_int(nh.ncid, nh.nodidv_varid, &nh.nodeidvaluesA[0])))
			ERR(retval);
	}

	/* These settings tell netcdf to write one timestep of data. (The
	setting of start[0] inside the loop below tells netCDF which
	timestep to write.) */
	nh.count[0] = 1;
	nh.count[1] = nh.nodeindex_len;
	nh.start[1] = 0; //fixed, whereas [0] changes with each export time step, starting with 0 and increasing

	nh.rec_index = 0;

	if (nh.putDfluxcomponents) {
		/* Write the vector-component IDs time-INDEPendent variable data. */
		/* SKIPPED because it's not necessary to have these simplistic 1,2,3 values explicitly in the data file!
		if ((retval = nc_put_var_int(nh.ncid, nh.Dflux_varid, nh.CZ_fixed_qprobsV.data())))
			ERR(retval);
		*/
		
		// (0=time; 1=vector-component index; 2=node-as-ordered index)
		nh.DFLUXcount[0] = 1;
		nh.DFLUXcount[1] = VECTC_FIXED_LEN;
		nh.DFLUXcount[2] = nh.nodeindex_len;
		nh.DFLUXstart[1] = 0; //fixed, whereas [0] changes with each export time step, starting with 0 and increasing
		nh.DFLUXstart[2] = 0; //fixed
	}

	if (nh.putnbnegconcs) {
		/* Write the Upper Bounds for Negative Concentration counting (NBNEGC) time-INDEPendent variable data. */
		if ((retval = nc_put_var_double(nh.ncid, nh.negcubnd_varid, &negCubounds[0])))
			ERR(retval);

		nh.NBNEGCstart[1] = 0; //fixed, whereas [0] changes with each export time step, starting with 0 and increasing
		nh.NBNEGCcount[0] = 1;
		nh.NBNEGCcount[1] = NbNegCuppBnds;
	}

	if (nh.nbcontentzones > 0) {
		/* Write the content zone IDs time-INDEPendent variable data. */
		if ((retval = nc_put_var_int(nh.ncid, nh.CZindex_varid, &nh.uniqsortedczidsA[0])))
			ERR(retval);
		/* These settings tell netcdf to write one timestep of data. (The
		setting of CZstart[0] inside the loop below tells netCDF which
		timestep to write.) */
		nh.CZcount[0] = 1;
		nh.CZcount[1] = nh.nbcontentzones;
		nh.CZstart[1] = 0; //fixed, whereas [0] changes with each export time step, starting with 0 and increasing

		if (nh.putconcquantiles) {
			/* Write the content zone IDs time-INDEPendent variable data. */
			if ((retval = nc_put_var_double(nh.ncid, nh.qprob_varid, nh.CZ_fixed_qprobsV.data())))
				ERR(retval);

			// (0=time; 1=CZ; 2=quantile aligned with probability qprob)
			nh.CZQcount[0] = 1;
			nh.CZQcount[1] = nh.nbcontentzones;
			nh.CZQcount[2] = HWnbProbDivs + 1;
			nh.CZQstart[1] = 0; //fixed, whereas [0] changes with each export time step, starting with 0 and increasing
			nh.CZQstart[2] = 0; //fixed
		}
	}

	if (nh.nbmonitzones > 0) {
		/* Write the MonitZone IDs time-INDEPendent variable data. */
		if ((retval = nc_put_var_int(nh.ncid, nh.MZindex_varid, &nh.uniqsortedmzidsA[0])))
			ERR(retval);
		/* These settings tell netcdf to write one timestep of data. (The
		setting of MZstart[0] inside the loop below tells netCDF which
		timestep to write.) */
		nh.MZcount[0] = 1;
		nh.MZcount[1] = nh.nbmonitzones;
		nh.MZstart[1] = 0; //fixed, whereas [0] changes with each export time step, starting with 0 and increasing
	}

	return NC_NOERR;
}

int netcdf_writedata(IfmDocument pDoc, ncdfheader &nh)
{
	/* Error handling. */
	int retval;

	int i;
	int nodi;

	/* Program variables to hold the data we will write out. We will only
	need enough space to hold one timestep of data; one record. */
	//int *nodeff_out;
	//double *xcoord_out;
	//double *ycoord_out;
	double *head_out;
	double *conc_out;
	double *mrbudg_out;

	//if (!nh.arefixednodes) {
	//	if (nh.putNodeFF) nodeff_out = new int[nh.nodeindex_len];
	//	if (nh.putX) xcoord_out = new double[nh.nodeindex_len];
	//	if (nh.putY) ycoord_out = new double[nh.nodeindex_len];
	//}
	if (nh.puthead) head_out = new double[nh.nodeindex_len];
	if (nh.putconc) conc_out = new double[nh.nodeindex_len];
	if (nh.putmratebudg) mrbudg_out = new double[nh.nodeindex_len];

	bool extsrcmrb = nh.PtocompmrbudgetA != nullptr; //external source of mass rate budget data

	double tmpDfvcX, tmpDfvcY, tmpheadval; //tmpnohead is used for detecting if the node is currently deactivated (i.e. out of active model domain)
	bool isvelocpres = IfmIsVelocityFieldPresent(pDoc) == True;
	bool nodeisactive; //used along with tmpheadval...

	for (i = 0; i < nh.nodeindex_len; i++) {
		nodi = nh.nodeindexA[i];
		//if (!nh.arefixednodes) {
		//	if (nh.putNodeFF) nodeff_out[i] = nodi;
		//	if (nh.putX) xcoord_out[i] = IfmGetX(pDoc, nodi);
		//	if (nh.putY) ycoord_out[i] = IfmGetY(pDoc, nodi);
		//}
		if (nh.puthead) head_out[i] = IfmGetResultsFlowHeadValue(pDoc, nodi);
		if (nh.putconc) conc_out[i] = IfmGetResultsTransportMassValue(pDoc, nodi);
		if (nh.putDfluxcomponents) {
			tmpheadval = nh.puthead ? head_out[i] : IfmGetResultsFlowHeadValue(pDoc, nodi); //i.e. preferably the fast way, else the slow way
			nodeisactive = !isnan(tmpheadval); //TODO-later: Consider moving outside of this IF, for more general use of that info.
			tmpDfvcX = (isvelocpres & nodeisactive) ? IfmGetResultsXVelocityValue(pDoc, nodi) : nodata_double_ncdf;
			tmpDfvcY = (isvelocpres & nodeisactive) ? IfmGetResultsYVelocityValue(pDoc, nodi) : nodata_double_ncdf;
			//tmpDfvcZ = isvelocpres ? IfmGetResultsZVelocityValue(pDoc, nodi) : nodata_double_ncdf;
			nh.DFLUX_2DvaluesA[0 * nh.nodeindex_len + i] = tmpDfvcX;
			nh.DFLUX_2DvaluesA[1 * nh.nodeindex_len + i] = tmpDfvcY;
		}
		if (nh.putmratebudg) mrbudg_out[i] = extsrcmrb ? nh.PtocompmrbudgetA[nodi] : IfmGetResultsTransportMassValue(pDoc, nodi);
	}

	double rectime_out = IfmGetAbsoluteSimulationTime(pDoc);

	if (nh.putnbnegconcs || nh.putglobextrconcs || nh.putvconcgrads) {
		if (!genFFtools_extract_nodal_param(pDoc, IfmP_CONC, nh.copyoftransientnodalconcsA)) {
			IfmWarning(pDoc, "[netcdf_writedata] Unsuccessful read of IfmP_CONC FF param. array, hence exports nbnegconcs | globextrconcs | vconcgrads may be strange for this time step.");
		}
	}
	if (nh.putvconcgrads) {
		netcdftools_compute_vconcgrads(pDoc, nh);
	}
	if (nh.putglobextrconcs) {
		netcdftools_compute_globextrconcs(pDoc, nh);
	}
	if (nh.putnbnegconcs) {
		netcdftools_compute_nbnegconcs(pDoc, nh);
	}

	/* Write the pretend data. This will write our surface pressure and
	surface temperature data. The arrays only hold one timestep worth
	of data. We will just rewrite the same data for each timestep. In
	a real application, the data would change between timesteps. */
	nh.start[0] = nh.rec_index;
	//write current time
	if ((retval = nc_put_var1_double(nh.ncid, nh.rectime_varid, nh.start,
		&rectime_out)))
		ERR(retval);

	/* write desired 'put' variables */
	//if (!nh.arefixednodes) {
	//	if (nh.putNodeFF) {
	//		//NodeFF:
	//		if ((retval = nc_put_vara_int(nh.ncid, nh.nodeff_varid, nh.start, nh.count,
	//			&nodeff_out[0])))
	//			ERR(retval);
	//	}
	//	if (nh.putX) {
	//		//X:
	//		if ((retval = nc_put_vara_double(nh.ncid, nh.xcoord_varid, nh.start, nh.count,
	//			&xcoord_out[0])))
	//			ERR(retval);
	//	}
	//	if (nh.putY) {
	//		//Y:
	//		if ((retval = nc_put_vara_double(nh.ncid, nh.ycoord_varid, nh.start, nh.count,
	//			&ycoord_out[0])))
	//			ERR(retval);
	//	}
	//}
	/* If nodes are not fixed, coordinate variables are written here (at each export time step). */
	if (!nh.arefixednodes) {
		if ((retval = netcdf_writecoordsdata(pDoc, nh))) ERR(retval);
	}
	if (nh.puthead) {
		//head:
		if ((retval = nc_put_vara_double(nh.ncid, nh.head_varid, nh.start, nh.count,
			&head_out[0])))
			ERR(retval);
	}
	if (nh.putconc) {
		//conc:
		if ((retval = nc_put_vara_double(nh.ncid, nh.conc_varid, nh.start, nh.count,
			&conc_out[0])))
			ERR(retval);
	}
	if (nh.putvconcgrads) {
		//vertical.conc.gradient:
		if ((retval = nc_put_vara_double(nh.ncid, nh.vconcgrad_varid, nh.start, nh.count,
			&nh.compvconcgradsA[0])))
			ERR(retval);
	}
	if (nh.putglobextrconcs) {
		//global minimum & maximum concentration:
		if ((retval = nc_put_var1_double(nh.ncid, nh.globminconc_varid, nh.start,
			&nh.comp_globminc)))
			ERR(retval);
		if ((retval = nc_put_var1_double(nh.ncid, nh.globmaxconc_varid, nh.start,
			&nh.comp_globmaxc)))
			ERR(retval);
	}
	if (nh.putnbnegconcs) {
		//Number of nodal conc. values below the upper bounds
		nh.NBNEGCstart[0] = nh.rec_index;
		if ((retval = nc_put_vara_int(nh.ncid, nh.nbnegc_varid, nh.NBNEGCstart, nh.NBNEGCcount, &nh.compnbnegconcsA[0])))
			ERR(retval);
	}
	if (nh.putDfluxcomponents) {
		//Darcy fluxes (2D vector components):
		nh.DFLUXstart[0] = nh.rec_index;
		if ((retval = nc_put_vara_double(nh.ncid, nh.Dflux_varid, nh.DFLUXstart, nh.DFLUXcount,
			&nh.DFLUX_2DvaluesA[0])))
			ERR(retval);
	}
	if (nh.putmratebudg) {
		//conc:
		if ((retval = nc_put_vara_double(nh.ncid, nh.mrbudg_varid, nh.start, nh.count,
			&mrbudg_out[0])))
			ERR(retval);
	}
	if (nh.nbcontentzones > 0) {
		/* First, compute the content summations, for each active type */
		int cti;
		IfmContentType ctypi;
		for (cti = 0; cti < NBCTYPES; cti++) {
			ctypi = (IfmContentType)cti;
			if (!nh.dothiscontenttype[cti]) continue;
			netcdftools_compute_contentzones(pDoc, nh, ctypi);
			//Programmer's note: Depends on the order of constant enums of IfmContentType definition in IFM document.h
		}
		/* Then, write the results to the netCDF file */
		if ((retval = netcdftools_write_CZsums(pDoc, nh))) ERR(retval);

		/* NEXT (optional), compute quantiles of elemental concentration data within each content zone */
		if (nh.putconcquantiles) {
			netcdftools_compute_CZ_quantiles(pDoc, nh);
			netcdftools_write_CZquantiles(pDoc, nh);
		}
	}
	if (nh.nbmonitzones > 0) {
		/* First, call the computation of the diagnostic values (for each MZ zone) */
		netcdftools_compute_MZ_indicators(pDoc, nh);
		/* Then, write the results to the netCDF file */
		if ((retval = netcdftools_write_MZindicators(pDoc, nh))) ERR(retval);
	}

	nh.rec_index++;

	//if (!nh.arefixednodes) {
	//	if (nh.putNodeFF) delete[] nodeff_out;
	//	if (nh.putX) delete[] xcoord_out;
	//	if (nh.putY) delete[] ycoord_out;
	//}
	if (nh.puthead) delete[] head_out;
	if (nh.putconc) delete[] conc_out;
	if (nh.putmratebudg) delete[] mrbudg_out;

	return NC_NOERR;
}

int netcdf_close_file(IfmDocument pDoc, ncdfheader &nh)
{
	/* Error handling. */
	int retval;
	int i;

	/* Close the file. */
	if ((retval = nc_close(nh.ncid)))
		ERR(retval);

	char txtbuffer[_MAX_PATH + 40];
	sprintf_s(txtbuffer, _MAX_PATH + 40, "[NetCDF] Data were successfully exported to '%s'.", nh.filepath);
	IfmInfo(pDoc, txtbuffer);

	delete[] nh.nodeifixedA;

	/* Deallocate adaptative fixed|transient coordinate related arrays */
	if (nh.putNodeFF) delete[] nh.Pnodeff_expN_outA;
	if (nh.putX) delete[] nh.Pxcoord_expN_outA;
	if (nh.putY) delete[] nh.Pycoord_expN_outA;
	if (nh.putYoftoprock) delete[] nh.Pyoftoprock_expN_outA;

	/* Deallocate the Dflux related arrays */
	if (nh.putDfluxcomponents) {
		delete[] nh.DFLUX_2DvaluesA;
	}

	/* Deallocate the Content related arrays */
	if (nh.nbcontentzones > 0) {
		delete[] nh.contentzoneidsA;
		delete[] nh.uniqsortedczidsA;
		for (i = 0; i < NBCTYPES; i++) {
			if (nh.dothiscontenttype[i]) delete[] nh.CZcomputedsums[i];
		}
		delete[] nh.CZcomputedsums;
		if (nh.isporosityloaded) delete[] nh.intern_steadyPorosA;
		if (nh.istotvolumloaded) delete[] nh.intern_steadyTotVolumA;
		
		if (nh.putconcquantiles) {
			for (i = 0; i < nh.nbcontentzones; i++) {
				nh.CZ_dynA_Eindeces[i].clear();
			}
			delete[] nh.CZ_dynA_Eindeces;
			delete[] nh.CZ_2DquantsA;
		}
	}
	/* Deallocate the MonitZone related arrays */
	if (nh.nbmonitzones > 0) {
		delete[] nh.monitzoneidsA;
		delete[] nh.uniqsortedmzidsA;
		delete[] nh.MZ_vCOMdA;
		delete[] nh.MZ_DPFdA;
		if (true) {
			for (i = 0; i < nh.nbmonitzones; i++) {
				nh.MZ_dynA_Eindeces[i].clear();
			}
			delete[] nh.MZ_dynA_Eindeces;
		}
		//TODO: COMPLETE the deallocs !!!
	}
	/* Deallocate the MonitZone related arrays */
	if (nh.putnbnegconcs || nh.putglobextrconcs || nh.putvconcgrads) {
		delete[] nh.copyoftransientnodalconcsA;
	}
	if (nh.putvconcgrads) {
		delete[] nh.copyoffixednodalXcoordsA;
		delete[] nh.copyoffixednodalYcoordsA;
		delete[] nh.compvconcgradsA;
	}

	return NC_NOERR;
}
