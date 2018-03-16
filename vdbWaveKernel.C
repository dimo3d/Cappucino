#include "vdbWaveKernel.h"
#include <limits.h>
#include <SYS/SYS_Math.h>


#include <UT/UT_Interrupt.h>

#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>

#include <GU/GU_Detail.h>
#include <GEO/GEO_PrimPoly.h>

#include <PRM/PRM_Include.h>
#include <CH/CH_LocalVariable.h>

#include <OP/OP_AutoLockInputs.h>

#include <GU/GU_PrimVDB.h>

#include <openvdb/openvdb.h>

using namespace VdbCappucino;

// label node inputs, 0 corresponds to first input, 1 to the second one
const char *
SOP_VdbWaveKernel::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "VDB";
	default: return "default";
	}
}

double mycubic(double interp)
{
	if (interp < 0) return 1.0;
	if (interp > 1) return 0;
	double squared = interp * interp;
	return 2 * squared * interp - 3 * squared + 1;
}
// define parameter for debug option
static PRM_Name debugPRM("debug", "Print debug information"); // internal name, UI name
static PRM_Name sigmaPRM("sigma", "sigma value");
static PRM_Name dimPRM("dim", "Dimension value");
static PRM_Name dkPRM("dk", "dk value");
static PRM_Name endkPRM("endk", "endk value");
static PRM_Default sigmaDefault(0.77);
static PRM_Default dkDefault(1.3);
static PRM_Default endkDefault(7.7);
static PRM_Default dimDefault(6);
															  // assign parameter to the interface, which is array of PRM_Template objects
PRM_Template SOP_VdbWaveKernel::myTemplateList[] =
{
	PRM_Template(PRM_TOGGLE, 1, &debugPRM, PRMzeroDefaults), // type (checkbox), size (one in our case, but rgb/xyz values would need 3), pointer to a PRM_Name describing the parameter name, default value (0 - disabled)
	PRM_Template(PRM_INT, 1, &dimPRM, &dimDefault),
	PRM_Template(PRM_FLT, 1, &dkPRM, &dkDefault),
	PRM_Template(PRM_FLT, 1, &endkPRM, &endkDefault),
	PRM_Template(PRM_FLT, 1, &sigmaPRM, &sigmaDefault),
	PRM_Template() // at the end there needs to be one empty PRM_Template object
};

// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbWaveKernel::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbWaveKernel(net, name, op);
}

SOP_VdbWaveKernel::SOP_VdbWaveKernel(OP_Network *net, const char *name, OP_Operator *op) : SOP_Node(net, name, op) {}

SOP_VdbWaveKernel::~SOP_VdbWaveKernel() {}

// function that does the actual job
OP_ERROR
SOP_VdbWaveKernel::cookMySop(OP_Context &context)
{
	// we must lock our inputs before we try to access their geometry, OP_AutoLockInputs will automatically unlock our inputs when we return
	OP_AutoLockInputs inputs(this);
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();
	float sigma = SIGMA(context.getTime());
	float dk = DK(context.getTime());
	float endk = ENDK(context.getTime());
	int dim = DIM(context.getTime());
	// duplicate our incoming geometry
	duplicateSource(0, context);

	// check for interrupt - interrupt scope closes automatically when 'progress' is destructed.
	UT_AutoInterrupt progress("Activating voxels...");
	
	
	GEO_PrimVDB* vdbPrim = NULL; // empty pointer to vdb primitive

								 // iterate over all incoming primitives and find the first one which is VDB
	for (GA_Iterator it(gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		GEO_Primitive* prim = gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<const GEO_PrimVDB *>(prim))
		{
			vdbPrim = dynamic_cast<GEO_PrimVDB *>(prim);
			break;
		}
	}

	// terminate if volume is not VDB
	if (!vdbPrim)
	{
		addError(SOP_MESSAGE, "First input must contain a VDB!");
		return error();
	}

	// volume primitives in different nodes in Houdini by default share the same volume tree (for memory optimization) this will make sure that we will have our own deep copy of volume tree which we can write to 
	vdbPrim->makeGridUnique();

	// get grid base pointer and cast it to float grid pointer
	openvdb::GridBase::Ptr vdbPtrBase = vdbPrim->getGridPtr();
	openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(vdbPtrBase);

	// get accessor to the float grid
	openvdb::FloatGrid::Accessor vdb_access = grid->getAccessor();

	// get a reference to transformation of the grid
	const openvdb::math::Transform &vdbGridXform = grid->transform();

	// loop over all the points by handle
	//double dk = 0.1;
	//double sigma = 1.0;
	double norm = 0;
	double startK = dk;
	//double endK = 15;

	for (double freq = startK; freq < endk; freq += dk)
		// the original iWave kernel
		norm += freq * freq * exp(-sigma * freq * freq);
	// Try to get the vdbs grid

	if (!grid) {
		addError(SOP_MESSAGE, "First input must contain a float grid!");
		return error();
	}
	openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
	// Compute the signed distance from the surface of the sphere of each
	// voxel within the bounding box and insert the value into the grid
	// if it is smaller in magnitude than the background value.
	openvdb::Coord ijk;
	float weight = 0;
	int &i = ijk[0], &j = ijk[1], &k = ijk[2];
	for (i = -dim; i < +dim; ++i) {
		const double x2 = openvdb::math::Pow2(i);
		for (j = -dim; j < +dim; ++j) {
			const double x2y2 = openvdb::math::Pow2(j) + x2;
			for (k = -dim; k < +dim; ++k) {
				if (progress.wasInterrupted())
					return error();
				const double r = openvdb::math::Sqrt(x2y2 + openvdb::math::Pow2(k));
				double kern = 0;
				if (r == 0) {
					accessor.setValue(ijk, 1.0f);
					continue;
				}
				for (double freq = startK; freq < endk; freq += dk)
				{
					double currentSinc = sin(r * freq) / (r * freq);
					kern += freq * freq * freq * exp(-sigma * freq * freq) * currentSinc;
				}

				double interp = mycubic(((r / dim) - 0.9) / 0.1);
				kern *= (interp / (M_PI*norm));
				weight += kern;
				accessor.setValue(ijk, kern);
			}
		}
	}
	i = 0;
	j = 0;
	k = 0;
	float sumOfY = 0.0f;
	for (j = -dim; j < +dim; ++j) { sumOfY += accessor.getValue(ijk); };
	j = 0;
	if (DEBUG()) {
		printf("weight %f, norm %f\n", weight, norm);
	}
	accessor.setValue(ijk, 2.0 - sumOfY);
	


	return error();
}
