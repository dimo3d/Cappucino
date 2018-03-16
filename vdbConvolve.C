#include "vdbConvolve.h"
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
SOP_VdbConvolve::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "VDB";
	default: return "default";
	}
}

// define parameter for debug option
static PRM_Name debugPRM("debug", "Print debug information"); // internal name, UI name


															  // assign parameter to the interface, which is array of PRM_Template objects
PRM_Template SOP_VdbConvolve::myTemplateList[] =
{
	PRM_Template(PRM_TOGGLE, 1, &debugPRM, PRMzeroDefaults), // type (checkbox), size (one in our case, but rgb/xyz values would need 3), pointer to a PRM_Name describing the parameter name, default value (0 - disabled)
	
	PRM_Template() // at the end there needs to be one empty PRM_Template object
};

// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbConvolve::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbConvolve(net, name, op);
}

SOP_VdbConvolve::SOP_VdbConvolve(OP_Network *net, const char *name, OP_Operator *op) : SOP_Node(net, name, op) {}

SOP_VdbConvolve::~SOP_VdbConvolve() {}

// function that does the actual job
OP_ERROR
SOP_VdbConvolve::cookMySop(OP_Context &context)
{
	// we must lock our inputs before we try to access their geometry, OP_AutoLockInputs will automatically unlock our inputs when we return
	OP_AutoLockInputs inputs(this);
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();
	int numberOfFoundVdbs = 0;
	// duplicate our incoming geometry
	duplicateSource(0, context);

	// check for interrupt - interrupt scope closes automatically when 'progress' is destructed.
	UT_AutoInterrupt progress("Activating voxels...");



	GEO_PrimVDB* vdbPrim = NULL; // empty pointer to vdb primitive
	openvdb::GridBase::Ptr color_baseGrid;
	openvdb::Vec3SGrid::Ptr grid = NULL;
	// iterate over all incoming primitives and find the first one which is VDB
	for (GA_Iterator it(gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		GEO_Primitive* prim = gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<GEO_PrimVDB *>(prim))
		{
			vdbPrim = dynamic_cast<GEO_PrimVDB *>(prim);
			if (vdbPrim->hasGrid()) {
				vdbPrim->makeGridUnique();
				color_baseGrid = vdbPrim->getGridPtr();
				grid = openvdb::gridPtrCast<openvdb::Vec3SGrid>(color_baseGrid);
				if (grid && (grid->getName() == "Cd")) {
					numberOfFoundVdbs += 1;
					break;
				}
			}
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
		
	
	// get a reference to transformation of the grid
	const openvdb::math::Transform &vdbGridXform = grid->transform();

	// get pointer to geometry from second input
	const GU_Detail *kernel_gdp = inputGeo(1);
	const GEO_PrimVDB* kernelPrim = NULL;
	openvdb::GridBase::ConstPtr kernel_baseGrid;
	openvdb::FloatGrid::ConstPtr kernel_grid;
	for (GA_Iterator it(kernel_gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		const GEO_Primitive* prim = kernel_gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<const GEO_PrimVDB *>(prim))
		{
			kernelPrim = dynamic_cast<const GEO_PrimVDB *>(prim);
			if (kernelPrim->hasGrid()) {

				kernel_baseGrid = kernelPrim->getConstGridPtr();
				kernel_grid = openvdb::gridConstPtrCast<openvdb::FloatGrid>(kernel_baseGrid);
				if (kernel_grid) {
					
					break;
				}
			}
		}
	}
	// Make sure we got a valid prim
	if ((!vdbPrim) || (!kernelPrim))
	{
		if (DEBUG())
		printf("Convolve Number of Found Vdbs");
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error();
	}
	if (!grid) {
		addError(SOP_MESSAGE, "Input geometry must contain a Color VDB");
		return error();
	};
	if (!kernel_grid) {
		addError(SOP_MESSAGE, "second input geometry must contain a Float VDB");
		return error();
	};
	openvdb::Vec3SGrid::ConstPtr grid_buffer = grid->deepCopy();
	openvdb::Vec3SGrid::ConstAccessor color_accessor = grid_buffer->getConstAccessor();

	struct Convolve {
		openvdb::FloatGrid::ConstPtr kernel_grid;
		openvdb::Vec3SGrid::ConstAccessor color_accessor;
		Convolve(openvdb::FloatGrid::ConstPtr g, openvdb::Vec3SGrid::ConstAccessor c) :kernel_grid(g), color_accessor(c) {}
		inline void operator()(const openvdb::Vec3SGrid::ValueOnIter iter) const {
			openvdb::Vec3f temp = openvdb::Vec3f(0.0f, 0.0f, 0.0f);
			float weight = 0;
			for (openvdb::FloatGrid::ValueOnCIter kernel_iter = kernel_grid->cbeginValueOn(); kernel_iter.test(); ++kernel_iter) {

				temp += color_accessor.getValue(iter.getCoord() + kernel_iter.getCoord())*kernel_iter.getValue();
				weight += kernel_iter.getValue();
			}

			iter.setValue(temp);
		}
	};


	// Iterate over all active values.
	openvdb::tools::foreach(grid->beginValueOn(), Convolve(kernel_grid, color_accessor), false, 0);
	
	return error();
}
