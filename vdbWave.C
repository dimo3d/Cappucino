#include "vdbWave.h"
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
SOP_VdbWave::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "VDB";
	case 1: return "Points where active voxels should be";
	default: return "default";
	}
}

// define parameter for debug option
static PRM_Name debugPRM("debug", "Print debug information"); // internal name, UI name
static PRM_Name alphaPRM("alpha", "alpha value");

static PRM_Name dtPRM("dt", "delta time value");
static PRM_Name gravityPRM("gravity", "gravity value");
static PRM_Default alphaDefault(0.77);
static PRM_Default dtDefault(1.3);
static PRM_Default gravityDefault(9.83);

															  // assign parameter to the interface, which is array of PRM_Template objects
PRM_Template SOP_VdbWave::myTemplateList[] =
{
	PRM_Template(PRM_TOGGLE, 1, &debugPRM, PRMzeroDefaults), // type (checkbox), size (one in our case, but rgb/xyz values would need 3), pointer to a PRM_Name describing the parameter name, default value (0 - disabled)
	PRM_Template(PRM_FLT, 1, &dtPRM, &dtDefault),
	PRM_Template(PRM_FLT, 1, &alphaPRM, &alphaDefault),
	PRM_Template(PRM_FLT, 1, &gravityPRM, &gravityDefault),
	PRM_Template() // at the end there needs to be one empty PRM_Template object
};

// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbWave::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbWave(net, name, op);
}

SOP_VdbWave::SOP_VdbWave(OP_Network *net, const char *name, OP_Operator *op) : SOP_Node(net, name, op) {}

SOP_VdbWave::~SOP_VdbWave() {}

// function that does the actual job
OP_ERROR
SOP_VdbWave::cookMySop(OP_Context &context)
{
	// we must lock our inputs before we try to access their geometry, OP_AutoLockInputs will automatically unlock our inputs when we return
	OP_AutoLockInputs inputs(this);
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();

	// duplicate our incoming geometry
	duplicateSource(0, context);

	// check for interrupt - interrupt scope closes automatically when 'progress' is destructed.
	UT_AutoInterrupt progress("Activating voxels...");
	float _gravity = GRAVITY(context.getTime());
	float _dt = DT(context.getTime());
	float _alpha = ALPHA(context.getTime());
	// Get the first VDB primitive in the geometry
	_gravity *= _dt * _dt;
	int numberOfFoundVdbs = 0;
	GEO_PrimVDB* vdbPrim = NULL;
	const GEO_PrimVDB* verticalDerivativePrim = NULL;
	openvdb::GridBase::Ptr color_baseGrid;
	openvdb::Vec3SGrid::Ptr grid;
	openvdb::Vec3SGrid::Ptr grid_old;
	openvdb::GridBase::ConstPtr verticalDerivative_baseGrid;
	openvdb::Vec3SGrid::ConstPtr verticalDerivative_grid;


	for (GA_Iterator it(gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		GEO_Primitive* prim = gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<GEO_PrimVDB *>(prim))
		{
			vdbPrim = dynamic_cast<GEO_PrimVDB *>(prim);
			if (vdbPrim->hasGrid()) {

				color_baseGrid = vdbPrim->getGridPtr();
				grid = openvdb::gridPtrCast<openvdb::Vec3SGrid>(color_baseGrid);
				if (grid && (grid->getName() == "Cd")) {
					vdbPrim->makeGridUnique();
					grid = openvdb::gridPtrCast<openvdb::Vec3SGrid>(color_baseGrid);
					numberOfFoundVdbs += 1;
					break;
				}
			}
		}
	}
	for (GA_Iterator it(gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		GEO_Primitive* prim = gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<GEO_PrimVDB *>(prim))
		{
			vdbPrim = dynamic_cast<GEO_PrimVDB *>(prim);
			if (vdbPrim->hasGrid()) {

				color_baseGrid = vdbPrim->getGridPtr();
				grid_old = openvdb::gridPtrCast<openvdb::Vec3SGrid>(color_baseGrid);
				if (grid_old && (grid_old->getName() == "Cd_old")) {
					numberOfFoundVdbs += 1;
					//vdbPrim->makeGridUnique();
					grid_old = openvdb::gridPtrCast<openvdb::Vec3SGrid>(color_baseGrid);
					break;
				}
			}
		}
	}
	const GU_Detail *verticalDerivative_gdp = inputGeo(1);
	for (GA_Iterator it(verticalDerivative_gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		const GEO_Primitive* prim = verticalDerivative_gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<const GEO_PrimVDB *>(prim))
		{
			verticalDerivativePrim = dynamic_cast<const GEO_PrimVDB *>(prim);
			if (verticalDerivativePrim->hasGrid()) {
				verticalDerivative_baseGrid = verticalDerivativePrim->getConstGridPtr();
				verticalDerivative_grid = openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(verticalDerivative_baseGrid);
				if (verticalDerivative_grid) {
					numberOfFoundVdbs += 1;
					break;
				}
			}
		}
	}
	// Make sure we got a valid prim
	if ((!vdbPrim) || (!verticalDerivativePrim))
	
	{
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error(); 
	}
	
	// Docs say to do this in case the grid is shared



	// Try to get the vdbs grid

	if (!grid) {
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error(); 
	}

	if (!verticalDerivative_grid) {
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error();
	}
	if (!grid_old) {
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error();
	}

	openvdb::Vec3SGrid::ConstAccessor height_old_accessor = grid_old->getConstAccessor();
	openvdb::Vec3SGrid::ConstAccessor verticalDerivative_accessor = verticalDerivative_grid->getConstAccessor();
	openvdb::Vec3SGrid::ConstPtr grid_buffer = grid->deepCopy();

	float adt = _alpha * _dt;


	struct iWaveOp {
		openvdb::Vec3SGrid::ConstAccessor verticalDerivative_accessor;
		openvdb::Vec3SGrid::ConstAccessor height_old_accessor;
		float _gravity;
		float adt;
		float adt2;
		float twoMinusAdt;

		iWaveOp(openvdb::Vec3SGrid::ConstAccessor v, openvdb::Vec3SGrid::ConstAccessor c, float grav, float adt) :verticalDerivative_accessor(v),
			height_old_accessor(c), _gravity(grav), adt(adt) {
			adt2 = 1.0 / (1.0 + adt);
			twoMinusAdt = 2.0 - adt;

		}
		void operator()(const openvdb::Vec3SGrid::ValueOnIter iter) const {
			//for (openvdb::Vec3SGrid::ValueOnIter iter = grid->beginValueOn(); iter.test(); ++iter) {
			openvdb::Vec3f height = iter.getValue();
			const openvdb::Vec3f heightOld = height_old_accessor.getValue(iter.getCoord());

			height *= twoMinusAdt;
			height -= heightOld;
			height -= _gravity * verticalDerivative_accessor.getValue(iter.getCoord());
			height *= adt2;

			iter.setValue(height);

		}
	};
	openvdb::tools::foreach(grid->beginValueOn(), iWaveOp(verticalDerivative_accessor, height_old_accessor, _gravity, adt), false, 1);
	grid_old->setTree(grid_buffer->deepCopy()->treePtr());
	return error();

}
