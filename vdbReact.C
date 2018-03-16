#include "vdbReact.h"
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
SOP_VdbReact::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "VDB";
	case 1: return "Points where active voxels should be";
	case 2: return "VDB ";
	default: return "default";
	}
}

// define parameter for debug option
static PRM_Name debugPRM("debug", "Print debug information"); // internal name, UI name
static PRM_Name feedPRM("feed", "feed value");
static PRM_Name killPRM("kill", "kill value");
static PRM_Name deltaPRM("delta", "delta value");
static PRM_Name diffratePRM("diffrate", "diffrate value");
static PRM_Default feedDefault(0.07);
static PRM_Default killDefault(0.07);
static PRM_Default deltaDefault(1.00);
static PRM_Default diffrateDefault(0.25);
															  // assign parameter to the interface, which is array of PRM_Template objects
PRM_Template SOP_VdbReact::myTemplateList[] =
{
	PRM_Template(PRM_TOGGLE, 1, &debugPRM, PRMzeroDefaults), // type (checkbox), size (one in our case, but rgb/xyz values would need 3), pointer to a PRM_Name describing the parameter name, default value (0 - disabled)
	PRM_Template(PRM_FLT, 1, &feedPRM, &feedDefault),
	PRM_Template(PRM_FLT, 1, &killPRM, &killDefault),
	PRM_Template(PRM_FLT, 1, &deltaPRM, &deltaDefault),
	PRM_Template(PRM_FLT, 1, &diffratePRM, &diffrateDefault),
	PRM_Template() // at the end there needs to be one empty PRM_Template object
};

// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbReact::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbReact(net, name, op);
}

SOP_VdbReact::SOP_VdbReact(OP_Network *net, const char *name, OP_Operator *op) : SOP_Node(net, name, op) {}

SOP_VdbReact::~SOP_VdbReact() {}

// function that does the actual job
OP_ERROR
SOP_VdbReact::cookMySop(OP_Context &context)
{
	// we must lock our inputs before we try to access their geometry, OP_AutoLockInputs will automatically unlock our inputs when we return
	OP_AutoLockInputs inputs(this);
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();
	
	// duplicate our incoming geometry
	duplicateSource(0, context);
	
	// check for interrupt - interrupt scope closes automatically when 'progress' is destructed.
	UT_AutoInterrupt progress("Activating voxels...");

	// get pointer to geometry from second input
	GEO_PrimVDB* vdbPrim = NULL;
	const GEO_PrimVDB* distancePrim = NULL;
	openvdb::GridBase::Ptr color_baseGrid;
	openvdb::Vec3SGrid::Ptr grid;
	openvdb::GridBase::ConstPtr distance_baseGrid;
	openvdb::FloatGrid::ConstPtr distance_grid;
	int numberOfFoundVdbs = 0;
	const GU_Detail *distance_gdp = inputGeo(1);
	// Get the first VDB primitive in the geometry

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
				if ((grid) && (grid->getName() == "Cd")) {
					numberOfFoundVdbs += 1;
					break;
				}
			}
		}
	}

	for (GA_Iterator it(distance_gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		const GEO_Primitive* prim = distance_gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<const GEO_PrimVDB *>(prim))
		{
			distancePrim = dynamic_cast<const GEO_PrimVDB *>(prim);
			if (distancePrim->hasGrid()) {

				distance_baseGrid = distancePrim->getConstGridPtr();
				distance_grid = openvdb::gridConstPtrCast<openvdb::FloatGrid>(distance_baseGrid);
				if (distance_grid) {
					numberOfFoundVdbs += 1;
					break;
				}
			}
		}
	}
	// Make sure we got a valid prim
	if ((!vdbPrim) || (!distancePrim))
	{
		printf("react Number of Found Vdbs %i\n", numberOfFoundVdbs);
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error();
	}



	// Try to get the vdbs grid

	if (!grid) {
		printf("react Number of Found Vdbs %i\n", numberOfFoundVdbs);
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error();
	}
	if (!distance_grid) {
		printf("react Number of Found Vdbs %i\n", numberOfFoundVdbs);
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error();
	}
	openvdb::Vec3SGrid::Ptr cb = grid->deepCopy();


	struct Convolve {
		openvdb::FloatGrid::ConstPtr distance_grid;
		openvdb::Vec3SGrid::ConstPtr cb;
		float feed;
		float kill;
		float delta;
		float diffrate;
		openvdb::Coord ijk[12] = { openvdb::Coord(-1,-1,0),
			openvdb::Coord(-1,+1,0),
			openvdb::Coord(+1,-1,0),
			openvdb::Coord(+1,+1,0),
			openvdb::Coord(-1,0,-1),
			openvdb::Coord(-1,0,+1),
			openvdb::Coord(+1,0,-1),
			openvdb::Coord(+1,0,+1),
			openvdb::Coord(0,-1,-1),
			openvdb::Coord(0,-1,+1),
			openvdb::Coord(0,+1,-1),
			openvdb::Coord(0,+1,+1) };

		Convolve(openvdb::FloatGrid::ConstPtr g, openvdb::Vec3SGrid::Ptr cb, float f, float k, float d, float df) :
			distance_grid(g), feed(f), kill(k), delta(d), diffrate(df), cb(cb) {

		}
		inline void operator()(const openvdb::Vec3SGrid::ValueOnIter iter) const {
			openvdb::Vec3SGrid::ConstAccessor color_accessor = cb->getConstAccessor();
			openvdb::Vec3f lap = openvdb::Vec3f(0.0f, 0.0f, 0.0f);
			openvdb::Vec3f uv = color_accessor.getValue(iter.getCoord());

			for (int i = 0; i < 12; i++) {
				lap += color_accessor.getValue(iter.getCoord() + ijk[i]);
			}
			lap -= (12 * uv);
			float lapx = (lap.x()<0.0f) ? 0.0f : lap.x();
			float lapy = (lap.y()<0.0f) ? 0.0f : lap.y();
			//lap /= 4.0f;
			float du = diffrate * lapx - uv.x()*uv.y()*uv.y() + feed*(1.0 - uv.x());
			float dv = diffrate * lapy + uv.x()*uv.y()*uv.y() - (feed + kill)*uv.y();
			float u = uv.x() + delta * du;
			float v = uv.y() + delta * dv;
			if (u<0.0f) u = 0.0f;
			if (v<0.0f) v = 0.0f;
			if (u>1.0f) u = 1.0f;
			if (v>1.0f) v = 1.0f;
			iter.setValue(openvdb::Vec3f(u, v, 0.0f));
		}
	};

	float feed = FEED(context.getTime());
	float kill = KILL(context.getTime());
	float delta = DELTA(context.getTime());
	float diffrate = DIFFRATE(context.getTime());
	// Iterate over all active values.
	openvdb::tools::foreach(grid->beginValueOn(), Convolve(distance_grid, cb, feed, kill, delta, diffrate), true, true);
	cb->clear();
	cb->pruneGrid();
	cb = 0;
	


	return error();
}
