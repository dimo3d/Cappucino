#include "vdbDivergence.h"
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
using ColorAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Color_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
using VelocityAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Velocity_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
using GradientAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Gradient_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;

// label node inputs, 0 corresponds to first input, 1 to the second one
const char *
SOP_VdbDivergence::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "VDB";
	case 1: return "velocity";
	case 2: return "gradient";
	default: return "default";
	}
}

// define parameter for debug option
static PRM_Name debugPRM("debug", "Print debug information"); // internal name, UI name
static PRM_Name dtPRM("dt", "delta time value");
static PRM_Default dtDefault(0.04);
															  // assign parameter to the interface, which is array of PRM_Template objects
PRM_Template SOP_VdbDivergence::myTemplateList[] =
{
	PRM_Template(PRM_TOGGLE, 1, &debugPRM, PRMzeroDefaults), // type (checkbox), size (one in our case, but rgb/xyz values would need 3), pointer to a PRM_Name describing the parameter name, default value (0 - disabled)
	PRM_Template(PRM_FLT, 1, &dtPRM, &dtDefault),
	PRM_Template() // at the end there needs to be one empty PRM_Template object
};

// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbDivergence::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbDivergence(net, name, op);
}

SOP_VdbDivergence::SOP_VdbDivergence(OP_Network *net, const char *name, OP_Operator *op) : SOP_Node(net, name, op) {}

SOP_VdbDivergence::~SOP_VdbDivergence() {}

// function that does the actual job
OP_ERROR
SOP_VdbDivergence::cookMySop(OP_Context &context)
{
	// we must lock our inputs before we try to access their geometry, OP_AutoLockInputs will automatically unlock our inputs when we return
	OP_AutoLockInputs inputs(this);
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();
	int numberOfFoundVdbs = 0;
	// duplicate our incoming geometry
	duplicateSource(0, context);
	float dt = DT();
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
	const GU_Detail *velocity_gdp = inputGeo(1);
	const GU_Detail *gradient_gdp = inputGeo(2);
	const GEO_PrimVDB* velocityPrim = NULL;
	const GEO_PrimVDB* gradientPrim = NULL;
	openvdb::GridBase::ConstPtr velocity_baseGrid;
	openvdb::Vec3SGrid::ConstPtr velocity_grid;
	openvdb::GridBase::ConstPtr gradient_baseGrid;
	openvdb::Vec3SGrid::ConstPtr gradient_grid;
	for (GA_Iterator it(velocity_gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		const GEO_Primitive* prim = velocity_gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<const GEO_PrimVDB *>(prim))
		{
			velocityPrim = dynamic_cast<const GEO_PrimVDB *>(prim);
			if (velocityPrim->hasGrid()) {

				velocity_baseGrid = velocityPrim->getConstGridPtr();
				velocity_grid = openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(velocity_baseGrid);
				if (velocity_grid) {
					
					break;
				}
			}
		}
	}
	for (GA_Iterator it(gradient_gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		const GEO_Primitive* prim = gradient_gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<const GEO_PrimVDB *>(prim))
		{
			gradientPrim = dynamic_cast<const GEO_PrimVDB *>(prim);
			if (gradientPrim->hasGrid()) {

				gradient_baseGrid = gradientPrim->getConstGridPtr();
				gradient_grid = openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(gradient_baseGrid);
				if (gradient_grid) {
					
					break;
				}
			}
		}
	}
	// Make sure we got a valid prim
	if ((!vdbPrim) || (!velocityPrim)||(!gradientPrim))
	{
		if (DEBUG())
		printf("Diverge Number of Found Vdbs");
		addError(SOP_MESSAGE, "Input geometry must contain all VDBs");
		return error();
	}
	if (!grid) {
		addError(SOP_MESSAGE, "Input geometry must contain a Color VDB");
		return error();
	};
	if (!velocity_grid) {
		addError(SOP_MESSAGE, "second input geometry must contain a velocity VDB");
		return error();
	};
	if (!gradient_grid) {
		addError(SOP_MESSAGE, "second input geometry must contain a gradient VDB");
		return error();
	};
	openvdb::Vec3SGrid::ConstPtr grid_buffer = grid->deepCopy();
	openvdb::Vec3SGrid::ConstAccessor color_accessor = grid_buffer->getConstAccessor();

	

	struct Diverge {
		openvdb::Vec3SGrid::ConstPtr color_grid;
		openvdb::Vec3SGrid::ConstPtr velocity_grid;
		openvdb::Vec3SGrid::ConstPtr gradient_grid;
		float dt;
		Diverge(openvdb::Vec3SGrid::ConstPtr g, 
				openvdb::Vec3SGrid::ConstPtr vel_g,
				openvdb::Vec3SGrid::ConstPtr grad_g,
			float d) :color_grid(g), velocity_grid(vel_g), gradient_grid(grad_g), dt(d)
				{}
	
		inline void operator()(const openvdb::Vec3SGrid::ValueOnIter iter) const {
			openvdb::Vec3f temp = openvdb::Vec3f(0.0f, 0.0f, 0.0f);
			std::unique_ptr<VelocityAccessor> velocityAccessor;
			std::unique_ptr<Velocity_fastSampler> velocity_fastSampler;
			std::unique_ptr<GradientAccessor> gradientAccessor;
			std::unique_ptr<Gradient_fastSampler> gradient_fastSampler;
			std::unique_ptr<ColorAccessor> colorAccessor;
			std::unique_ptr<Color_fastSampler> color_fastSampler;
			velocityAccessor.reset(new VelocityAccessor(velocity_grid->getConstAccessor()));
			
			gradientAccessor.reset(new GradientAccessor(gradient_grid->getConstAccessor()));
			gradient_fastSampler.reset(new Gradient_fastSampler(*gradientAccessor, gradient_grid->transform()));
			velocity_fastSampler.reset(new Velocity_fastSampler(*velocityAccessor, velocity_grid->transform()));
			colorAccessor.reset(new ColorAccessor(color_grid->getConstAccessor()));
			color_fastSampler.reset(new Color_fastSampler(*colorAccessor, color_grid->transform()));
			float weight = 0;
			openvdb::Vec3f normal = gradient_fastSampler->wsSample(color_grid->transform().indexToWorld(iter.getCoord()));
			normal.normalize();
			openvdb::Vec3f surface_velocity_0;// = velocity_fastSampler->wsSample(color_grid->transform().indexToWorld(iter.getCoord()));
			openvdb::Vec3f surface_velocity_1; //= velocity - velocity.projection(normal);
			//divergence
			surface_velocity_0 = (velocity_fastSampler->wsSample(color_grid->transform().indexToWorld(iter.getCoord() + openvdb::Coord(1, 0, 0))));
			surface_velocity_0  = surface_velocity_0 - surface_velocity_0.projection(normal);
			surface_velocity_1 = (velocity_fastSampler->wsSample(color_grid->transform().indexToWorld(iter.getCoord() + openvdb::Coord(-1, 0, 0))));
			surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
			float dudx = surface_velocity_0.x() - surface_velocity_1.x();
			surface_velocity_0 = (velocity_fastSampler->wsSample(color_grid->transform().indexToWorld(iter.getCoord() + openvdb::Coord(0, 1, 0))));
			surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
			surface_velocity_1 = (velocity_fastSampler->wsSample(color_grid->transform().indexToWorld(iter.getCoord() + openvdb::Coord(0, -1, 0))));
			surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
			float dvdy = surface_velocity_0.y() - surface_velocity_1.y();
			surface_velocity_0 = (velocity_fastSampler->wsSample(color_grid->transform().indexToWorld(iter.getCoord() + openvdb::Coord(0, 0, 1))));
			surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
			surface_velocity_1 = (velocity_fastSampler->wsSample(color_grid->transform().indexToWorld(iter.getCoord() + openvdb::Coord(0, 0, -1))));
			surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
			float dwdz = surface_velocity_0.z() - surface_velocity_1.z();
			float divergence = (dudx + dvdy + dwdz) / (2.0f * color_grid->voxelSize().x());
			iter.setValue(iter.getValue() + iter.getValue() * divergence * dt);
		}
	};


	// Iterate over all active values.
	openvdb::tools::foreach(grid->beginValueOn(), Diverge(grid,velocity_grid,gradient_grid, dt), false, 0);
	
	return error();
}
