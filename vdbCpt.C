#include "vdbCpt.h"
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
SOP_VdbCpt::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "Input VDB";
	case 1: return "Closest Point Field";
	case 2: return "Distance VDB ";
	default: return "default";
	}
}

// define parameter for debug option
static PRM_Name debugPRM("debug", "Print debug information"); // internal name, UI name
static PRM_Name doworldPRM("doworldpos", "create Vector to World Pos"); 
static PRM_Name maxcellsPRM("maxcells", "max cells value");
static PRM_Default maxcellsDefault(3.33);
															  // assign parameter to the interface, which is array of PRM_Template objects
PRM_Template SOP_VdbCpt::myTemplateList[] =
{
	PRM_Template(PRM_TOGGLE, 1, &debugPRM, PRMzeroDefaults), // type (checkbox), size (one in our case, but rgb/xyz values would need 3), pointer to a PRM_Name describing the parameter name, default value (0 - disabled)
	PRM_Template(PRM_FLT, 1, &maxcellsPRM, &maxcellsDefault),
	PRM_Template(PRM_TOGGLE, 1, &doworldPRM, PRMzeroDefaults),
	PRM_Template() // at the end there needs to be one empty PRM_Template object
};

// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbCpt::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbCpt(net, name, op);
}

SOP_VdbCpt::SOP_VdbCpt(OP_Network *net, const char *name, OP_Operator *op) : SOP_Node(net, name, op) {}

SOP_VdbCpt::~SOP_VdbCpt() {}

// function that does the actual job
OP_ERROR
SOP_VdbCpt::cookMySop(OP_Context &context)
{
	// we must lock our inputs before we try to access their geometry, OP_AutoLockInputs will automatically unlock our inputs when we return
	OP_AutoLockInputs inputs(this);
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();
	
	// duplicate our incoming geometry
	duplicateSource(0, context);
	float maxCells = MAXCELLS(context.getTime());
	// check for interrupt - interrupt scope closes automatically when 'progress' is destructed.
	UT_AutoInterrupt progress("Activating voxels...");

	// get pointer to geometry from second input
	const GU_Detail *cpm_gdp = inputGeo(1);
	const GU_Detail *dist_gdp = inputGeo(2);
	// Get the first VDB primitive in the geometry

	int numberOfFoundVdbs = 0;
	GEO_PrimVDB* vdbPrim = NULL;
	const GEO_PrimVDB* cpmPrim = NULL;
	const GEO_PrimVDB* distPrim = NULL;
	openvdb::GridBase::Ptr color_baseGrid;
	openvdb::GridBase::ConstPtr cpm_baseGrid;
	openvdb::GridBase::ConstPtr dist_baseGrid;
	openvdb::Vec3SGrid::Ptr grid;
	openvdb::Vec3SGrid::Ptr grid_old;
	openvdb::Vec3SGrid::ConstPtr cpm_grid;
	openvdb::FloatGrid::ConstPtr dist_grid;

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

	for (GA_Iterator it(gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		GEO_Primitive* prim = gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<GEO_PrimVDB *>(prim))
		{
			vdbPrim = dynamic_cast<GEO_PrimVDB *>(prim);
			if (vdbPrim->hasGrid()) {
				vdbPrim->makeGridUnique();
				color_baseGrid = vdbPrim->getGridPtr();
				grid_old = openvdb::gridPtrCast<openvdb::Vec3SGrid>(color_baseGrid);
				if (grid_old && (grid_old->getName() == "Cd_old")) {
					numberOfFoundVdbs += 1;

					break;
				}
			}
		}
	}
	for (GA_Iterator it(cpm_gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		const GEO_Primitive* prim = cpm_gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<const GEO_PrimVDB *>(prim))
		{
			cpmPrim = dynamic_cast<const GEO_PrimVDB *>(prim);
			if (cpmPrim->hasGrid()) {
				//vdbPrim->makeGridUnique();
				cpm_baseGrid = cpmPrim->getConstGridPtr();
				cpm_grid = openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(cpm_baseGrid);
				if (cpm_grid) {
					numberOfFoundVdbs += 1;
					break;
				}
			}
		}
	}
	for (GA_Iterator it(dist_gdp->getPrimitiveRange()); !it.atEnd(); it.advance())
	{
		const GEO_Primitive* prim = dist_gdp->getGEOPrimitive(it.getOffset());
		if (dynamic_cast<const GEO_PrimVDB *>(prim))
		{
			distPrim = dynamic_cast<const GEO_PrimVDB *>(prim);
			if (distPrim->hasGrid()) {
				//vdbPrim->makeGridUnique();
				dist_baseGrid = distPrim->getConstGridPtr();
				dist_grid = openvdb::gridConstPtrCast<openvdb::FloatGrid>(dist_baseGrid);
				if (dist_grid) {
					numberOfFoundVdbs += 1;
					break;
				}
			}
		}
	}

	// Make sure we got a valid prim
	if ((!vdbPrim) || (!cpmPrim) || (!distPrim))
	{
		printf("cpt Number of Found Vdbs %i\n", numberOfFoundVdbs);
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error();
	}

	// Docs say to do this in case the grid is shared



	// Try to get the vdbs grid


	if ((!grid) || (!cpm_grid) || (!grid_old) || (!dist_grid)){
		addError(SOP_MESSAGE, "Input geometry must contain a VDB");
		return error();
	}
	//openvdb::Vec3SGrid::Ptr grid_buffer = grid->copy();
	//openvdb::Vec3SGrid::Ptr grid_buffer_old = grid_old->copy();



	class ClosestPointOp {
	private:
		openvdb::Vec3SGrid::Ptr grid;
		openvdb::Vec3SGrid::Ptr color_grid;
		openvdb::Vec3SGrid::ConstPtr cpm_grid;
		openvdb::FloatGrid::ConstPtr dist_grid;
		using CpmAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
		using ColorAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
		using DistAccessor = typename openvdb::FloatGrid::ConstAccessor;
		using Cpm_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
		using Color_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
		using Dist_fastSampler = openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler>;
		float maxCells;
		int doworld;
	public:
		ClosestPointOp(openvdb::Vec3SGrid::Ptr  g,
			openvdb::Vec3SGrid::Ptr color_g,
			openvdb::Vec3SGrid::ConstPtr cpm_g,
			openvdb::FloatGrid::ConstPtr dist_g,
			int dw,
			float mC) :grid(g), maxCells(mC), color_grid(color_g), cpm_grid(cpm_g), dist_grid(dist_g), doworld(dw) {
		};
	

		inline void operator()(const openvdb::Vec3SGrid::ValueOnIter iter) const {
			std::unique_ptr<DistAccessor> distAccessor;
			std::unique_ptr<Dist_fastSampler> dist_fastSampler;
			std::unique_ptr<CpmAccessor> cpmAccessor;
			std::unique_ptr<Cpm_fastSampler> cpm_fastSampler; 
			std::unique_ptr<ColorAccessor> colorAccessor;
			std::unique_ptr<Color_fastSampler> color_fastSampler;
			distAccessor.reset(new DistAccessor(dist_grid->getConstAccessor()));
			dist_fastSampler.reset(new Dist_fastSampler(*distAccessor, dist_grid->transform()));
			cpmAccessor.reset(new CpmAccessor(cpm_grid->getConstAccessor()));
			cpm_fastSampler.reset(new Cpm_fastSampler(*cpmAccessor, cpm_grid->transform()));
			colorAccessor.reset(new ColorAccessor(color_grid->getConstAccessor()));
			color_fastSampler.reset(new Color_fastSampler(*colorAccessor, color_grid->transform()));

			
			
			
			const float distance = fabs(dist_fastSampler->wsSample(grid->transform().indexToWorld(iter.getCoord())));
			//const float distance = fabs(distAccessor->getValue(iter.getCoord()));
			if (distance>(maxCells*grid->voxelSize().x()))
			{
				iter.setValue(openvdb::Vec3f(0, 0, 0));
				iter.setValueOff();
			}
			else {
				const openvdb::Vec3f closestPoint = cpm_fastSampler->wsSample(grid->transform().indexToWorld(iter.getCoord()));
				//const openvdb::Vec3f closestPoint = cpmAccessor->getValue(iter.getCoord());
				//printf("closest Point %f, %f, %f\n",closestPoint.x(),closestPoint.y(),closestPoint.z());
				openvdb::Vec3f col;
				if (doworld == 0) {
					col = (color_fastSampler->wsSample(closestPoint));
				}
				else {
					col = closestPoint - grid->transform().indexToWorld(iter.getCoord()) ;
				}
				float x = col.x();
				float y = col.y();
				//if (x>1.0f) x = 1.0f;
				//if (x<-1.0f) x = -1.0f;
				//if (y>1.0f) y = 1.0f;
				//if (y<-1.0f) y = -1.0f;
				//if ((x<0.00001f)&&(x>-0.00001f)) x = 0.0f;
				//if ((y<0.00001f)&&(y>-0.00001f)) y = 0.0f;
				if ((x == 0.0f) && (y == 0.0f)) {
					iter.setValue(openvdb::Vec3f(0, 0, 0));
					iter.setValueOff();
				}
				else {
					iter.setValue(openvdb::Vec3f(x, y, col.z()));

				}
			}
		}
	};
	openvdb::tools::foreach(grid->beginValueOn(), ClosestPointOp(grid, grid, cpm_grid, dist_grid, DOWORLDCOORDS(), maxCells), true, true);
	//openvdb::tools::foreach(grid_old->beginValueOn(), ClosestPointOp(grid_old, grid_old, cpm_grid, dist_grid, maxCells), true, false);


	//grid_buffer->clear();
	//grid_buffer_old->clear();

	return error();
}
