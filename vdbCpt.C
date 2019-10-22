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
#include <Utils.h>
#include <ParmFactory.h>
using namespace VdbCappucino;
namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;
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



// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbCpt::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbCpt(net, name, op);
}

SOP_VdbCpt::SOP_VdbCpt(OP_Network *net, const char *name, OP_Operator *op) : openvdb_houdini::SOP_NodeVDB(net, name, op) {}

SOP_VdbCpt::~SOP_VdbCpt() {}

// function that does the actual job
OP_ERROR
SOP_VdbCpt::cookMySop(OP_Context &context)
{
	
	float maxCells = MAXCELLS(context.getTime());

	class ClosestPointOp {
	private:
		openvdb::Vec3SGrid::Ptr grid;
		openvdb::Vec3SGrid::Ptr color_grid;
		openvdb::Vec3SGrid::ConstPtr cpm_grid;
		openvdb::FloatGrid::ConstPtr dist_grid;
		using CpmAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
		using ColorAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
		using DistAccessor = typename openvdb::FloatGrid::ConstAccessor;
		using Cpm_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::PointSampler>;
		using Color_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::PointSampler>;
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
					iter.setValue(col);
				
			}
		}
	};

	try {
		hutil::ScopedInputLock lock(*this, context);
		duplicateSourceStealable(0, context);

		const fpreal time = context.getTime();

		hvdb::Interrupter boss("CPT ");


		UT_String GroupStr;
		evalString(GroupStr, "group", 0, time);
		const GA_PrimitiveGroup* Group = matchGroup(*gdp, GroupStr.toStdString());


		const GU_Detail* cptGdp = inputGeo(1, context);
		const GU_Detail* distGdp = inputGeo(2, context);

		UT_String cptGroupStr;
		evalString(cptGroupStr, "cptGroup", 0, time);
		const GA_PrimitiveGroup* cptGroup = matchGroup(const_cast<GU_Detail&>(*cptGdp), cptGroupStr.toStdString());

		UT_String distGroupStr;
		evalString(distGroupStr, "distanceGroup", 0, time);
		const GA_PrimitiveGroup* distGroup = matchGroup(const_cast<GU_Detail&>(*distGdp), distGroupStr.toStdString());

		bool processedVDB = false;
		//get gradient
		hvdb::VdbPrimCIterator gIt(cptGdp, cptGroup);
		const GU_PrimVDB *cptPrim = *gIt;

		if (!cptPrim) {
			addError(SOP_MESSAGE, "Missing cpt grid");
			return error();
		}
		if (cptPrim->getStorageType() != UT_VDB_VEC3F) {
			addError(SOP_MESSAGE, "Expected cpt grid to be of type Vec3f");
			return error();
		}

		//get gradient
		hvdb::VdbPrimCIterator dIt(distGdp, distGroup);
		const GU_PrimVDB *distPrim = *dIt;

		if (!distPrim) {
			addError(SOP_MESSAGE, "Missing distance grid");
			return error();
		}
		if (distPrim->getStorageType() != UT_VDB_FLOAT) {
			addError(SOP_MESSAGE, "Expected gradient grid to be of type Float");
			return error();
		}

		openvdb::GridBase::ConstPtr cpt_baseGrid;
		openvdb::Vec3SGrid::ConstPtr cpt_grid;
		cpt_baseGrid = cptPrim->getConstGridPtr();
		cpt_grid = openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(cpt_baseGrid);
		if (!cpt_grid) {
			addError(SOP_MESSAGE, "Missing cpt grid");
			return error();
		}

		openvdb::GridBase::ConstPtr dist_baseGrid;
		openvdb::FloatGrid::ConstPtr dist_grid;
		dist_baseGrid = distPrim->getConstGridPtr();
		dist_grid = openvdb::gridConstPtrCast<openvdb::FloatGrid>(dist_baseGrid);
		if (!dist_grid) {
			addError(SOP_MESSAGE, "Missing distance grid");
			return error();
		}

		//process
		for (hvdb::VdbPrimIterator vdbIt(gdp, Group); vdbIt; ++vdbIt) {

			if (boss.wasInterrupted()) break;

			if (vdbIt->getGrid().type() == openvdb::Vec3fGrid::gridType()) {

				processedVDB = true;

				vdbIt->makeGridUnique();

				//openvdb::Vec3fGrid& grid =
				openvdb::Vec3fGrid::Ptr grid = openvdb::gridPtrCast<openvdb::Vec3fGrid>(vdbIt->getGridPtr());
				openvdb::Vec3fGrid::Ptr grid_copy = grid->deepCopy();

				//openvdb::tools::foreach(grid->beginValueOn(), Diverge(grid->transform(),velocity_grid,gradient_grid, dt), false, 0);
				openvdb::tools::foreach(grid->beginValueOn(), ClosestPointOp(grid, grid_copy, cpt_grid, dist_grid, DOWORLDCOORDS(), maxCells), true, false);
				grid_copy->clear();
				grid->pruneGrid();
			}
		}

		if (!processedVDB && !boss.wasInterrupted()) {
			addWarning(SOP_MESSAGE, "No Vec3f VDBs found.");
		}
	}
	catch (std::exception& e) {
		addError(SOP_MESSAGE, e.what());
	}

	return error();




}
