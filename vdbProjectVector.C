#include "vdbProjectVector.h"
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
#include <Utils.h>
#include <openvdb/openvdb.h>
#include <ParmFactory.h>
using namespace VdbCappucino;
namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;

// label node inputs, 0 corresponds to first input, 1 to the second one
const char *
SOP_VdbProjectVector::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "velocity";
	case 1: return "gradient of distance";
	default: return "default";
	}
}


// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbProjectVector::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbProjectVector(net, name, op);
}

SOP_VdbProjectVector::SOP_VdbProjectVector(OP_Network *net, const char *name, OP_Operator *op) : openvdb_houdini::SOP_NodeVDB(net, name, op) {}

SOP_VdbProjectVector::~SOP_VdbProjectVector() {}

// function that does the actual job
OP_ERROR
SOP_VdbProjectVector::cookMySop(OP_Context &context)
{
	try {
	hutil::ScopedInputLock lock(*this, context);
	duplicateSourceStealable(0, context);

	const fpreal time = context.getTime();

	hvdb::Interrupter boss("Project Vector");

	
	UT_String velocityGroupStr;
	evalString(velocityGroupStr, "velocitygroup", 0, time);
	const GA_PrimitiveGroup* velocityGroup = matchGroup(*gdp, velocityGroupStr.toStdString());


	const GU_Detail* gradientGdp = inputGeo(1, context);
	const GU_Detail* divergenceGdp = inputGeo(2, context);

	UT_String gradientGroupStr;
	evalString(gradientGroupStr, "gradientgroup", 0, time);
	const GA_PrimitiveGroup* gradientGroup = matchGroup(const_cast<GU_Detail&>(*gradientGdp), gradientGroupStr.toStdString());

	
	bool processedVDB = false;
	//get gradient
	hvdb::VdbPrimCIterator gIt(gradientGdp, gradientGroup);
	const GU_PrimVDB *gradientPrim = *gIt;

	if (!gradientPrim) {
		addError(SOP_MESSAGE, "Missing gradient grid");
		return error();
	}
	if (gradientPrim->getStorageType() != UT_VDB_VEC3F) {
		addError(SOP_MESSAGE, "Expected gradient grid to be of type Vec3f");
		return error();
	}
	openvdb::GridBase::ConstPtr gradient_baseGrid;
	openvdb::Vec3SGrid::ConstPtr gradient_grid;
	gradient_baseGrid = gradientPrim->getConstGridPtr();
	gradient_grid = openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(gradient_baseGrid);
	if (!gradient_grid) {
		addError(SOP_MESSAGE, "Missing gradient grid");
		return error();
	}



	//process
	for (hvdb::VdbPrimIterator vdbIt(gdp, velocityGroup); vdbIt; ++vdbIt) {

		if (boss.wasInterrupted()) break;

		if (vdbIt->getGrid().type() == openvdb::Vec3fGrid::gridType()) {

			processedVDB = true;

			vdbIt->makeGridUnique();

			//openvdb::Vec3fGrid& grid =
			openvdb::Vec3fGrid::Ptr velocity_grid = openvdb::gridPtrCast<openvdb::Vec3fGrid>(vdbIt->getGridPtr());

			openvdb::Vec3fGrid::Ptr targetGrid = velocity_grid;// velocity_grid->deepCopy();

			std::string gridName = velocity_grid->getName();
			// Iterate over all active values.
	
			openvdb::tools::foreach(velocity_grid->beginValueOn(), ProjectVectorToSurface(gradient_grid, velocity_grid->transform()), true, 0);
	
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

