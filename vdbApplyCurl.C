#include "vdbApplyCurl.h"
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
SOP_VdbApplyCurl::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "velocity";
	case 1: return "external velocity";
	case 2: return "surface gradient";
	default: return "default";
	}
}


// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbApplyCurl::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbApplyCurl(net, name, op);
}

SOP_VdbApplyCurl::SOP_VdbApplyCurl(OP_Network *net, const char *name, OP_Operator *op) : openvdb_houdini::SOP_NodeVDB(net, name, op) {}

SOP_VdbApplyCurl::~SOP_VdbApplyCurl() {}

// function that does the actual job
OP_ERROR
SOP_VdbApplyCurl::cookMySop(OP_Context &context)
{
	try {
	hutil::ScopedInputLock lock(*this, context);
	duplicateSourceStealable(0, context);

	const fpreal time = context.getTime();

	hvdb::Interrupter boss("Apply Curl");

	
	UT_String velocityGroupStr;
	UT_String exvelGroupStr;
	UT_String gradientGroupStr;
	evalString(velocityGroupStr, "velocitygroup", 0, time);
	evalString(exvelGroupStr, "exvelgroup", 0, time);
	evalString(gradientGroupStr, "gradientgroup", 0, time);

	const GU_Detail* velocityGdp = inputGeo(0, context);
	const GU_Detail* exvelGdp = inputGeo(1, context);
	const GU_Detail* gradientGdp = inputGeo(2, context);

	
	const GA_PrimitiveGroup* velocityGroup = matchGroup(*gdp, velocityGroupStr.toStdString());
	const GA_PrimitiveGroup* exvelGroup = matchGroup(const_cast<GU_Detail&>(*exvelGdp), exvelGroupStr.toStdString());
	const GA_PrimitiveGroup* gradientGroup = matchGroup(const_cast<GU_Detail&>(*gradientGdp), gradientGroupStr.toStdString());
	
	bool processedVDB = false;
	//get exvel
	hvdb::VdbPrimCIterator gIt(exvelGdp, exvelGroup);
	const GU_PrimVDB *exvelPrim = *gIt;

	if (!exvelPrim) {
		addError(SOP_MESSAGE, "Missing exvel grid");
		return error();
	}
	if (exvelPrim->getStorageType() != UT_VDB_VEC3F) {
		addError(SOP_MESSAGE, "Expected exvel grid to be of type Vec3f");
		return error();
	}
	openvdb::GridBase::ConstPtr exvel_baseGrid;
	openvdb::Vec3SGrid::ConstPtr exvel_grid;
	exvel_baseGrid = exvelPrim->getConstGridPtr();
	exvel_grid = openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(exvel_baseGrid);
	if (!exvel_grid) {
		addError(SOP_MESSAGE, "Missing exvel grid");
		return error();
	}

	//get gradient
	hvdb::VdbPrimCIterator ggIt(gradientGdp, gradientGroup);
	const GU_PrimVDB *gradientPrim = *ggIt;

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
			
			// transform Grids
			openvdb::Vec3SGrid::Ptr
				transformed_gradient_grid = openvdb::Vec3SGrid::create();
			openvdb::Vec3SGrid::Ptr
				transformed_exvel_grid = openvdb::Vec3SGrid::create();
		
			// Get the source and target grids' index space to world space transforms.
			const openvdb::math::Transform
				&targetXform = velocity_grid->transform(),
				&source_gradient_Xform = gradient_grid->transform(),
				&source_exvel_Xform = exvel_grid->transform();
			transformed_gradient_grid->setTransform(velocity_grid->transformPtr());
			transformed_exvel_grid->setTransform(velocity_grid->transformPtr());
			// Compute a source grid to target grid transform.
			// (For this example, we assume that both grids' transforms are linear,
			// so that they can be represented as 4 x 4 matrices.)
			openvdb::Mat4R xform_gradient =
				source_gradient_Xform.baseMap()->getAffineMap()->getMat4() *
				targetXform.baseMap()->getAffineMap()->getMat4().inverse();
			openvdb::Mat4R xform_exvel =
				source_exvel_Xform.baseMap()->getAffineMap()->getMat4() *
				targetXform.baseMap()->getAffineMap()->getMat4().inverse();
			// Create the transformer.
			openvdb::tools::GridTransformer transformer_gradient(xform_gradient);
			openvdb::tools::GridTransformer transformer_exvel(xform_exvel);
			
			// Resample using trilinear interpolation.
			transformer_gradient.transformGrid<openvdb::tools::BoxSampler, openvdb::Vec3SGrid>(
				*gradient_grid, *transformed_gradient_grid);
			transformer_exvel.transformGrid<openvdb::tools::BoxSampler, openvdb::Vec3SGrid>(
				*exvel_grid, *transformed_exvel_grid);
					
			
			// Iterate over all active values.
			openvdb::tools::foreach(velocity_grid->beginValueOn(), applyJacobiMatrix( transformed_gradient_grid, transformed_exvel_grid), true, 0);
	
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

