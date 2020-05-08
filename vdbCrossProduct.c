
#include "vdbCrossProduct.h"



using namespace VdbCappucino;
namespace hvdb = openvdb_houdini;
namespace hutil = houdini_utils;


// label node inputs, 0 corresponds to first input, 1 to the second one
const char *
SOP_VdbCrossProduct::inputLabel(unsigned idx) const
{
	switch (idx) {
	case 0: return "vectorfield A";
	case 1: return "vectorfield B";
	default: return "default";
	}
};




// constructors, destructors, usually there is no need to really modify anything here, the constructor's job is to ensure the node is put into the proper network
OP_Node *
SOP_VdbCrossProduct::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
	return new SOP_VdbCrossProduct(net, name, op);
}

SOP_VdbCrossProduct::SOP_VdbCrossProduct(OP_Network *net, const char *name, OP_Operator *op) : openvdb_houdini::SOP_NodeVDB(net, name, op) {}

SOP_VdbCrossProduct::~SOP_VdbCrossProduct() {}


// function that does the actual job
OP_ERROR
SOP_VdbCrossProduct::cookMySop(OP_Context &context)
{
	try {
		hutil::ScopedInputLock lock(*this, context);
		duplicateSourceStealable(0, context);

		const fpreal time = context.getTime();

		hvdb::Interrupter boss("CrossProduct");


		UT_String aGroupStr;
		evalString(aGroupStr, "agroup", 0, time);
		const GA_PrimitiveGroup* aGroup = matchGroup(*gdp, aGroupStr.toStdString());


		const GU_Detail* bGdp = inputGeo(1, context);
	

		UT_String bGroupStr;
		evalString(bGroupStr, "bgroup", 0, time);
		const GA_PrimitiveGroup* bGroup = matchGroup(const_cast<GU_Detail&>(*bGdp), bGroupStr.toStdString());


		bool processedVDB = false;
		//get b
		hvdb::VdbPrimCIterator gIt(bGdp, bGroup);
		const GU_PrimVDB *bPrim = *gIt;

		if (!bPrim) {
			addError(SOP_MESSAGE, "Missing b grid");
			return error();
		}
		if (bPrim->getStorageType() != UT_VDB_VEC3F) {
			addError(SOP_MESSAGE, "Expected b grid to be of type Vec3f");
			return error();
		}
		openvdb::GridBase::ConstPtr b_baseGrid;
		openvdb::Vec3SGrid::ConstPtr b_grid;
		b_baseGrid = bPrim->getConstGridPtr();
		b_grid = openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(b_baseGrid);
		if (!b_grid) {
			addError(SOP_MESSAGE, "Missing b grid");
			return error();
		}



		//process
		for (hvdb::VdbPrimIterator vdbIt(gdp, aGroup); vdbIt; ++vdbIt) {

			if (boss.wasInterrupted()) break;

			if (vdbIt->getGrid().type() == openvdb::Vec3fGrid::gridType()) {

				processedVDB = true;

				vdbIt->makeGridUnique();

				//openvdb::Vec3fGrid& grid =
				openvdb::Vec3fGrid::Ptr a_grid = openvdb::gridPtrCast<openvdb::Vec3fGrid>(vdbIt->getGridPtr());

				openvdb::Vec3fGrid::Ptr targetGrid = a_grid;// a_grid->deepCopy();

				std::string gridName = a_grid->getName();
				// Iterate over all active values.

				openvdb::tools::foreach(a_grid->beginValueOn(), CrossProduct(b_grid, a_grid->transform()), true, 0);

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

