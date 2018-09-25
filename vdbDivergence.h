#pragma once
#include <SOP/SOP_Node.h>
#include <SOP_NodeVDB.h>
#include <openvdb/Grid.h>
#include <openvdb/tools/GridOperators.h>

#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
#include <Diverge.h>
namespace VdbCappucino {
	

	

	class SOP_VdbDivergence : public openvdb_houdini::SOP_NodeVDB
	{
	public:
		// node contructor for HDK
		static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

		// parameter array for Houdini UI
		static PRM_Template myTemplateList[];

	protected:
		// constructor, destructor
		SOP_VdbDivergence(OP_Network *net, const char *name, OP_Operator *op);

		virtual ~SOP_VdbDivergence();

		// labeling node inputs in Houdini UI
		virtual const char *inputLabel(unsigned idx) const;

		// main function that does geometry processing
		virtual OP_ERROR cookMySop(OP_Context &context);

	private:
		// helper function for returning value of parameter
		int DEBUG() { return evalInt("debug", 0, 0); }
		float DT() { return evalFloat("dt", 0, 0); }

	};
	

	
}
