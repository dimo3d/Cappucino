#pragma once
#include <SOP/SOP_Node.h>
#include <openvdb/tools/ValueTransformer.h>
namespace VdbCappucino {
	class SOP_VdbWave : public SOP_Node
	{
	public:
		// node contructor for HDK
		static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

		// parameter array for Houdini UI
		static PRM_Template myTemplateList[];

	protected:
		// constructor, destructor
		SOP_VdbWave(OP_Network *net, const char *name, OP_Operator *op);

		virtual ~SOP_VdbWave();

		// labeling node inputs in Houdini UI
		virtual const char *inputLabel(unsigned idx) const;

		// main function that does geometry processing
		virtual OP_ERROR cookMySop(OP_Context &context);

	private:
		// helper function for returning value of parameter
		int DEBUG() { return evalInt("debug", 0, 0); }
		
		fpreal DIM(int t) { return evalInt("dim", 0, t); }
		fpreal DT(fpreal t) { return evalFloat("dt", 0, t); }
		fpreal GRAVITY(fpreal t) { return evalFloat("gravity", 0, t); }
		fpreal ALPHA(fpreal t) { return evalFloat("alpha", 0, t); }
	};


}
