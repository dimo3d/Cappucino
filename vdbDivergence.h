#pragma once
#include <SOP/SOP_Node.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/Grid.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
namespace VdbCappucino {
	
	

	class SOP_VdbDivergence : public SOP_Node
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
	
	using VelocityAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
	using Velocity_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
	using GradientAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
	using Gradient_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
	struct Diverge {
		//openvdb::Vec3SGrid::ConstPtr color_grid;
		const openvdb::math::Transform targetTransform;
		openvdb::Vec3SGrid::ConstPtr velocity_grid;
		openvdb::Vec3SGrid::ConstPtr gradient_grid;
		
		Diverge(const openvdb::math::Transform t,
			openvdb::Vec3SGrid::ConstPtr vel_g,
			openvdb::Vec3SGrid::ConstPtr grad_g,
			float d) :targetTransform(t), velocity_grid(vel_g), gradient_grid(grad_g)
		{}

		inline void operator()(const openvdb::FloatGrid::ValueOnIter iter) const {
			openvdb::Vec3f temp = openvdb::Vec3f(0.0f, 0.0f, 0.0f);
			std::unique_ptr<VelocityAccessor> velocityAccessor;
			std::unique_ptr<Velocity_fastSampler> velocity_fastSampler;
			std::unique_ptr<GradientAccessor> gradientAccessor;
			std::unique_ptr<Gradient_fastSampler> gradient_fastSampler;
			velocityAccessor.reset(new VelocityAccessor(velocity_grid->getConstAccessor()));
			gradientAccessor.reset(new GradientAccessor(gradient_grid->getConstAccessor()));
			gradient_fastSampler.reset(new Gradient_fastSampler(*gradientAccessor, gradient_grid->transform()));
			velocity_fastSampler.reset(new Velocity_fastSampler(*velocityAccessor, velocity_grid->transform()));
			float weight = 0;
			openvdb::Vec3f normal = gradient_fastSampler->wsSample(targetTransform.indexToWorld(iter.getCoord()));
			normal.normalize();
			openvdb::Vec3f surface_velocity_0;
			openvdb::Vec3f surface_velocity_1;
			//divergence
			surface_velocity_0 = (velocity_fastSampler->wsSample(targetTransform.indexToWorld(iter.getCoord() + openvdb::Coord(1, 0, 0))));
			surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
			surface_velocity_1 = (velocity_fastSampler->wsSample(targetTransform.indexToWorld(iter.getCoord() + openvdb::Coord(-1, 0, 0))));
			surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
			float dudx = surface_velocity_0.x() - surface_velocity_1.x();
			surface_velocity_0 = (velocity_fastSampler->wsSample(targetTransform.indexToWorld(iter.getCoord() + openvdb::Coord(0, 1, 0))));
			surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
			surface_velocity_1 = (velocity_fastSampler->wsSample(targetTransform.indexToWorld(iter.getCoord() + openvdb::Coord(0, -1, 0))));
			surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
			float dvdy = surface_velocity_0.y() - surface_velocity_1.y();
			surface_velocity_0 = (velocity_fastSampler->wsSample(targetTransform.indexToWorld(iter.getCoord() + openvdb::Coord(0, 0, 1))));
			surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
			surface_velocity_1 = (velocity_fastSampler->wsSample(targetTransform.indexToWorld(iter.getCoord() + openvdb::Coord(0, 0, -1))));
			surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
			float dwdz = surface_velocity_0.z() - surface_velocity_1.z();
			float divergence = -(dudx + dvdy + dwdz) / (2.0f * targetTransform.voxelSize().x());

			iter.setValue( divergence);
		}
	};
}
