#pragma once
//#include <SOP_NodeVDB.h>
#include <openvdb/Grid.h>
#include <openvdb/tools/GridOperators.h>

#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
using VelocityAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Velocity_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
using GradientAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Gradient_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;

class Diverge {
private:
	const openvdb::math::Transform targetTransform;
	openvdb::Vec3SGrid::ConstPtr velocity_grid;
	openvdb::Vec3SGrid::ConstPtr gradient_grid;
public:
	Diverge(const openvdb::math::Transform t,
		openvdb::Vec3SGrid::ConstPtr vel_g,
		openvdb::Vec3SGrid::ConstPtr grad_g) :targetTransform(t), velocity_grid(vel_g), gradient_grid(grad_g)
	{};
	inline void operator()(
		const openvdb::Vec3SGrid::ValueOnCIter& iter,
		openvdb::FloatGrid::Accessor& accessor) const
	{
		openvdb::Vec3f temp = openvdb::Vec3f(0.0f, 0.0f, 0.0f);
		std::unique_ptr<VelocityAccessor> velocityAccessor;
		std::unique_ptr<Velocity_fastSampler> velocity_fastSampler;
		std::unique_ptr<GradientAccessor> gradientAccessor;
		std::unique_ptr<Gradient_fastSampler> gradient_fastSampler;
		velocityAccessor.reset(new VelocityAccessor(velocity_grid->getConstAccessor()));
		gradientAccessor.reset(new GradientAccessor(gradient_grid->getConstAccessor()));
		gradient_fastSampler.reset(new Gradient_fastSampler(*gradientAccessor, gradient_grid->transform()));
		float weight = 0;
		openvdb::Vec3f normal = gradient_fastSampler->wsSample(targetTransform.indexToWorld(iter.getCoord()));
		normal.normalize();
		openvdb::Vec3f surface_velocity_0;
		openvdb::Vec3f surface_velocity_1;
		//divergence
		surface_velocity_0 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(1, 0, 0))));
		surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
		surface_velocity_1 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(-1, 0, 0))));
		surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
		float dudx = surface_velocity_0.x() - surface_velocity_1.x();
		surface_velocity_0 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(0, 1, 0))));
		surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
		surface_velocity_1 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(0, -1, 0))));
		surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
		float dvdy = surface_velocity_0.y() - surface_velocity_1.y();
		surface_velocity_0 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(0, 0, 1))));
		surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
		surface_velocity_1 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(0, 0, -1))));
		surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
		float dwdz = surface_velocity_0.z() - surface_velocity_1.z();
		float divergence = (dudx + dvdy + dwdz) / (2.0f * targetTransform.voxelSize().x());
		accessor.setValue(iter.getCoord(), divergence);
	}


};