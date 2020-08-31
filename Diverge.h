#pragma once
//#include <SOP_NodeVDB.h>
#include <openvdb/Grid.h>
#include <openvdb/tools/GridOperators.h>

#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
using VelocityAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Velocity_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>;
using GradientAccessor = typename openvdb::Vec3SGrid::ConstAccessor;
using Gradient_fastSampler = openvdb::tools::GridSampler<openvdb::Vec3SGrid::ConstAccessor, openvdb::tools::QuadraticSampler>;

struct CrossProduct {
	openvdb::Vec3SGrid::ConstPtr b_grid;
	openvdb::math::Transform trans;
	CrossProduct(
		openvdb::Vec3SGrid::ConstPtr b_g, openvdb::math::Transform tr
	) : b_grid(b_g), trans(tr)
	{}

	inline void operator()(const openvdb::Vec3SGrid::ValueOnIter iter) const {
		std::unique_ptr<GradientAccessor> bAccessor;
		std::unique_ptr<Gradient_fastSampler> b_fastSampler;

		bAccessor.reset(new GradientAccessor(b_grid->getConstAccessor()));
		b_fastSampler.reset(new Gradient_fastSampler(*bAccessor, b_grid->transform()));


		openvdb::Vec3f b_Vec = b_fastSampler->wsSample(trans.indexToWorld(iter.getCoord()));
		openvdb::Vec3f a_Vec = iter.getValue();
		
		iter.setValue(a_Vec.cross(b_Vec));
	}
};
struct ProjectVectorToSurface {
	openvdb::Vec3SGrid::ConstPtr grad_grid;
	openvdb::math::Transform trans;
	bool keepLength;
	ProjectVectorToSurface(
		openvdb::Vec3SGrid::ConstPtr grad_g, openvdb::math::Transform tr, bool kL = true
	) : grad_grid(grad_g), trans(tr), keepLength(kL)
	{}

	inline void operator()(const openvdb::Vec3SGrid::ValueOnIter iter) const {
		std::unique_ptr<GradientAccessor> gradientAccessor;
		std::unique_ptr<Gradient_fastSampler> gradient_fastSampler;

		gradientAccessor.reset(new GradientAccessor(grad_grid->getConstAccessor()));
		gradient_fastSampler.reset(new Gradient_fastSampler(*gradientAccessor, grad_grid->transform()));


		openvdb::Vec3f normal = gradient_fastSampler->wsSample(trans.indexToWorld(iter.getCoord()));
		openvdb::Vec3f oldVelocity = iter.getValue();
		float length = oldVelocity.length();
		normal.normalize();
		openvdb::Vec3f projected_Velocity = oldVelocity - oldVelocity.projection(normal);
		if (keepLength)
			{ 
			projected_Velocity.normalize();
			projected_Velocity *= length;

			}
		iter.setValue(projected_Velocity);
	}
};
struct applyJacobiMatrix {
	openvdb::Vec3SGrid::ConstPtr grad_grid;
	openvdb::Vec3SGrid::ConstPtr exvel_grid;
	
	applyJacobiMatrix(
		openvdb::Vec3SGrid::ConstPtr grad_g, openvdb::Vec3SGrid::ConstPtr exvel_g
	) : grad_grid(grad_g), exvel_grid(exvel_g)
	{}
	inline void operator()(const openvdb::Vec3SGrid::ValueOnIter iter) const {
		

		// Compute the value of the grid at ijk via nearest-neighbor (zero-order)
		// interpolation.
		const openvdb::Vec3R ijk= iter.getCoord().asVec3d();
		openvdb::Vec3R normal = openvdb::tools::PointSampler::sample(grad_grid->tree(), ijk);
		
		openvdb::Vec3R surface_velocity_0;
		openvdb::Vec3R surface_velocity_1;
		//divergence
		surface_velocity_0 = openvdb::tools::PointSampler::sample(exvel_grid->tree(), ijk + openvdb::Vec3R(1, 0, 0));
		//surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
		surface_velocity_1 = openvdb::tools::PointSampler::sample(exvel_grid->tree(), ijk + openvdb::Vec3R(-1, 0, 0));
		//surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
		double dudx = surface_velocity_0.x() - surface_velocity_1.x();
		double dvdx = surface_velocity_0.y() - surface_velocity_1.y();
		double dwdx = surface_velocity_0.z() - surface_velocity_1.z();
		surface_velocity_0 = openvdb::tools::PointSampler::sample(exvel_grid->tree(), ijk + openvdb::Vec3R(0, 1, 0));
		//surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
		surface_velocity_1 = openvdb::tools::PointSampler::sample(exvel_grid->tree(), ijk + openvdb::Vec3R(0, -1, 0));
		//surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
		double dudy = surface_velocity_0.x() - surface_velocity_1.x();
		double dvdy = surface_velocity_0.y() - surface_velocity_1.y();
		double dwdy = surface_velocity_0.z() - surface_velocity_1.z();
		surface_velocity_0 = openvdb::tools::PointSampler::sample(exvel_grid->tree(), ijk + openvdb::Vec3R(0, 0, 1));
		//surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
		surface_velocity_1 = openvdb::tools::PointSampler::sample(exvel_grid->tree(), ijk + openvdb::Vec3R(0, 0, -1));
		//surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
		double dudz = surface_velocity_0.x() - surface_velocity_1.x();
		double dvdz = surface_velocity_0.y() - surface_velocity_1.y();
		double dwdz = surface_velocity_0.z() - surface_velocity_1.z();
		//openvdb::math::Mat3<double> jacobi = openvdb::math::Mat3< double >(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz);
		
		openvdb::Vec3R external_velocity = openvdb::tools::PointSampler::sample(exvel_grid->tree(), ijk);
		openvdb::Vec3R oldVelocity = iter.getValue();
		openvdb::Vec3R newVelocity = openvdb::Vec3R(
			dudx*oldVelocity.x() + dudy*oldVelocity.y() + dudz*oldVelocity.z(),
			dvdx*oldVelocity.x() + dvdy*oldVelocity.y() + dvdz*oldVelocity.z(),
			dwdx*oldVelocity.x() + dwdy*oldVelocity.y() + dwdz*oldVelocity.z());
		iter.setValue(newVelocity);
	}
};
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
		double dudx = surface_velocity_0.x() - surface_velocity_1.x();
		surface_velocity_0 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(0, 1, 0))));
		surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
		surface_velocity_1 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(0, -1, 0))));
		surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
		double dvdy = surface_velocity_0.y() - surface_velocity_1.y();
		surface_velocity_0 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(0, 0, 1))));
		surface_velocity_0 = surface_velocity_0 - surface_velocity_0.projection(normal);
		surface_velocity_1 = (velocityAccessor->getValue((iter.getCoord() + openvdb::Coord(0, 0, -1))));
		surface_velocity_1 = surface_velocity_1 - surface_velocity_1.projection(normal);
		double dwdz = surface_velocity_0.z() - surface_velocity_1.z();
		double divergence = (dudx + dvdy + dwdz) / (2.0f * targetTransform.voxelSize().x());
		accessor.setValue(iter.getCoord(), divergence);
	}


};