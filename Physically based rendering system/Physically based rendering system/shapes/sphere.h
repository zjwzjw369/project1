#pragma once
#ifndef PBRS_SHAPES_SPHERE_H
#define PBRS_SHAPES_SPHERE_H

#include "../core/pbrs.h"
#include "../core/shape.h"


namespace pbrs {
	class Sphere :public Shape{
	public:
		Sphere(const Transform *ObjectToWorld, const Transform *WorldToObject,
			bool reverseOrientation, Float radius, Float zMin, Float zMax,
			Float phiMax): Shape(ObjectToWorld, WorldToObject, reverseOrientation),
			radius(radius), zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
			zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
			thetaMin(std::acos(Clamp(zMin / radius, -1, 1))),
			thetaMax(std::acos(Clamp(zMax / radius, -1, 1))),
			phiMax(Radians(Clamp(phiMax, 0, 360))) { }
		Bounds3f Sphere::ObjectBound() const {
			return Bounds3f(Point3f(-radius, -radius, zMin),
				Point3f(radius, radius, zMax));
		}
		bool Intersect(const Ray &r, Float *tHit,SurfaceInteraction *isect, bool testAlphaTexture) const;
	private:
		const Float radius;
		const Float zMin, zMax;
		const Float thetaMin, thetaMax, phiMax;
	};

}



















#endif // !PBRS_SHAPES_SPHERE_H
