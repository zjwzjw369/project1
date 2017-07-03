#pragma once

#ifndef PBRS_CORE_SHAPE_H
#define PBRS_CORE_SHAPE_H

#include "pbrs.h"
#include "transform.h"
namespace pbrs {
	class Shape {
	public:
		//³éÏóÀà
		Shape(const Transform *ObjectToWorld,
			const Transform *WorldToObject, bool reverseOrientation)
			: ObjectToWorld(ObjectToWorld), WorldToObject(WorldToObject),
			reverseOrientation(reverseOrientation),
			transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {
		}
		virtual bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture = true) const = 0;
		virtual Bounds3f ObjectBound() const = 0;
		virtual bool IntersectP(const Ray &ray, bool testAlphaTexture = true)const;
		virtual Float Area() const = 0;
		Bounds3f WorldBound() const {
			return (*ObjectToWorld)(ObjectBound());
		}
		const Transform *ObjectToWorld, *WorldToObject;
		const bool reverseOrientation;
		const bool transformSwapsHandedness;
	};
}







































#endif // !PBRS_CORE_SHAPE_H
