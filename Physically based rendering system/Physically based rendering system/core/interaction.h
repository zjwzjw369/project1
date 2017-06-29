#pragma once
#ifndef PBRS_CORE_INTERACTION_H
#define PBRS_CORE_INTERACTION_H

#include "pbrs.h"
#include "geometry.h"
#include "transform.h"
namespace pbrs {
	struct  Interaction {
		Interaction() :time(0) {}
		Interaction(const Point3f &p, const Normal3f &n, const Vector3f &pError,const Vector3f &wo, Float time) :p(p), time(time), pError(pError), wo(wo), n(n) {}
		bool IsSurfaceInteraction() const { return n != Normal3f(); }
		
		Normal3f n;//For interactions on surfaces, n stores the surface normal at the point.
		Vector3f wo;// negative ray direction
		Vector3f pError;
		Point3f p;
		Float time;
		//MediumInterface mediumInterface;
	};
	class SurfaceInteraction :public Interaction {
	public:
		SurfaceInteraction() {}
		SurfaceInteraction(const Point3f &p,
			const Vector3f &pError, const Point2f &uv, const Vector3f &wo,
			const Vector3f &dpdu, const Vector3f &dpdv,
			const Normal3f &dndu, const Normal3f &dndv,
			Float time, const Shape *shape);
		void SetShadingGeometry(const Vector3f &dpdus,
			const Vector3f &dpdvs, const Normal3f &dndus,
			const Normal3f &dndvs, bool orientationIsAuthoritative);
		Point2f uv;
		Vector3f dpdu, dpdv;
		Normal3f dndu, dndv;
		struct {
			Normal3f n;
			Vector3f dpdu, dpdv;
			Normal3f dndu, dndv;
		}shading;
		const Shape *shape = nullptr;
	};//class SurfaceInteraction 

	

	


}//namespace pbrs






















#endif // !PBRS_CORE_INTERACTION_H
