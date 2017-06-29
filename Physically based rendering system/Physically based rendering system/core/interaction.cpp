#include "interaction.h"

namespace pbrs {
	void SurfaceInteraction::SetShadingGeometry(const Vector3f &dpdus,
		const Vector3f &dpdvs, const Normal3f &dndus,
		const Normal3f &dndvs, bool orientationIsAuthoritative) {
		//Compute shading.n for SurfaceInteraction
		/*
		shading.n = Normalize((Normal3f)Cross(dpdus, dpdvs));
		if (shape && (shape->reverseOrientation ^
		shape->transformSwapsHandedness))
		shading.n = -shading.n;
		if (orientationIsAuthoritative)
		n = Faceforward(n, shading.n);
		else
		shading.n = Faceforward(shading.n, n);
		*/
		shading.dpdu = dpdus;
		shading.dpdv = dpdvs;
		shading.dndu = dndus;
		shading.dndv = dndvs;
	}
	SurfaceInteraction::SurfaceInteraction(const Point3f &p,
		const Vector3f &pError, const Point2f &uv, const Vector3f &wo,
		const Vector3f &dpdu, const Vector3f &dpdv,
		const Normal3f &dndu, const Normal3f &dndv,
		Float time, const Shape *shape) :Interaction(p, Normal3f(Normalize(Cross(dpdu, dpdv))), pError, wo,
			time), uv(uv), dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv), shape(shape) {
		shading.n = n;
		shading.dpdu = dpdu;
		shading.dpdv = dpdv;
		shading.dndu = dndu;
		shading.dndv = dndv;
		//if (shape && (shape->reverseOrientation ^
		//	shape->transformSwapsHandedness)) {
		//	n *= -1;
		//	shading.n *= -1;
		//}
	}
}