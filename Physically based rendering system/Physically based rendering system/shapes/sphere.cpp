#include "sphere.h"
#include "../core/efloat.h"
namespace pbrs {
	bool Sphere::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture) const {
		Float phi;
		Point3f pHit;
		//将Ray变换到物体坐标
		Vector3f oErr, dErr;//光线的起始点和方向的误差
		Ray ray = (*WorldToObject)(r, &oErr, &dErr);
		//计算at^2 + bt + c = 0的系数a,b,c
		EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
		EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
		EFloat a = dx * dx + dy * dy + dz * dz;
		EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
		EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);
		//解二次方程式求t
		EFloat t0, t1;
		if (!Quadratic(a, b, c, &t0, &t1))
			return false;
		if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0)
			return false;
		EFloat tShapeHit = t0;
		if (tShapeHit.LowerBound() <= 0) {
			tShapeHit = t1;
			if (tShapeHit.UpperBound() > ray.tMax)
				return false;
		}
		//计算球与光线的交点位置和φ
		pHit = ray((Float)tShapeHit);
		pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
		//255
		if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
		phi = std::atan2(pHit.y, pHit.x);
		if (phi < 0) phi += 2 * Pi;
		//测试球交点是否正确

		//找球的二次表达式的交点

		//计算球的交点的误差范围

		//初始化SurfaceInteraction的参数信息

		//更新tHit

		return true;
	}
}