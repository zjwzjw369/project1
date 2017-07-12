#include "sphere.h"
#include "../core/efloat.h"
namespace pbrs {
	bool Sphere::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture) const {
		Float phi;
		Point3f pHit;
		//��Ray�任����������
		Vector3f oErr, dErr;//���ߵ���ʼ��ͷ�������
		Ray ray = (*WorldToObject)(r, &oErr, &dErr);
		//����at^2 + bt + c = 0��ϵ��a,b,c
		EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
		EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
		EFloat a = dx * dx + dy * dy + dz * dz;
		EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
		EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);
		//����η���ʽ��t
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
		//����������ߵĽ���λ�úͦ�
		pHit = ray((Float)tShapeHit);
		pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
		//255
		if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
		phi = std::atan2(pHit.y, pHit.x);
		if (phi < 0) phi += 2 * Pi;
		//�����򽻵��Ƿ���ȷ

		//����Ķ��α��ʽ�Ľ���

		//������Ľ������Χ

		//��ʼ��SurfaceInteraction�Ĳ�����Ϣ

		//����tHit

		return true;
	}
}