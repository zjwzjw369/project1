#include "sphere.h"

namespace pbrs {
	bool Sphere::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture) const {
		Float phi;
		Point3f pHit;
		//��Ray�任����������
		Vector3f oErr, dErr;//���ߵ���ʼ��ͷ�������
		Ray ray = (*WorldToObject)(r, &oErr, &dErr);
		//���������ϵ��

		//����η���ʽ��t

		//����������ߵĽ���ͦ�

		//�����򽻵��Ƿ���ȷ

		//����Ķ��α��ʽ�Ľ���

		//������Ľ������Χ

		//��ʼ��SurfaceInteraction�Ĳ�����Ϣ

		//����tHit

		return true;
	}
}