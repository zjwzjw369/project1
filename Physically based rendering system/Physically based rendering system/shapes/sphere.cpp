#include "sphere.h"

namespace pbrs {
	bool Sphere::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture) const {
		Float phi;
		Point3f pHit;
		//将Ray变换到物体坐标
		Vector3f oErr, dErr;//光线的起始点和方向的误差
		Ray ray = (*WorldToObject)(r, &oErr, &dErr);
		//计算二次球系数

		//解二次方程式求t

		//计算球与光线的焦点和φ

		//测试球交点是否正确

		//找球的二次表达式的交点

		//计算球的交点的误差范围

		//初始化SurfaceInteraction的参数信息

		//更新tHit

		return true;
	}
}