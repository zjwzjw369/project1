#include "ports.h"

using namespace pbrs;

int main() {
	//Interaction i;
	//std::cout << i.p<<std::endl;
	Ray r(Point3f(0.0, 0.0, 0.0), Vector3f(1.0, 0.0, 0.0));
	Bounds3f b(Point3f(1.0, 1.0, 3.0), Point3f(-1.0, -2.0, -3.0));
	Point3f center;
	Float rid;
	b.BoundingSphere(&center, &rid);
	Float f[4][4];
	memset(f, 0, 16 * sizeof(Float));
	f[0][0] = 1;
	f[1][1] = 1;
	f[2][2] = 1;
	f[3][3] = 1;
	Point3<Float> p;
	Point3f p222;
	Normal3f n(1.0, 1.0, 1.0);
	Vector3f v(2.0, 2.0, 2.0);
	Faceforward(n, v);
	Matrix4x4 m;
	Transform t(f);
	
	std::cout << t.HasScale() << std::endl;
	std::cout <<m<<" "<<t(Vector3f(1.0,2.0,3.0));
	system("pause");
	return 0;
}