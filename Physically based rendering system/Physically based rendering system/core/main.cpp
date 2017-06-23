#include "ports.h"

using namespace pbrs;

int main() {
	Ray r(Vector3f(0.0, 0.0, 0.0), Vector3f(1.0, 0.0, 0.0));
	Bounds3f b(Vector3f(1.0, 1.0, 3.0), Vector3f(-1.0, -2.0, -3.0));
	Vector3f center;
	Float rid;
	b.BoundingSphere(&center, &rid);
	Float f[4][4];
	memset(f, 0, 16 * sizeof(Float));
	f[0][0] = 1;
	f[1][1] = 1;
	f[2][2] = 1;
	f[3][3] = 1;
	
	Matrix4x4 m;
	Transform t(f);
	std::cout << t.HasScale() << std::endl;
	std::cout <<m<<" "<<f[1][2];
	system("pause");
	return 0;
}