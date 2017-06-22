#include "ports.h"

using namespace pbrs;

int main() {
	Ray r(Vector3f(0.0, 0.0, 0.0), Vector3f(1.0, 0.0, 0.0));
	Bounds3f b(Vector3f(1.0, 1.0, 3.0), Vector3f(-1.0, -2.0, -3.0));
	Vector3f center;
	float rid;
	b.BoundingSphere(&center, &rid);
	std::cout <<center<<""<<rid;

	system("pause");
	return 0;
}