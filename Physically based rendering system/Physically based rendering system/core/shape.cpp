#include "shape.h"
#include "interaction.h"
namespace pbrs {

	bool Shape::IntersectP(const Ray &ray, bool testAlphaTexture)const {
		return Intersect(ray, nullptr, nullptr, testAlphaTexture);
		return false;
	}
}