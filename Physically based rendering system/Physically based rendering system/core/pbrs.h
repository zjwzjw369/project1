#pragma once
#ifndef PBRS_CORE_PBRT_H
#define PBRS_CORE_PBRT_H

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <assert.h>

namespace pbrs{
	
	inline float Clamp(float x, float low = 0, float high = 1) {
		return (x < high) ? ((x > low) ? x : low) : high;
	}
	const float INF = std::numeric_limits<float>::infinity();

}// namespace pbrs






#endif // !PBRS_CORE_PBRT_H
