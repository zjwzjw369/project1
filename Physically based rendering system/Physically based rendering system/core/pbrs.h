#pragma once
#ifndef PBRS_CORE_PBRT_H
#define PBRS_CORE_PBRT_H

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <assert.h>
#define MachineEpsilon (std::numeric_limits<Float>::epsilon() * 0.5)
namespace pbrs{
	class Shape;
	template<typename T>
	class Vector3;
	template<typename T>
	class Point3;
	template<typename T>
	class Normal3;
	typedef double Float;
	static Float Pi = 3.14159265358979323846;
	inline Float Clamp(Float x, Float low = 0, Float high = 1) {
		return (x < high) ? ((x > low) ? x : low) : high;
	}
	inline Float gamma(int n) {
		return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
	}
	class SurfaceInteraction;
	const Float INF = std::numeric_limits<Float>::infinity();
	inline Float Radians(Float deg) { return (Pi / 180)*deg; }
}// namespace pbrs






#endif // !PBRS_CORE_PBRT_H
