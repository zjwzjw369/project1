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
	static constexpr Float MaxFloat = std::numeric_limits<Float>::max();
	static constexpr Float Infinity = std::numeric_limits<Float>::infinity();
	static Float Pi = 3.14159265358979323846;
	inline Float Clamp(Float x, Float low = 0, Float high = 1) {
		return (x < high) ? ((x > low) ? x : low) : high;
	}
	inline Float gamma(int n) {
		return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
	}
	inline uint32_t FloatToBits(float f) {
		uint32_t ui;
		memcpy(&ui, &f, sizeof(float));
		return ui;
	}
	inline float BitsToFloat(uint32_t ui) {
		float f;
		memcpy(&f, &ui, sizeof(uint32_t));
		return f;
	}
	inline uint64_t DoubleToBits(double lf) {
		uint64_t ui;
		memcpy(&ui, &lf, sizeof(double));
		return ui;
	}
	inline double BitsToDouble(uint64_t ui) {
		double lf;
		memcpy(&lf, &ui, sizeof(uint64_t));
		return lf;
	}
	inline float NextFloatUp(float v) {
		if (std::isinf(v) && v > 0.)
			return v;
		if (v == -0.f)
			v = 0.f;
		uint32_t ui = FloatToBits(v);
		if (v >= 0)++ui;
		else --ui;
		return BitsToFloat(ui);
	}
	inline float NextFloatDown(float v) {
		if (std::isinf(v) && v < 0.) return v;
		if (v == 0.f) v = -0.f;
		uint32_t ui = FloatToBits(v);
		if (v > 0)
			--ui;
		else
			++ui;
		return BitsToFloat(ui);
	}
	inline double NextFloatUp(double v, int delta = 1) {
		if (std::isinf(v) && v > 0.) return v;
		if (v == -0.f) v = 0.f;
		uint64_t ui = FloatToBits(v);
		if (v >= 0.)
			ui += delta;
		else
			ui -= delta;
		return BitsToFloat(ui);
	}

	inline double NextFloatDown(double v, int delta = 1) {
		if (std::isinf(v) && v < 0.) return v;
		if (v == 0.f) v = -0.f;
		uint64_t ui = FloatToBits(v);
		if (v > 0.)
			ui -= delta;
		else
			ui += delta;
		return BitsToFloat(ui);
	}

	class SurfaceInteraction;
	const Float INF = std::numeric_limits<Float>::infinity();
	inline Float Radians(Float deg) { return (Pi / 180)*deg; }
}// namespace pbrs






#endif // !PBRS_CORE_PBRT_H
