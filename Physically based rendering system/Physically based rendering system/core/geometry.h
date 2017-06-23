#pragma once
#ifndef PBRS_CORE_GEOMETRY
#define PBRS_CORE_GEOMETRY

#include "pbrs.h"
namespace pbrs {

	template<typename T> class Vector3 {

	public:
		Vector3(T a = 0) : x(a), y(a), z(a) {}
		Vector3(T x, T y, T z) : x(x), y(y), z(z) {}
		Vector3(const Vector3<T> &v) : x(v.x), y(v.y), z(v.z) {}

		inline bool HasNaNs() const {
			return std::isnan(x) || isnan(y) || isnan(z);
		}

		inline Vector3<T> &operator=(const Vector3<T> &v) {
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}

		inline Vector3<T> operator-() const {
			return Vector3<T>(-x, -y, -z);
		}
		inline Vector3<T> operator+(const Vector3<T> &v) const {
			return Vector3<T>(x + v.x, y + v.y, z + v.z);
		}
		inline Vector3<T> operator-(const Vector3<T> &v) const {
			return Vector3<T>(x - v.x, y - v.y, z - v.z);
		}
		inline Vector3<T> operator*(const Vector3<T> &v) const {
			return Vector3<T>(x * v.x, y * v.y, z * v.z);
		}
		inline Vector3<T> operator/(const Vector3<T> &v) const {
			return Vector3<T>(x / v.x, y / v.y, z / v.z);
		}
		inline Vector3<T> operator+(T a) const {
			return Vector3<T>(x + a, y + a, z + a);
		}
		inline Vector3<T> operator-(T a) const {
			return Vector3<T>(x - a, y - a, z - a);
		}
		inline Vector3<T> operator*(T a) const {
			return Vector3<T>(x * a, y * a, z * a);
		}
		inline Vector3<T> operator/(T a) const {
			const T inv_a = 1.0f / a;
			return Vector3<T>(x * inv_a, y * inv_a, z * inv_a);
		}
		friend inline Vector3<T> operator+(T a, const Vector3<T> &v) {
			return Vector3<T>(v.x + a, v.y + a, v.z + a);
		}
		friend inline Vector3<T> operator-(T a, const Vector3<T> &v) {
			return Vector3<T>(v.x - a, v.y - a, v.z - a);
		}
		friend inline Vector3<T> operator*(T a, const Vector3<T> &v) {
			return Vector3<T>(v.x * a, v.y * a, v.z * a);
		}
		friend inline Vector3<T> operator/(T a, const Vector3<T> &v) {
			const T inv_a = 1.0f / a;
			return Vector3<T>(v.x * inv_a, v.y * inv_a, v.z * inv_a);
		}

		inline Vector3<T> &operator+=(const Vector3<T> &v) {
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}
		inline Vector3<T> &operator-=(const Vector3<T> &v) {
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}
		inline Vector3<T> &operator*=(const Vector3<T> &v) {
			x *= v.x;
			y *= v.y;
			z *= v.z;
			return *this;
		}
		inline Vector3<T> &operator/=(const Vector3<T> &v) {
			x /= v.x;
			y /= v.y;
			z /= v.z;
			return *this;
		}
		inline Vector3<T> &operator+=(T a) {
			x += a;
			y += a;
			z += a;
			return *this;
		}
		inline Vector3<T> &operator-=(T a) {
			x -= a;
			y -= a;
			z -= a;
			return *this;
		}
		inline Vector3<T> &operator*=(T a) {
			x *= a;
			y *= a;
			z *= a;
			return *this;
		}
		inline Vector3<T> &operator/=(T a) {
			const T inv_a = 1.0f / a;
			x *= inv_a;
			y *= inv_a;
			z *= inv_a;
			return *this;
		}

		inline T Dot(const Vector3<T> &v) const {
			return x * v.x + y * v.y + z * v.z;
		}
		inline Vector3<T> Cross(const Vector3<T> &v) const {
			return Vector3<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
		}

		inline bool operator==(const Vector3<T> &v) const {
			return x == v.x && y == v.y && z == v.z;
		}
		inline bool operator!=(const Vector3<T> &v) const {
			return x != v.x || y != v.y || z != v.z;
		}
		inline bool operator<(const Vector3<T> &v) const {
			return x < v.x && y < v.y && z < v.z;
		}
		inline bool operator<=(const Vector3<T> &v) const {
			return x <= v.x && y <= v.y && z <= v.z;
		}
		inline bool operator>(const Vector3<T> &v) const {
			return x > v.x && y > v.y && z > v.z;
		}
		inline bool operator>=(const Vector3<T> &v) const {
			return x >= v.x && y >= v.y && z >= v.z;
		}

		friend inline Vector3<T> Sqrt(const Vector3<T> &v) {
			return Vector3<T>(sqrt(v.x), sqrt(v.y), sqrt(v.z));
		}
		friend inline Vector3<T> Pow(const Vector3<T> &v, T a) {
			return Vector3<T>(pow(v.x, a), pow(v.y, a), pow(v.z, a));
		}
		friend inline Vector3<T> Abs(const Vector3<T> &v) {
			return Vector3<T>(abs(v.x), abs(v.y), abs(v.z));
		}
		friend inline Vector3<T> Min(const Vector3<T> &v1, const Vector3<T> &v2) {
			return Vector3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
		}
		friend inline Vector3<T> Max(const Vector3<T> &v1, const Vector3<T> &v2) {
			return Vector3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
		}
		friend inline Vector3<T> Round(const Vector3<T> &v) {
			return Vector3<T>(std::round(v.x), std::round(v.y), std::round(v.z));
		}
		friend inline Vector3<T> Floor(const Vector3<T> &v) {
			return Vector3<T>(std::floor(v.x), std::floor(v.y), std::floor(v.z));
		}
		friend inline Vector3<T> Ceil(const Vector3<T> &v) {
			return Vector3<T>(std::ceil(v.x), std::ceil(v.y), std::ceil(v.z));
		}
		friend inline Vector3<T> Trunc(const Vector3<T> &v) {
			return Vector3<T>(std::trunc(v.x), std::trunc(v.y), std::trunc(v.z));
		}
		friend inline Vector3<T> Clamp(const Vector3<T> &v, T low = 0, T high = 1) {
			return Vector3<T>(pbrs::Clamp(v.x, low, high), pbrs::Clamp(v.y, low, high), pbrs::Clamp(v.z, low, high));
		}
		friend inline Vector3<T> Lerp(T a, const Vector3<T> &v1, const Vector3<T> &v2) {
			return v1 + a * (v2 - v1);
		}
		friend inline Vector3<T> Permute(const Vector3<T> &v, int x, int y, int z) {
			return Vector3<T>(v[x], v[y], v[z]);
		}

		inline T operator[](size_t i) const {
			return raw[i];
		}
		inline T &operator[](size_t i) {
			return raw[i];
		}

		inline int MinDimension() const {
			return (x < y && x < z) ? 0 : ((y < z) ? 1 : 2);
		}
		inline int MaxDimension() const {
			return (x > y && x > z) ? 0 : ((y > z) ? 1 : 2);
		}
		inline T Min() const {
			return (x < y && x < z) ? x : ((y < z) ? y : z);
		}
		inline T Max() const {
			return (x > y && x > z) ? x : ((y > z) ? y : z);
		}

		inline T Norm2_squared() const {
			return x * x + y * y + z * z;
		}
		inline T Norm2() const {
			return sqrt(Norm2_squared());
		}
		inline Vector3<T> &Normalize() {
			const T a = 1 / Norm2();
			x *= a;
			y *= a;
			z *= a;
			return *this;
		}

		friend inline std::ostream &operator<<(std::ostream& os, const Vector3<T> &v) {
			os << '[' << v.x << ' ' << v.y << ' ' << v.z << ']';
			return os;
		}

		union {
			struct {
				T x, y, z;
			};
			T raw[3];
		};
	};//class Vector
	typedef Vector3<Float> Vector3f;
	typedef Vector3<int> Vector3i;


	class Ray{
	public:
		Ray() :tMax(INF), time(0.0f) {}
		Ray(const Vector3f &o, const Vector3f &d, Float tMax = INF, Float time = 0.0f) :o(0), d(d), tMax(tMax), time(time) {}
		inline Vector3f getPoint(const Float t) const { return o + t*d; }
		Vector3f operator()(Float t) const { return o + t*d; }
		Vector3f o,d;
		mutable Float tMax;
		Float time;
		//const Medium *medium;
	};//class ray

	template<typename T> class Bounds3 {//AABB:axis-aligned bounding boxes
	public:
		Bounds3() {
			T minNum = std::numeric_limits<T>::lowest();
			T maxNum = std::numeric_limits<T>::max();
			pMin = Vector3<T>(maxNum, maxNum, maxNum);
			pMax = Vector3<T>(minNum, minNum, minNum);
		}
		explicit Bounds3(const Vector3<T> &p) :pMin(p), pMax(p) {}
		Bounds3(const Vector3<T> &p1, const Vector3<T> &p2):
			pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)),
			pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)){}
		const Vector3<T> operator[](int i) const { if (i == 0)return pMin; else return pMax; }
		Vector3<T> &operator[](int i) { if (i == 0)return pMin; else return pMax; }
		Vector3<T> Corner(int corner) const {
			return Vector3<T>((*this)[(corner & 1)].x,
				(*this)[(corner & 2) ? 1 : 0].y,
				(*this)[(corner & 4) ? 1 : 0].z);
		}
		Vector3<T> Diagonal() const { return pMax - pMin; }
		T SurfaceArea() const {
			Vector3<T> d = Diagonal();
			return 2 * (d.x*d.y + d.x*d.z + d.y*d.z);
		}
		T volume() const {
			Vector3<T> d = Diagonal();
			return d.x*d.y*d.z;
		}
		int MaximumExtent() const {
			Vector3<T> d = Diagonal();
			if (d.x > d.z&&d.x > d.z)
				return 0;
			else if (d.y > d.z)
				return 1;
			else
				return 2;
		}
		Vector3<T> Lerp(const Vector3f &t)const {
			return Vector3<T>(::Lerp(t.x, pMin.x, pMax.x), ::Lerp(t.y, pMin.y, pMax.y), ::Lerp(t.z, pMin.z, pMax.z));
		}
		Vector3<T> Offset(const Vector3<T> &p) const{
			Vector3<T> o = p - pMin;
			if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
			if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
			if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
			return o;
		}
		void BoundingSphere(Vector3<T> *center,Float *radius) const{
			*center = (pMin + pMax) / 2;
			*radius = Inside(*center, *this) ? (pMax - *center).Norm2() : 0;
		}
		Vector3<T> pMin, pMax;
	};//class Bounds

	template<typename T> 
	Bounds3<T> Union(const Bounds3<T> &b, const Vector3<T> & p) {
		return Bounds3<T>(Point3<T>(std::min(b.pMin.x, p.x), std::min(b.pMin.y, p.y), std::min(b.pMin.z, p.z)),
			Point3<T>(std::max(b.pMax.x, p.x), std::max(b.pMax.y, p.y), std::max(b.pMax.z, p.z)));
	}
	template<typename T>
	Bounds3<T> Union(const Bounds3<T> &b1, const Bounds3<T> &b2) {
		return Bounds3<T>(Point3<T>(std::min(b1.pMin.x, b2.pMin.x), std::min(b1.pMin.y, b2.pMin.y), std::min(b1.pMin.z, b2.pMin.z)),
			Point3<T>(std::max(b1.pMax.x, b2.pMax.x), std::max(b1.pMax.y, b2.pMax.y), std::max(b1.pMax.z, b2.pMax.z)));
	}
	template<typename T> //返回重叠部分
	Bounds3<T> Intersect(const Bounds3<T> &b1, const Bounds3<T> &b2) {
		return Bounds3<T>(Point3<T>(std::max(b1.pMin.x, b2.pMin.x),std::max(b1.pMin.y, b2.pMin.y),std::max(b1.pMin.z, b2.pMin.z)),
			Point3<T>(std::min(b1.pMax.x, b2.pMax.x),std::min(b1.pMax.y, b2.pMax.y),std::min(b1.pMax.z, b2.pMax.z)));
	}
	template<typename T> //返回是否重叠
	bool Overlaps(const Bounds3<T> &b1,Bounds3<T> &b2) {
		bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
		bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
		bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
		return (x && y && z);
	}
	template<typename T>
	bool Inside(const Vector3<T> &p,const Bounds3<T> &b) {
		return (p.x >= b.pMin.x && p.x <= b.pMax.x &&p.y >= b.pMin.y && p.y <= b.pMax.y &&p.z >= b.pMin.z && p.z <= b.pMax.z);
	}
	template<typename T>//不包含上界，用于整形包围
	bool InsideExclusive(const Vector3<T> &p, const Bounds3<T> &b) {
		return (p.x >= b.pMin.x && p.x < b.pMax.x &&p.y >= b.pMin.y && p.y < b.pMax.y &&p.z >= b.pMin.z && p.z < b.pMax.z);
	}
	template<typename T, typename U >
	inline Bounds3<T> Expand(const Bounds3<T> &b, U delta) {
		return Bounds3<T>(b.pMin - Vector3<T>(delta, delta, delta), b.pMax + Vector3<T>(delta, delta, delta))
	}

	typedef Bounds3<Float> Bounds3f;
	typedef Bounds3<int> Bounds3i;

}// namespace pbrs

#endif // !PBRS_CORE_GEOMETRY