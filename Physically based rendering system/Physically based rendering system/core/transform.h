#pragma once
#ifndef PBRS_CORE_TRANSFORM_H
#define PBRS_CORE_TRANSFORM_H

#include "pbrs.h"
#include "geometry.h"

namespace pbrs {
	struct Matrix4x4{
		// Matrix4x4 Public Methods
		Matrix4x4() {
			m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
			m[0][1] = m[0][2] = m[0][3] = m[1][0] = m[1][2] = m[1][3] = m[2][0] =
				m[2][1] = m[2][3] = m[3][0] = m[3][1] = m[3][2] = 0.f;
		}
		Matrix4x4(Float mat[4][4]);
		Matrix4x4(Float t00, Float t01, Float t02, Float t03, Float t10, Float t11,
			Float t12, Float t13, Float t20, Float t21, Float t22, Float t23,
			Float t30, Float t31, Float t32, Float t33);
		bool operator==(const Matrix4x4 &m2) const {
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					if (m[i][j] != m2.m[i][j]) return false;
			return true;
		}
		bool operator!=(const Matrix4x4 &m2) const {
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					if (m[i][j] != m2.m[i][j]) return true;
			return false;
		}
		friend Matrix4x4 Transpose(const Matrix4x4 &);
		void Print(FILE *f) const {
			fprintf(f, "[ ");
			for (int i = 0; i < 4; ++i) {
				fprintf(f, "  [ ");
				for (int j = 0; j < 4; ++j) {
					fprintf(f, "%f", m[i][j]);
					if (j != 3) fprintf(f, ", ");
				}
				fprintf(f, " ]\n");
			}
			fprintf(f, " ] ");
		}
		static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
			Matrix4x4 r;
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
					m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
			return r;
		}
		friend Matrix4x4 Inverse(const Matrix4x4 &);

		friend std::ostream &operator<<(std::ostream &os, const Matrix4x4 &m) {
			os << "[  " << m.m[0][0] << ", " << m.m[0][1] << ", " << m.m[0][2] << ", " << m.m[0][3] << "  " << std::endl << "   "
				<< m.m[1][0] << ", " << m.m[1][1] << ", " << m.m[1][2] << ", " << m.m[1][3] << std::endl <<"  "
				" " << m.m[2][0] << ", " << m.m[2][1] << ", " << m.m[2][2] << ", " << m.m[2][3] << " " << std::endl << "   "
				<<m.m[3][0]<<", "<<m.m[3][1]<<", "<<m.m[3][2]<<", "<<m.m[3][3]<<"]";
			return os;
		}

		Float m[4][4];
	};//class Matrix4x4
	
	class Transform
	{
	public:
		Transform() {};
		Transform(const Float mat[4][4]) {
			m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
				mat[1][0], mat[1][1], mat[1][2], mat[1][3],
				mat[2][0], mat[2][1], mat[2][2], mat[2][3],
				mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
			mInv = Inverse(m);
		}
		Transform(const Matrix4x4 &m) :m(m), mInv(Inverse(m)) {}
		Transform(const Matrix4x4 &m,const Matrix4x4 &mInv) :m(m), mInv(mInv) {}
		friend Transform Inverse(const Transform &t) {
			return Transform(t.mInv, t.m);
		}
		friend Transform Transpose(const Transform &t) {
			return Transform(Transpose(t.m), Transpose(t.mInv));
		}
		bool operator==(const Transform &t) const {
			return t.m == m&&t.mInv == mInv;
		}
		bool operator!=(const Transform &t) const {
			return t.m != m || t.mInv != mInv;
		}
		bool operator<(const Transform &t2) const {
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j) {
					if (m.m[i][j] < t2.m.m[i][j]) return true;
					if (m.m[i][j] > t2.m.m[i][j]) return false;
				}
			return false;
		}

		template <typename T>
		inline Vector3<T> operator()(const Point<T> &p) const {
			T x = p.x, y = p.y, z = p.z;
			T xp = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3];
			T yp = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3];
			T zp = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3];
			T wp = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];
			if (wp == 1)//��pΪ����ʱ wp����0
				return Vector3<T>(xp, yp, zp);
			else
				return Vector3<T>(xp, yp, zp) / wp;
		}
		template <typename T>
		inline Vector3<T> operator()(const Vector3<T> &v) const {
			T x = v.x, y = v.y, z = v.z;
			return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
				m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
				m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
		}

		bool HasScale() const {
			Float la2 = (*this)(Vector3f(1, 0, 0)).Norm2_squared();
			Float lb2 = (*this)(Vector3f(0, 1, 0)).Norm2_squared();
			Float lc2 = (*this)(Vector3f(0, 0, 1)).Norm2_squared();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
			return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
		}
		bool IsIdentity() const {
			return (m.m[0][0] == 1.f && m.m[0][1] == 0.f && m.m[0][2] == 0.f &&
				m.m[0][3] == 0.f && m.m[1][0] == 0.f && m.m[1][1] == 1.f &&
				m.m[1][2] == 0.f && m.m[1][3] == 0.f && m.m[2][0] == 0.f &&
				m.m[2][1] == 0.f && m.m[2][2] == 1.f && m.m[2][3] == 0.f &&
				m.m[3][0] == 0.f && m.m[3][1] == 0.f && m.m[3][2] == 0.f &&
				m.m[3][3] == 1.f);
		}
		const Matrix4x4 &GetMatrix() const { return m; }
		const Matrix4x4 &GetInverseMatrix() const { return mInv; }
	private:
		Matrix4x4 m, mInv;
	};
	Transform translate(const Vector3f &delta);
	Transform Scale(Float x, Float y, Float z);
	Transform RotateX(Float theta);
	Transform RotateY(Float theta);
	Transform RotateZ(Float theta);
	Transform LookAt(const Vector3f &pos, const Vector3f &look, const Vector3f &up);
}









#endif // !PBRS_CORE_TRANSFORM_H
