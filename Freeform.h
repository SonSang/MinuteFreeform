/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_FREEFORM_H__
#define __MN_FREEFORM_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "MinuteUtils/utils.h"
#include <vector>
#include <memory>

namespace MN {
	using BasisVector = std::vector<Real>;
	using KnotVector = std::vector<Real>;
	const Binomial Bin16 = Binomial::create(16);

	/* 
	 * Base class of all freeform entities. e.g. Bezier, B-Spline
	 * @Dimension	: Dimension of this freeform object. It equals to #variables it depends on
	 *				   e.g) Bezier curve has one variable [ t ], so its dimension is one
	 * @T			: Type of point. It could be a vector with 3 values, or planar points with 2 values
	 */
	template<int dimension, typename T>
	class Freeform {	
	};

	// ============================================================= Freeform 2d
	template<int dimension>
	class Freeform2d : public Freeform<dimension, Vec2> {
	public:
	};

	// ============================================================= Freeform 2d curve
	class Freeform2dc : public Freeform2d<1> {
	public:
		using ControlPoints = std::vector<Vec2>;
		using Ptr = std::shared_ptr<Freeform2dc>;
	protected:
		ControlPoints cpts;
		Domain	domain = Domain::create(0, 1);
		int		degree;

		Freeform2dc() = default;
		inline static Vec2 tensorProduct(const BasisVector& basis, const ControlPoints& tensor) {
			Vec2 vec{ 0, 0 };

			int num = (int)basis.size();
			for (int r = 0; r < num; r++)
				vec += tensor[r] * basis[r];
			return vec;
		}
	public:
		inline ControlPoints& getCpts() noexcept {
			return cpts;
		}
		inline const ControlPoints& getCptsC() const noexcept {
			return cpts;
		}
		inline void setCpts(const ControlPoints& cpts) noexcept {
			this->cpts = cpts;
		}

		inline void setDomain(const Domain& domain) noexcept {
			this->domain = domain;
		}
		inline Domain& getDomain() noexcept {
			return domain;
		}
		inline const Domain& getDomainC() const noexcept {
			return domain;
		}

		inline void setDegree(int degree) noexcept {
			this->degree = degree;
		}
		inline int getDegree() const noexcept {
			return degree;
		}

		inline bool validateParam(Real t) const noexcept {
			return domain.has(t);
		}
		inline bool validateDomain(const Domain& domain) const noexcept {
			return this->domain.has(domain);
		}

		static inline Vec2 normal(const Vec2& firstD, const Vec2& secondD) {
			Real factor = Vec2::factorize(secondD, firstD, false);
			return secondD - firstD * factor;
		}
		static inline Real curvature(const Vec2& firstD, const Vec2& secondD) {
			Real x1, x2, y1, y2;
			x1 = firstD[0];
			x2 = secondD[0];
			y1 = firstD[1];
			y2 = secondD[1];

			Real denom = fabs(x1 * y2 - x2 * y1);
			Real nom = sqrt(x1 * x1 + y1 * y1);
			nom = nom * nom * nom;
			return denom / nom;
		}

		inline virtual Vec2 evaluate(Real t) const {
			return Vec2();
		}
		inline virtual Vec2 differentiate(Real t, int order) const {
			return Vec2();
		}
		inline virtual Vec2 normal(Real t) const {
			Vec2 firstD, secondD;
			firstD = differentiate(t, 1);
			secondD = differentiate(t, 2);
			return normal(firstD, secondD);
		}
		inline virtual Real curvature(Real t) const {
			Vec2 firstD, secondD;
			firstD = differentiate(t, 1);
			secondD = differentiate(t, 2);

			return curvature(firstD, secondD);
		}
	};

	// =============================================================  Freeform 2d surface
	class Freeform2ds : public Freeform2d<2> {
	public:
		using ControlPoints = std::vector<std::vector<Vec2>>;
		using Ptr = std::shared_ptr<Freeform2ds>;
	protected:
		ControlPoints cpts;
		Domain	uDomain = Domain::create(0, 1);
		Domain	vDomain = Domain::create(0, 1);
		int		uDegree;
		int		vDegree;

		Freeform2ds() = default;
		inline static Vec2 tensorProduct(const BasisVector& left, const ControlPoints& tensor, const BasisVector& right) {
			Vec2 vec{ 0, 0 };

			int row = (int)left.size();
			int col = (int)right.size();

			for (int r = 0; r < row; r++)
				for (int c = 0; c < col; c++)
					vec += tensor[r][c] * left[r] * right[c];
			return vec;
		}
	public:
		inline ControlPoints& getCpts() noexcept {
			return cpts;
		}
		inline const ControlPoints& getCptsC() const noexcept {
			return cpts;
		}
		inline void setCpts(const ControlPoints& cpts) noexcept {
			this->cpts = cpts;
		}

		// Direction : 0 for U, 1 for V
		inline void setDomain(int dir, Real beg, Real end) noexcept {
			if (dir == 0)
				uDomain.set(beg, end);
			else
				vDomain.set(beg, end);
		}
		inline void setDomain(int dir, const Domain& domain) noexcept {
			if (dir == 0)
				uDomain = domain;
			else
				vDomain = domain;
		}
		inline const Domain& getDomainC(int dir) const noexcept {
			if (dir == 0)
				return uDomain;
			else
				return vDomain;
		}

		inline void setDegree(int dir, int degree) noexcept {
			if (dir == 0)
				uDegree = degree;
			else
				vDegree = degree;
		}
		inline int getDegree(int dir) const noexcept {
			if (dir == 0)
				return uDegree;
			else
				return vDegree;
		}

		inline bool validateParam(Real u, Real v) const noexcept {
			return uDomain.has(u) && vDomain.has(v);
		}
		inline bool validateDomain(const Domain& uDomain, const Domain& vDomain) const noexcept {
			return this->uDomain.has(uDomain) && this->vDomain.has(vDomain);
		}

		inline virtual Vec2 evaluate(double u, double v) const {
			return Vec2();
		}
		inline virtual Vec2 differentiate(double u, double v, int u_order, int v_order) const {
			return Vec2();
		}
	};

	// ============================================================= Freeform 3d
	template<int dimension>
	class Freeform3d : public Freeform<dimension, Vec3> {
	public:
	};

	// ============================================================= Freeform 3d curve
	class Freeform3dc : public Freeform3d<1> {
	public:
		using ControlPoints = std::vector<Vec3>;
		using Ptr = std::shared_ptr<Freeform3dc>;
	protected:
		ControlPoints cpts;
		Domain	domain = Domain::create(0, 1);
		int		degree;

		Freeform3dc() = default;
		inline static Vec3 tensorProduct(const BasisVector& basis, const ControlPoints& tensor) {
			Vec3 vec{ 0, 0, 0 };

			int num = (int)basis.size();
			for (int r = 0; r < num; r++)
				vec += tensor[r] * basis[r];
			return vec;
		}
	public:
		inline ControlPoints& getCpts() noexcept {
			return cpts;
		}
		inline const ControlPoints& getCptsC() const noexcept {
			return cpts;
		}
		inline void setCpts(const ControlPoints& cpts) noexcept {
			this->cpts = cpts;
		}

		inline void setDomain(const Domain& domain) noexcept {
			this->domain = domain;
		}
		inline Domain& getDomain() noexcept {
			return domain;
		}
		inline const Domain& getDomainC() const noexcept {
			return domain;
		}

		inline void setDegree(int degree) noexcept {
			this->degree = degree;
		}
		inline int getDegree() const noexcept {
			return degree;
		}

		inline bool validateParam(Real t) const noexcept {
			return domain.has(t);
		}
		inline bool validateDomain(const Domain& domain) const noexcept {
			return this->domain.has(domain);
		}

		inline static Vec3 normal(const Vec3& firstD, const Vec3& secondD) {
			Real factor = Vec3::factorize(secondD, firstD);
			return secondD - firstD * factor;
		}
		inline static Real curvature(const Vec3& firstD, const Vec3& secondD) {
			Vec3 cross = firstD.cross(secondD);
			Real denom = cross.len();
			Real nom = firstD.len();
			return denom / (nom * nom * nom);
		}

		inline virtual Vec3 evaluate(Real t) const {
			return Vec3();
		}
		inline virtual Vec3 differentiate(Real t, int order) const {
			return Vec3();
		}
		inline virtual Vec3 normal(Real t) const {
			Vec3 firstD = differentiate(t, 1);
			Vec3 secondD = differentiate(t, 2);
			return normal(firstD, secondD);
		}
		inline virtual Real curvature(Real t) const {
			Vec3 firstD, secondD;
			firstD = differentiate(t, 1);
			secondD = differentiate(t, 2);

			return curvature(firstD, secondD);
		}
	};

	// ============================================================= Freeform 3d surface
	class Freeform3ds : public Freeform3d<2> {
	public:
		using ControlPoints = std::vector<std::vector<Vec3>>;
		using Ptr = std::shared_ptr<Freeform3ds>;
	protected:
		ControlPoints cpts;
		Domain	uDomain = Domain::create(0, 1);
		Domain	vDomain = Domain::create(0, 1);
		int		uDegree;
		int		vDegree;

		Freeform3ds() = default;
		inline static Vec3 tensorProduct(const BasisVector& left, const ControlPoints& tensor, const BasisVector& right) {
			Vec3 vec{ 0, 0, 0 };

			int row = (int)left.size();
			int col = (int)right.size();

			for (int r = 0; r < row; r++)
				for (int c = 0; c < col; c++)
					vec += tensor[r][c] * left[r] * right[c];
			return vec;
		}
	public:
		inline ControlPoints& getCpts() noexcept {
			return cpts;
		}
		inline const ControlPoints& getCptsC() const noexcept {
			return cpts;
		}
		inline void setCpts(const ControlPoints& cpts) noexcept {
			this->cpts = cpts;
		}

		// Direction : 0 for U, 1 for V
		inline void setDomain(int dir, double beg, double end) noexcept {
			if (dir == 0)
				uDomain.set(beg, end);
			else
				vDomain.set(beg, end);
		}
		inline void setDomain(int dir, const Domain& domain) noexcept {
			if (dir == 0)
				uDomain = domain;
			else
				vDomain = domain;
		}
		inline const Domain& getDomainC(int dir) const noexcept {
			if (dir == 0)
				return uDomain;
			else
				return vDomain;
		}

		inline void setDegree(int dir, int degree) noexcept {
			if (dir == 0)
				uDegree = degree;
			else
				vDegree = degree;
		}
		inline int getDegree(int dir) const noexcept {
			if (dir == 0)
				return uDegree;
			else
				return vDegree;
		}

		inline bool validateParam(double u, double v) const noexcept {
			return uDomain.has(u) && vDomain.has(v);
		}
		inline bool validateDomain(const Domain& uDomain, const Domain& vDomain) const noexcept {
			return this->uDomain.has(uDomain) && this->vDomain.has(vDomain);
		}

		inline virtual Vec3 evaluate(double u, double v) const {
			return Vec3::zero();
		}
		inline virtual Vec3 differentiate(double u, double v, int u_order, int v_order) const {
			return Vec3::zero();
		}
		inline virtual Vec3 normal(double u, double v) const {
			auto Su = differentiate(u, v, 1, 0);
			auto Sv = differentiate(u, v, 0, 1);
			auto n = Su.cross(Sv);
			n.normalize();
			return n;
		}
	};

	// ============================================================= Freeform 3d volume
	class Freeform3dv : public Freeform3d<3> {
	public:
		using ControlPoints = std::vector<std::vector<std::vector<Vec3>>>;
		using BasisVector = std::vector<Real>;
		using Ptr = std::shared_ptr<Freeform3dv>;
	protected:
		ControlPoints cpts;
		Domain	uDomain = Domain::create(0, 1);
		Domain	vDomain = Domain::create(0, 1);
		Domain  wDomain = Domain::create(0, 1);
		int		uDegree;
		int		vDegree;
		int		wDegree;

		Freeform3dv() = default;
		inline static Vec3 tensorProduct(const BasisVector& basisU, const BasisVector& basisV, const BasisVector& basisW, const ControlPoints& tensor) {
			Vec3 vec{ 0, 0, 0 };

			int uSize = (int)basisU.size();
			int vSize = (int)basisV.size();
			int wSize = (int)basisW.size();

			for (int u = 0; u < uSize; u++)
				for (int v = 0; v < vSize; v++)
					for (int w = 0; w < wSize; w++)
						vec += tensor[u][v][w] * basisU[u] * basisV[v] * basisW[w];
			return vec;
		}
	public:
		inline ControlPoints& getCpts() noexcept {
			return cpts;
		}
		inline const ControlPoints& getCptsC() const noexcept {
			return cpts;
		}
		inline void setCpts(const ControlPoints& cpts) noexcept {
			this->cpts = cpts;
		}

		// Direction : 0 for U, 1 for V, 2 for W
		inline void setDomain(int dir, Real beg, Real end) noexcept {
			if (dir == 0)
				uDomain = Domain::create(beg, end);
			else if (dir == 1)
				vDomain = Domain::create(beg, end);
			else
				wDomain = Domain::create(beg, end);
		}
		inline void setDomain(int dir, const Domain& domain) noexcept {
			if (dir == 0)
				uDomain = domain;
			else if (dir == 1)
				vDomain = domain;
			else
				wDomain = domain;
		}
		inline const Domain& getDomainC(int dir) const noexcept {
			if (dir == 0)
				return uDomain;
			else if (dir == 1)
				return vDomain;
			else
				return wDomain;
		}

		inline void setDegree(int dir, int degree) noexcept {
			if (dir == 0)
				uDegree = degree;
			else if (dir == 1)
				vDegree = degree;
			else
				wDegree = degree;
		}
		inline int getDegree(int dir) const noexcept {
			if (dir == 0)
				return uDegree;
			else if (dir == 1)
				return vDegree;
			else
				return wDegree;
		}

		inline bool validateParam(Real u, Real v, Real w) const noexcept {
			return uDomain.has(u) && vDomain.has(v) && wDomain.has(w);
		}
		inline bool validateDomain(const Domain& uDomain, const Domain& vDomain, const Domain& wDomain) const noexcept {
			return this->uDomain.has(uDomain) && this->vDomain.has(vDomain) && this->wDomain.has(wDomain);
		}

		inline virtual Vec3 evaluate(Real u, Real v, Real w) const {
			return Vec3();
		}
		inline virtual Vec3 differentiate(Real u, Real v, Real w, int uOrder, int vOrder, int wOrder) const {
			return Vec3();
		}
	};

	// Bezier
	class Bezier {
	public:
		inline static void calBasisVector(Real t, int degree, BasisVector& basis) {
			Real t_1 = 1.0 - t;
			std::vector<Real> Ts;
			std::vector<Real> T_1s;
			Ts.resize(degree + 1);
			T_1s.resize(degree + 1);
			Ts[0] = 1.0;
			T_1s[0] = 1.0;
			for (int i = 1; i < degree + 1; i++) {
				Ts[i] = Ts[i - 1] * t;
				T_1s[i] = T_1s[i - 1] * t_1;
			}

			basis.resize(degree + 1);
			for (int i = 0; i < degree + 1; i++)
				basis[i] = Bin16.at(degree, i) * T_1s[degree - i] * Ts[i];
		}
	};

	// Bspline
	class Bspline {
	public:
		inline static KnotVector createOpenUniformKnotVector(int degree, int cptsNum, const Domain& domain = Domain::create(0, 1)) {
			KnotVector knotVector;
			int
				length = cptsNum + degree + 1;
			Real
				val = domain.beg(),
				base = cptsNum - degree;
			for (int i = 0; i < degree + 1; i++)
				knotVector.push_back(val);
			for (int i = 0; i < cptsNum - degree - 1; i++) {
				val = (domain.width() * (i + 1) / base) + domain.beg();
				knotVector.push_back(val);
			}
			val = domain.end();
			for (int i = 0; i < degree + 1; i++)
				knotVector.push_back(val);
			return knotVector;
		}
	};
}

#endif