#ifndef __MN_BEZIER_SURFACE_2D_H__
#define __MN_BEZIER_SURFACE_2D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"
#include "../Curve/BezierCurve2d.h"
#include <memory>

namespace MN {
	class BezierSurface2d : public Freeform2ds {
	private:
		BezierSurface2d() = default;

		// Control points for calculating derivatives up to 2nd order
		// @WARNING : Degree multiplication is already done in those matrices e.g) When 3rd degree Bezier is differentiated, 6 must be multiplied. 
		ControlPoints derivMatU;
		ControlPoints derivMatV;
		ControlPoints derivMatUU;
		ControlPoints derivMatUV;
		ControlPoints derivMatVV;
		ControlPoints derivMatUUU;
		ControlPoints derivMatUUV;
		ControlPoints derivMatUVV;
		ControlPoints derivMatVVV;

		static void subdivideCpts(const std::vector<Vec2>& cpts, Real t, std::vector<Vec2>& lower, std::vector<Vec2>& upper);
	public:
		using Ptr = std::shared_ptr<BezierSurface2d>;
		const static Binomial binomial;

		// @buildMat : Option for building derivMats in creation time
		static BezierSurface2d create();
		static BezierSurface2d create(int uDegree, int vDegree, const ControlPoints& cpts, bool buildMat = true);
		static Ptr createPtr(int uDegree, int vDegree, const ControlPoints& cpts, bool buildMat = true);

		void updateDerivMat();	// Update deriv matrices with current control points
		virtual Vec2 evaluate(Real u, Real v) const;
		virtual Vec2 differentiate(Real u, Real v, int uOrder, int vOrder) const;

		void uSubdivide(Real u, BezierSurface2d& lower, BezierSurface2d& upper) const;
		void vSubdivide(Real v, BezierSurface2d& lower, BezierSurface2d& upper) const;
		Ptr subdivide(const Domain& uSubdomain, const Domain& vSubdomain) const;

		void uIsoCurve(Real u, BezierCurve2d& curve) const;
		void vIsoCurve(Real v, BezierCurve2d& curve) const;

		inline const ControlPoints& getDerivMatU() const noexcept {
			return derivMatU;
		}
		inline const ControlPoints& getDerivMatV() const noexcept {
			return derivMatV;
		}
		inline const ControlPoints& getDerivMatUU() const noexcept {
			return derivMatUU;
		}
		inline const ControlPoints& getDerivMatUV() const noexcept {
			return derivMatUV;
		}
		inline const ControlPoints& getDerivMatVV() const noexcept {
			return derivMatVV;
		}
		inline const ControlPoints& getDerivMatUUU() const noexcept {
			return derivMatUUU;
		}
		inline const ControlPoints& getDerivMatUUV() const noexcept {
			return derivMatUUV;
		}
		inline const ControlPoints& getDerivMatUVV() const noexcept {
			return derivMatUVV;
		}
		inline const ControlPoints& getDerivMatVVV() const noexcept {
			return derivMatVVV;
		}
	};
}

#endif