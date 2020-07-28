#ifndef __MN_BEZIER_CURVE_3D_H__
#define __MN_BEZIER_CURVE_3D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"

namespace MN {
	class BezierCurve3d : public Freeform3dc {
	private:
		BezierCurve3d() = default;

		// Control points for calculating derivatives up to 3rd order
		// @WARNING : Degree multiplication is already done in those matrices e.g) When 3rd degree Bezier is differentiated once, 6 must be multiplied. 
		ControlPoints derivMatT;
		ControlPoints derivMatTT;
		ControlPoints derivMatTTT;

		static void subdivideCpts(const ControlPoints& cpts, Real t, ControlPoints& lower, ControlPoints& upper);
	public:
		using Ptr = std::shared_ptr<BezierCurve3d>;

		// @buildMat : Option for building derivMats in creation time
		static BezierCurve3d create();
		static BezierCurve3d create(int degree, const ControlPoints& cpts, bool buildMat = true);
		static Ptr createPtr(int degree, const ControlPoints& cpts, bool buildMat = true);

		void updateDerivMat();	// Update deriv matrices with current control points

		virtual Vec3 evaluate(Real t) const;
		virtual Vec3 differentiate(Real t, int order) const;
		void subdivide(Real t, BezierCurve3d& lower, BezierCurve3d& upper) const;
		Ptr subdivide(const Domain& subdomain) const;

		inline const ControlPoints& getDerivMatT() const noexcept {
			return derivMatT;
		}
		inline const ControlPoints& getDerivMatTT() const noexcept {
			return derivMatTT;
		}
		inline const ControlPoints& getDerivMatTTT() const noexcept {
			return derivMatTTT;
		}
	};
}

#endif