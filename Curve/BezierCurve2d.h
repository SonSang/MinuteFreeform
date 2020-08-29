/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_BEZIER_CURVE_2D_H__
#define __MN_BEZIER_CURVE_2D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"

namespace MN {
	class BezierCurve2d : public Freeform2dc {
	private:
		BezierCurve2d() = default;

		// Control points for calculating derivatives up to 3rd order
		// @WARNING : Degree multiplication is already done in those matrices e.g) When 3rd degree Bezier is differentiated once, 6 must be multiplied. 
		ControlPoints derivMatT;
		ControlPoints derivMatTT;
		ControlPoints derivMatTTT;

		static void subdivideCpts(const ControlPoints& cpts, Real t, ControlPoints& lower, ControlPoints& upper);
	public:
		using Ptr = std::shared_ptr<BezierCurve2d>;

		// @buildMat : Option for building derivMats in creation time
		static BezierCurve2d create();
		static BezierCurve2d create(int degree, const ControlPoints& cpts, bool buildMat = true);
		static Ptr createPtr(int degree, const ControlPoints& cpts, bool buildMat = true);

		void updateDerivMat();	// Update deriv matrices with current control points

		virtual Vec2 evaluate(Real t) const;
		virtual Vec2 differentiate(Real t, int order) const;
		void subdivide(Real t, BezierCurve2d& lower, BezierCurve2d& upper) const;
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