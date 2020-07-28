#ifndef __MN_BEZIER_VOLUME_3D_H__
#define __MN_BEZIER_VOLUME_3D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"
#include <memory>

namespace MN {
	class BezierVolume3d : public Freeform3dv {
	private:
		BezierVolume3d() = default;

		// Control points for calculating derivatives up to 3rd order
		// @WARNING : Degree multiplication is already done in those matrices e.g) When 3rd degree Bezier is differentiated, 6 must be multiplied. 
		ControlPoints derivMatU;
		ControlPoints derivMatV;
		ControlPoints derivMatW;

		ControlPoints derivMatUU;
		ControlPoints derivMatUV;
		ControlPoints derivMatUW;
		ControlPoints derivMatVV;
		ControlPoints derivMatVW;
		ControlPoints derivMatWW;

		ControlPoints derivMatUUU;
		ControlPoints derivMatUUV;
		ControlPoints derivMatUUW;
		ControlPoints derivMatUVV;
		ControlPoints derivMatUVW;
		ControlPoints derivMatUWW;
		ControlPoints derivMatVVV;
		ControlPoints derivMatVVW;
		ControlPoints derivMatVWW;
		ControlPoints derivMatWWW;

		static void subdivideCpts(const std::vector<Vec3>& cpts, Real t, std::vector<Vec3>& lower, std::vector<Vec3>& upper);

	public:
		using Ptr = std::shared_ptr<BezierVolume3d>;
		const static Binomial binomial;

		// @buildMat : Option for building derivMats in creation time
		static BezierVolume3d create();
		static BezierVolume3d create(int uDegree, int vDegree, int wDegree, const ControlPoints& cpts, bool buildMat = true);
		static Ptr createPtr(int uDegree, int vDegree, int wDegree, const ControlPoints& cpts, bool buildMat = true);

		void updateDerivMat();		// Update deriv matrices with current control points
		virtual Vec3 evaluate(Real u, Real v, Real w) const;
		virtual Vec3 differentiate(Real u, Real v, Real w, int uOrder, int vOrder, int wOrder) const;

		void uSubdivide(Real u, BezierVolume3d& lower, BezierVolume3d& upper, bool buildMat = true) const;
		void vSubdivide(Real v, BezierVolume3d& lower, BezierVolume3d& upper, bool buildMat = true) const;
		void wSubdivide(Real w, BezierVolume3d& lower, BezierVolume3d& upper, bool buildMat = true) const;
		Ptr subdivide(const Domain& uSubdomain, const Domain& vSubdomain, const Domain& wSubdomain) const;

		inline const ControlPoints& getDerivMatU() const noexcept {
			return derivMatU;
		}
		inline const ControlPoints& getDerivMatV() const noexcept {
			return derivMatV;
		}
		inline const ControlPoints& getDerivMatW() const noexcept {
			return derivMatW;
		}

		inline const ControlPoints& getDerivMatUU() const noexcept {
			return derivMatUU;
		}
		inline const ControlPoints& getDerivMatUV() const noexcept {
			return derivMatUV;
		}
		inline const ControlPoints& getDerivMatUW() const noexcept {
			return derivMatUW;
		}
		inline const ControlPoints& getDerivMatVV() const noexcept {
			return derivMatVV;
		}
		inline const ControlPoints& getDerivMatVW() const noexcept {
			return derivMatVW;
		}
		inline const ControlPoints& getDerivMatWW() const noexcept {
			return derivMatWW;
		}

		inline const ControlPoints& getDerivMatUUU() const noexcept {
			return derivMatUUU;
		}
		inline const ControlPoints& getDerivMatUUV() const noexcept {
			return derivMatUUV;
		}
		inline const ControlPoints& getDerivMatUUW() const noexcept {
			return derivMatUUW;
		}
		inline const ControlPoints& getDerivMatUVV() const noexcept {
			return derivMatUVV;
		}
		inline const ControlPoints& getDerivMatUVW() const noexcept {
			return derivMatUVW;
		}
		inline const ControlPoints& getDerivMatUWW() const noexcept {
			return derivMatUWW;
		}
		inline const ControlPoints& getDerivMatVVV() const noexcept {
			return derivMatVVV;
		}
		inline const ControlPoints& getDerivMatVVW() const noexcept {
			return derivMatVVW;
		}
		inline const ControlPoints& getDerivMatVWW() const noexcept {
			return derivMatVWW;
		}
		inline const ControlPoints& getDerivMatWWW() const noexcept {
			return derivMatWWW;
		}
	};
}

#endif