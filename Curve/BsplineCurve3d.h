/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_BSPLINE_CURVE_3D_H__
#define __MN_BSPLINE_CURVE_3D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"
#include "BezierCurve3d.h"

namespace MN {
	class BsplineCurve3d : public Freeform3dc {
	public:
		// Bspline curve is divided into a number of patches, which is represented by Bezier curve
		class Patch {
		public:
			// Subdomain : Domain that this patch occupies in original Bspline curve
			Domain subdomain = Domain::create(0, 1);
			BezierCurve3d::Ptr curve = nullptr;
		};
	public:
		using Ptr = std::shared_ptr<BsplineCurve3d>;
		using PatchVector = std::vector<Patch>;

		KnotVector knotVector;
		PatchVector patchVector;
	private:
		void insertKnot(Real knot);
		void insertKnotFull();
	public:
		static BsplineCurve3d create(int degree, const KnotVector& knot, const ControlPoints& cpts);
		static Ptr createPtr(int degree, const KnotVector& knot, const ControlPoints& cpts);

		void setKnotVector(const KnotVector& knotVector) noexcept;
		KnotVector& getKnotVector() noexcept;
		const KnotVector& getKnotVectorC() const noexcept;

		PatchVector& getPatchVector() noexcept;
		const PatchVector& getPatchVectorC() const noexcept;

		virtual Vec3 evaluate(Real t) const;
		virtual Vec3 differentiate(Real t, int order) const;

		void updatePatches();
	};
}

#endif