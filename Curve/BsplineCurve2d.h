/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_BSPLINE_CURVE_2D_H__
#define __MN_BSPLINE_CURVE_2D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"
#include "BezierCurve2d.h"

namespace MN {
	class BsplineCurve2d : public Freeform2dc {
	public:
		// Bspline curve is divided into a number of patches, which is represented by Bezier curves
		class Patch {
		public:
			// Subdomain : Domain that this patch occupies in original Bspline curve
			Domain subdomain = Domain::create(0, 1);
			BezierCurve2d::Ptr curve = nullptr;
		};
	public:
		using Ptr = std::shared_ptr<BsplineCurve2d>;
		using PatchVector = std::vector<Patch>;

		KnotVector knotVector;
		PatchVector patchVector;
	private:
		void insertKnot(Real knot);
		void insertKnotFull();
	public:
		static BsplineCurve2d create(int degree, const KnotVector& knot, const ControlPoints& cpts);
		static Ptr createPtr(int degree, const KnotVector& knot, const ControlPoints& cpts);

		void setKnotVector(const KnotVector& knotVector) noexcept;
		KnotVector& getKnotVector() noexcept;
		const KnotVector& getKnotVectorC() const noexcept;

		PatchVector& getPatchVector() noexcept;
		const PatchVector& getPatchVectorC() const noexcept;

		virtual Vec2 evaluate(Real t) const;
		virtual Vec2 differentiate(Real t, int order) const;

		void updatePatches();
	};
}

#endif