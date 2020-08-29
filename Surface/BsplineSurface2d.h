/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_BSPLINE_SURFACE_2D_H__
#define __MN_BSPLINE_SURFACE_2D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"
#include "BezierSurface2d.h"
#include <vector>

namespace MN {
	// 2d
	class BsplineSurface2d : public Freeform2ds {
	public:
		// Bspline surface 2d is divided into a number of patches, which is represented by Bezier surface 2d
		class Patch {
		public:
			// Subdomain : Domain that this patch occupies in original Bspline surface
			Domain2 subdomain;
			BezierSurface2d::Ptr patch = nullptr;

			bool domainHas(Real u, Real v) const noexcept;
			bool domainHas(const Real2& param) const noexcept;
		};
	private:
		BsplineSurface2d() = default;

		void insertKnot(int direction, Real knot);
		void insertKnotFull(int direction);
	public:
		using Ptr = std::shared_ptr<BsplineSurface2d>;
		// These patches are created in this Bspline surface's creation process, by knot insertion
		std::vector<Patch> patches;
		KnotVector uKnot;
		KnotVector vKnot;

		static BsplineSurface2d create(int uDegree, int vDegree, const KnotVector& uKnot, const KnotVector& vKnot, const ControlPoints& cpts);
		static Ptr createPtr(int uDegree, int vDegree, const KnotVector& uKnot, const KnotVector& vKnot, const ControlPoints& cpts);

		void updatePatches();
		virtual Vec2 evaluate(Real u, Real v) const;
		virtual Vec2 differentiate(Real u, Real v, int uOrder, int vOrder) const;
	};
}

#endif