#ifndef __MN_BSPLINE_SURFACE_3D_H__
#define __MN_BSPLINE_SURFACE_3D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"
#include "BezierSurface3d.h"
#include <vector>

namespace MN {
	class BsplineSurface3d : public Freeform3ds {
	public:
		// Bspline surface is divided into a number of patches, which is represented by Bezier surface
		class Patch {
		public:
			// Subdomain : Domain that this patch occupies in original Bspline surface
			Domain uSubdomain = Domain::create(0, 1);
			Domain vSubdomain = Domain::create(0, 1);
			BezierSurface3d::Ptr patch = nullptr;

			bool domainHas(double u, double v) const noexcept;

			// Check if this patch's domain meets given domain or not
			// Assume given domains are not 0-width, and exclude edge cases
			bool domainMeet(Domain& uDomain, Domain& vDomain) const noexcept;
		};
	private:
		void insertKnot(int direction, double knot);
		void insertKnotFull(int direction);
	public:
		using Ptr = std::shared_ptr<BsplineSurface3d>;
		using KnotVector = std::vector<double>;
		// These patches are created in this Bspline surface's creation process, by knot insertion
		std::vector<Patch> patches;
		KnotVector uKnot;
		KnotVector vKnot;

		static BsplineSurface3d create(int uDegree, int vDegree, const KnotVector& uKnot, const KnotVector& vKnot, const ControlPoints& cpts);
		static Ptr createPtr(int uDegree, int vDegree, const KnotVector& uKnot, const KnotVector& vKnot, const ControlPoints& cpts);

		void updatePatches();
		virtual Vec3 evaluate(double u, double v) const;
		virtual Vec3 differentiate(double u, double v, int uOrder, int vOrder) const;
	};
}

#endif