/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_BSPLINE_VOLUME_3D_H__
#define __MN_BSPLINE_VOLUME_3D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"
#include "BezierVolume3d.h"
#include <memory>

namespace MN {
	class BsplineVolume3d : public Freeform3dv {
	public:
		// Bspline volume is divided into a number of patches, which is represented by Bezier volume
		class Patch {
		public:
			// Subdomain : Domain that this patch occupies in original Bspline volume
			Domain uSubdomain = Domain::create(0, 1);
			Domain vSubdomain = Domain::create(0, 1);
			Domain wSubdomain = Domain::create(0, 1);
			BezierVolume3d::Ptr patch = nullptr;

			inline bool domainHas(Real u, Real v, Real w) const noexcept {
				return uSubdomain.has(u) && vSubdomain.has(v) && wSubdomain.has(w);
			}

			// Check if this patch's domain meets given domain or not
			// Assume given domains are not 0-width, and exclude edge cases
			inline bool domainMeet(const Domain& uDomain, const Domain& vDomain, const Domain& wDomain) const noexcept {
				return uSubdomain.has(uDomain) && vSubdomain.has(vDomain) && wSubdomain.has(wDomain);
			}
		};
	private:
		// @direction : 0 for U, 1 for V, 2 for W
		void insertKnot(int direction, Real knot);
		void insertKnotFull(int direction);
		BsplineVolume3d() = default;
	public:
		using Ptr = std::shared_ptr<BsplineVolume3d>;
		using KnotVector = std::vector<Real>;
		// These patches are created in this Bspline volume's creation process, by knot insertion
		std::vector<Patch> patches;
		KnotVector uKnot;
		KnotVector vKnot;
		KnotVector wKnot;

		// @buildMat : Option for building derivMats in creation time
		static BsplineVolume3d create(int uDegree, int vDegree, int wDegree, const KnotVector& uKnot, const KnotVector& vKnot, const KnotVector& wKnot, const ControlPoints& cpts);
		static Ptr createPtr(int uDegree, int vDegree, int wDegree, const KnotVector& uKnot, const KnotVector& vKnot, const KnotVector& wKnot, const ControlPoints& cpts);

		void updatePatches();
		virtual Vec3 evaluate(Real u, Real v, Real w) const;
		virtual Vec3 differentiate(Real u, Real v, Real w, int uOrder, int vOrder, int wOrder) const;
	};
}

#endif