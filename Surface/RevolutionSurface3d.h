/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_REVOLUTION_SURFACE_3D_H__
#define __MN_REVOLUTION_SURFACE_3D_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Freeform.h"
#include "../Curve/BsplineCurve2d.h"
#include <memory>

namespace MN {
	class RevolutionSurface3d : public Freeform3ds {
	private:
		RevolutionSurface3d() = default;

		// Profile curve on XY plane that is rotated around Y axis to form surface
		Freeform2dc::Ptr profile;
	public:
		using Ptr = std::shared_ptr<RevolutionSurface3d>;

		// @buildMat : Option for building derivMats in creation time
		static RevolutionSurface3d create(const Freeform2dc::Ptr profile);
		static Ptr createPtr(const Freeform2dc::Ptr profile);

		Freeform2dc::Ptr getProfile() const;
		void setProfile(const Freeform2dc::Ptr profile);

		// u : Parameter for profile curve evaluation
		// v : Parameter for rotation around Y axis
		virtual Vec3 evaluate(Real u, Real v) const;
		virtual Vec3 differentiate(Real u, Real v, int uOrder, int vOrder) const;
	};
}

#endif