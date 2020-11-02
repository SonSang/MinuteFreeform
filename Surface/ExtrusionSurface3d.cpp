#include "ExtrusionSurface3d.h"

namespace MN {
	ExtrusionSurface3d ExtrusionSurface3d::create(const Freeform2dc::Ptr profile) {
		ExtrusionSurface3d surface;
		surface.profile = profile;
		surface.uDomain = profile->getDomain();
		surface.vDomain.set(0, 1);
		return surface;
	}
	ExtrusionSurface3d::Ptr ExtrusionSurface3d::createPtr(const Freeform2dc::Ptr profile) {
		return std::make_shared<ExtrusionSurface3d>(create(profile));
	}

	Freeform2dc::Ptr ExtrusionSurface3d::getProfile() const {
		return profile;
	}
	void ExtrusionSurface3d::setProfile(const Freeform2dc::Ptr profile) {
		this->profile = profile;
		uDomain = profile->getDomain();
	}

	Vec3 ExtrusionSurface3d::evaluate(Real u, Real v) const {
		Vec2 pt = profile->evaluate(u);
		return { pt[0], pt[1], v };
	}
	Vec3 ExtrusionSurface3d::differentiate(Real u, Real v, int uOrder, int vOrder) const {
		if (uOrder == 0 && vOrder == 0)
			return evaluate(u, v);
		else if (uOrder == 1 && vOrder == 0) {
			auto curveDeriv = profile->differentiate(u, 1);
			return { curveDeriv[0], curveDeriv[1], 0 };
		}
		else if (uOrder == 0 && vOrder == 1) {
			return { 0, 0, 1 };
		}
		else if (uOrder == 2 && vOrder == 0) {
			auto curveDeriv = profile->differentiate(u, 2);
			return { curveDeriv[0], curveDeriv[1], 0 };
		}
		else if (uOrder == 1 && vOrder == 1) {
			return { 0, 0, 0 };
		}
		else if (uOrder == 0 && vOrder == 2) {
			return { 0, 0, 0 };
		}
		else if (uOrder == 3 && vOrder == 0) {
			auto curveDeriv = profile->differentiate(u, 3);
			return { curveDeriv[0], curveDeriv[1], 0 };
		}
		else if (uOrder == 2 && vOrder == 1) {
			return { 0, 0, 0 };
		}
		else if (uOrder == 1 && vOrder == 2) {
			return { 0, 0, 0 };
		}
		else if (uOrder == 0 && vOrder == 3) {
			return { 0, 0, 0 };
		}
		else
			throw(std::runtime_error("Extrusion surface differentiation is only allowed up to 2nd derivatives"));
	}
}