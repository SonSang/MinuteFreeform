#include "ExtrusionSurface3d.h"

namespace MN {
	ExtrusionSurface3d ExtrusionSurface3d::create(const std::vector<Freeform2dc::Ptr>& profile) {
		ExtrusionSurface3d surface;
		surface.profile = profile;
		return surface;
	}
	ExtrusionSurface3d::Ptr ExtrusionSurface3d::createPtr(const std::vector<Freeform2dc::Ptr>& profile) {
		return std::make_shared<ExtrusionSurface3d>(create(profile));
	}

	const std::vector<Freeform2dc::Ptr>& ExtrusionSurface3d::getProfile() const {
		return profile;
	}
	void ExtrusionSurface3d::setProfile(const std::vector<Freeform2dc::Ptr>& profile) {
		this->profile = profile;
	}

	// u : Parameter for profile curve evaluation
	// v : Parameter for extrusion along Z axis
	Vec3 ExtrusionSurface3d::evaluate(int curveID, Real u, Real v) const {
		Vec2 pt = profile[curveID]->evaluate(u);
		return { pt[0], pt[1], v };
	}
	Vec3 ExtrusionSurface3d::differentiate(int curveID, Real u, Real v, int uOrder, int vOrder) const {
		if (uOrder == 0 && vOrder == 0)
			return evaluate(curveID, u, v);
		else if (uOrder == 1 && vOrder == 0) {
			auto curveDeriv = profile[curveID]->differentiate(u, 1);
			return { curveDeriv[0], curveDeriv[1], 0 };
		}
		else if (uOrder == 0 && vOrder == 1) {
			return { 0, 0, 1 };
		}
		else if (uOrder == 2 && vOrder == 0) {
			auto curveDeriv = profile[curveID]->differentiate(u, 2);
			return { curveDeriv[0], curveDeriv[1], 0 };
		}
		else if (uOrder == 1 && vOrder == 1) {
			return { 0, 0, 0 };
		}
		else if (uOrder == 0 && vOrder == 2) {
			return { 0, 0, 0 };
		}
		else if (uOrder == 3 && vOrder == 0) {
			auto curveDeriv = profile[curveID]->differentiate(u, 3);
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