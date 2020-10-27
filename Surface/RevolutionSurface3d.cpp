#include "RevolutionSurface3d.h"

namespace MN {
	RevolutionSurface3d RevolutionSurface3d::create(const Freeform2dc::Ptr profile) {
		RevolutionSurface3d surf;
		surf.setDomain(0, profile->getDomain());
		surf.setDomain(1, Domain::create(0, PI20));
		surf.profile = profile;
		return surf;
	}
	RevolutionSurface3d::Ptr RevolutionSurface3d::createPtr(const Freeform2dc::Ptr profile) {
		return std::make_shared<RevolutionSurface3d>(create(profile));
	}

	static inline Vec3 rotY(const Vec2& v, Real radian) {
		return { v[0] * cos(radian), v[1], v[0] * sin(radian) };
	}
	Vec3 RevolutionSurface3d::evaluate(Real u, Real v) const {
		auto curvePoint = profile->evaluate(u);
		return rotY(curvePoint, v);
	}
	Vec3 RevolutionSurface3d::differentiate(Real u, Real v, int uOrder, int vOrder) const {
		if (uOrder == 0 && vOrder == 0)
			return evaluate(u, v);
		else if (uOrder == 1 && vOrder == 0) {
			auto curveDeriv = profile->differentiate(u, 1);
			return rotY(curveDeriv, v);
		}
		else if (uOrder == 0 && vOrder == 1) {
			auto curvePoint = profile->evaluate(u);
			return { curvePoint[0] * -sin(v), 0, curvePoint[0] * cos(v) };
		}
		else if (uOrder == 2 && vOrder == 0) {
			auto curveDeriv = profile->differentiate(u, 2);
			return rotY(curveDeriv, v);
		}
		else if (uOrder == 1 && vOrder == 1) {
			auto curveDeriv = profile->differentiate(u, 1);
			return { curveDeriv[0] * -sin(v), 0, curveDeriv[0] * cos(v) };
		}
		else if (uOrder == 0 && vOrder == 2) {
			auto curvePoint = profile->evaluate(u);
			return { curvePoint[0] * -cos(v), 0, curvePoint[0] * -sin(v) };
		}
		else if (uOrder == 3 && vOrder == 0) {
			auto curveDeriv = profile->differentiate(u, 3);
			return rotY(curveDeriv, v);
		}
		else if (uOrder == 2 && vOrder == 1) {
			auto curveDeriv = profile->differentiate(u, 2);
			return { curveDeriv[0] * -sin(v), 0, curveDeriv[0] * cos(v) };
		}
		else if (uOrder == 1 && vOrder == 2) {
			auto curveDeriv = profile->differentiate(u, 1);
			return { curveDeriv[0] * -cos(v), 0, curveDeriv[0] * -sin(v) };
		}
		else if (uOrder == 0 && vOrder == 3) {
			auto curvePoint = profile->evaluate(u);
			return { curvePoint[0] * sin(v), 0, curvePoint[0] * -cos(v) };
		}
		else
			throw(std::runtime_error("Bezier surface differentiation is only allowed up to 2nd derivatives"));
	}
}