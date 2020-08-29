/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#include "BezierCurve2d.h"

namespace MN {
	// BezierCurve2d
	const static BezierCurve2d empty = BezierCurve2d::create(0, {});
	void BezierCurve2d::subdivideCpts(const ControlPoints& cpts, Real t, ControlPoints& lower, ControlPoints& upper) {
		// De Casteljou's algorithm
		int size = (int)cpts.size();
		lower.resize(size);
		upper.resize(size);

		double t1 = 1 - t;

		ControlPoints copy = cpts;
		lower[0] = copy.front();
		upper[size - 1] = copy.back();

		int cnt = 1;
		for (int i = size - 1; i > 0; i--) {
			ControlPoints tmp;
			tmp.resize(i);
			for (int j = 0; j < i; j++)
				tmp[j] = copy[j] * t1 + copy[j + 1] * t;
			copy = tmp;

			lower[cnt] = copy.front();
			upper[size - cnt - 1] = copy.back();
			cnt++;
		}
	}
	BezierCurve2d BezierCurve2d::create() {
		return empty;
	}
	BezierCurve2d BezierCurve2d::create(int degree, const ControlPoints& cpts, bool buildMat) {
		BezierCurve2d curve;
		curve.setDegree(degree);
		curve.setCpts(cpts);
		curve.setDomain(Domain::create(0, 1));
		if (buildMat)
			curve.updateDerivMat();
		return curve;
	}
	BezierCurve2d::Ptr BezierCurve2d::createPtr(int degree, const ControlPoints& cpts, bool buildMat) {
		return std::make_shared<BezierCurve2d>(create(degree, cpts, buildMat));
	}
	void BezierCurve2d::updateDerivMat() {
		int num = (int)cpts.size();

		if (degree > 0) {
			// derivMatT
			derivMatT.resize(num - 1);
			for (int i = 0; i < num - 1; i++)
				derivMatT[i] = (cpts[i + 1] - cpts[i]) * degree;
		}
		if (degree > 1) {
			// derivMatTT
			derivMatTT.resize(num - 2);
			for (int i = 0; i < num - 2; i++)
				derivMatTT[i] = (derivMatT[i + 1] - derivMatT[i]) * (degree - 1);
		}
		if (degree > 2) {
			// derivMatTTT
			derivMatTTT.resize(num - 3);
			for (int i = 0; i < num - 3; i++)
				derivMatTTT[i] = (derivMatTT[i + 1] - derivMatTT[i]) * (degree - 2);
		}
	}
	Vec2 BezierCurve2d::evaluate(Real t) const {
		BasisVector basis;
		Bezier::calBasisVector(t, degree, basis);
		return tensorProduct(basis, cpts);
	}
	Vec2 BezierCurve2d::differentiate(Real t, int order) const {
		BasisVector basis;
		if (order == 0)
			return evaluate(t);
		else if (order == 1) {
			if (degree < 1) throw(std::runtime_error("Invalid differerntiation order for bezier curve 2d"));
			Bezier::calBasisVector(t, degree - 1, basis);
			return tensorProduct(basis, derivMatT);
		}
		else if (order == 2) {
			if (degree < 2) throw(std::runtime_error("Invalid differerntiation order for bezier curve 2d"));
			Bezier::calBasisVector(t, degree - 2, basis);
			return tensorProduct(basis, derivMatTT);
		}
		else if (order == 3) {
			if (degree < 3) throw(std::runtime_error("Invalid differerntiation order for bezier curve 2d"));
			Bezier::calBasisVector(t, degree - 3, basis);
			return tensorProduct(basis, derivMatTTT);
		}
		else
			throw(std::runtime_error("Bezier curve differentiation is only allowed up to 3rd derivatives"));
	}
	void BezierCurve2d::subdivide(Real t, BezierCurve2d& lower, BezierCurve2d& upper) const {
		ControlPoints lowerCpts, upperCpts;
		subdivideCpts(cpts, t, lowerCpts, upperCpts);

		lower = BezierCurve2d::create(degree, lowerCpts, false);
		upper = BezierCurve2d::create(degree, upperCpts, false);
	}
	BezierCurve2d::Ptr BezierCurve2d::subdivide(const Domain& subdomain) const {
		BezierCurve2d tmpLower, tmpUpper;
		subdivide(subdomain.beg(), tmpLower, tmpUpper);

		Real t = (subdomain.width()) / (1.0 - subdomain.beg());
		tmpUpper.subdivide(t, tmpLower, tmpUpper);
		tmpLower.updateDerivMat();

		return std::make_shared<BezierCurve2d>(tmpLower);
	}
}