/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#include "BsplineCurve3d.h"
#include <map>

namespace MN {
	void BsplineCurve3d::insertKnot(Real knot) {
		int
			k = 0,
			knotSize = (int)knotVector.size();

		// 1. Find [k] such that holds [ knot ] in between : knotVector[k] <= knot < knotVector[k+1]
		k = -1;
		for (int i = 0; i < knotSize; i++)
		{
			if (i == knotSize - 1)
				break;
			if (knotVector[i] <= knot && knotVector[i + 1] > knot)
				k = i;
		}
		if (k == -1)
			throw(std::runtime_error("Invalid knot value for knot insertion"));

		int
			nSize = (int)cpts.size() + 1;
		Real
			* a = nullptr;
		std::vector<Real>
			alpha(nSize, 0.0);

		// 2. Build alpha vector 
		for (int i = 0; i < nSize; i++) {
			a = &alpha[i];
			Real
				kvA = knotVector[i],
				kvB = knotVector[i + degree];
			if (i <= k - degree)
				*a = 1.0;
			else if (i >= k + 1)
				*a = 0.0;
			else
				*a = (knot - kvA) / (kvB - kvA);
		}

		// 3. Create new knot vector
		int
			nKnotSize = knotSize + 1;
		KnotVector
			nKnotVector;
		for (int i = 0; i < nKnotSize; i++) {
			if (i < k + 1)
				nKnotVector.push_back(knotVector[i]);
			else if (i == k + 1)
				nKnotVector.push_back(knot);
			else
				nKnotVector.push_back(knotVector[i - 1]);
		}
		knotVector = nKnotVector;

		// 4. Create new control points
		ControlPoints nCpts;
		nCpts.resize(nSize);
		for (int i = 0; i < nSize; i++) {
			Vec3 tmp0 = { 0.0, 0.0, 0.0 }, tmp1 = { 0.0, 0.0, 0.0 };
			if (i != nSize - 1)
				tmp0 = cpts[i] * alpha[i];
			if (i > 0)
				tmp1 = cpts[i - 1] * (1 - alpha[i]);
			nCpts[i] = tmp0 + tmp1;
		}
		cpts = nCpts;
	}
	void BsplineCurve3d::insertKnotFull() {
		std::map<Real, int> insertionTime;
		Real beg = knotVector.front();
		Real end = knotVector.back();
		Real prevKnot = beg;
		int multiplicity = 0;
		for (int i = 1; i < knotVector.size(); i++) {
			if (knotVector[i] == prevKnot)
				multiplicity++;
			else {
				if (prevKnot == beg)
					insertionTime.insert({ prevKnot, degree - multiplicity });
				else
					insertionTime.insert({ prevKnot, degree - 1 - multiplicity });
				multiplicity = 0;
			}
			prevKnot = knotVector[i];
		}

		for (auto& insertion : insertionTime) {
			for (int i = 0; i < insertion.second; i++)
				insertKnot(insertion.first);
		}
	}
	BsplineCurve3d BsplineCurve3d::create(int degree, const KnotVector& knot, const ControlPoints& cpts) {
		BsplineCurve3d curve;
		curve.setDegree(degree);
		curve.setKnotVector(knot);
		curve.setDomain(Domain::create(knot.front(), knot.back()));
		curve.setCpts(cpts);
		curve.updatePatches();
		return curve;
	}
	BsplineCurve3d::Ptr BsplineCurve3d::createPtr(int degree, const KnotVector& knot, const ControlPoints& cpts) {
		return std::make_shared<BsplineCurve3d>(create(degree, knot, cpts));
	}
	void BsplineCurve3d::setKnotVector(const KnotVector& knotVector) noexcept {
		this->knotVector = knotVector;
	}
	KnotVector& BsplineCurve3d::getKnotVector() noexcept {
		return knotVector;
	}
	const KnotVector& BsplineCurve3d::getKnotVectorC() const noexcept {
		return knotVector;
	}
	BsplineCurve3d::PatchVector& BsplineCurve3d::getPatchVector() noexcept {
		return patchVector;
	}
	const BsplineCurve3d::PatchVector& BsplineCurve3d::getPatchVectorC() const noexcept {
		return patchVector;
	}
	Vec3 BsplineCurve3d::evaluate(Real t) const {
		for (const auto& patch : patchVector) {
			if (patch.subdomain.has(t)) {
				Real nt = (t - patch.subdomain.beg()) / patch.subdomain.width();
				return patch.curve->evaluate(nt);
			}
		}
		throw(std::runtime_error("Invalid parameter for Bspline curve 2d evaluation"));
	}
	Vec3 BsplineCurve3d::differentiate(Real t, int order) const {
		for (const auto& patch : patchVector) {
			if (patch.subdomain.has(t)) {
				Real width;
				width = patch.subdomain.width();
				double nt = (t - patch.subdomain.beg()) / width;
				auto vec = patch.curve->differentiate(nt, order);

				if (order == 1)
					vec /= width;
				else if (order == 2)
					vec /= (width * width);
				else if (order == 3)
					vec /= (width * width * width);
				return vec;
			}
		}
		throw(std::runtime_error("Invalid parameter for Bspline curve 2d differentiation"));
	}
	void BsplineCurve3d::updatePatches() {
		// Insert knots full
		insertKnotFull();

		patchVector.clear();

		int patchNum = 0;
		std::vector<Real> uniqueKnots = { knotVector[0] };

		Real prevKnot = knotVector[0];
		for (auto knot : knotVector) {
			if (knot != prevKnot) {
				patchNum++;
				prevKnot = knot;
				uniqueKnots.push_back(knot);
			}
		}
		patchVector.reserve(patchNum);

		int indexA = 0, indexB = degree;
		for (int i = 0; i < patchNum; i++) {
			Domain subdomain = Domain::create(uniqueKnots[i], uniqueKnots[i + 1]);

			ControlPoints bezCpts;
			bezCpts.resize(degree + 1);
			for (int m = indexA; m <= indexB; m++)
				bezCpts[m - indexA] = cpts[m];

			Patch patch;
			patch.subdomain = subdomain;
			patch.curve = BezierCurve3d::createPtr(degree, bezCpts);

			patchVector.push_back(patch);
			indexA += degree;
			indexB += degree;
		}
	}
}