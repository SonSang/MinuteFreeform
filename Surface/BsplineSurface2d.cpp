#include "BsplineSurface2d.h"
#include <map>

namespace MN {
	bool BsplineSurface2d::Patch::domainHas(Real u, Real v) const noexcept {
		return subdomain.a.has(u) && subdomain.b.has(v);
	}
	bool BsplineSurface2d::Patch::domainHas(const Real2& param) const noexcept {
		return domainHas(param.first, param.second);
	}

	void BsplineSurface2d::insertKnot(int direction, Real knot) {
		auto& knotVector = (direction == 0 ? uKnot : vKnot);
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
			degree = (direction == 0 ? uDegree : vDegree),
			rowSize = (int)cpts.size(),
			colSize = (int)cpts[0].size(),
			nSize = (direction == 0) ? rowSize + 1 : colSize + 1;
		double
			* a = nullptr;
		std::vector<double>
			alpha(nSize, 0.0);

		// 2. Build alpha vector 
		for (int i = 0; i < nSize; i++) {
			a = &alpha[i];
			double
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
		bool
			expandRow = (direction == 1);	// Expand row ?
		int
			nRowSize = (expandRow == true) ? rowSize : rowSize + 1,
			nColSize = (expandRow == true) ? colSize + 1 : colSize;

		ControlPoints nCpts;
		nCpts.resize(nRowSize);
		for (int i = 0; i < nRowSize; i++) {
			nCpts[i].resize(nColSize);
			if (expandRow) {
				for (int j = 0; j < nColSize; j++) {
					Vec2 tmp0 = { 0.0, 0.0 }, tmp1 = { 0.0, 0.0 };
					if (j != colSize)
						tmp0 = cpts[i][j] * alpha[j];
					if (j > 0)
						tmp1 = cpts[i][j - 1] * (1 - alpha[j]);
					nCpts[i][j] = tmp0 + tmp1;
				}
			}
			else
			{
				for (int j = 0; j < nColSize; j++) {
					Vec2 tmp0 = { 0.0, 0.0 }, tmp1 = { 0.0, 0.0 };
					if (i > 0)
						tmp0 = cpts[i - 1][j] * (1 - alpha[i]);
					if (i != rowSize)
						tmp1 = cpts[i][j] * alpha[i];
					nCpts[i][j] = tmp0 + tmp1;
				}
			}
		}
		cpts = nCpts;
	}
	void BsplineSurface2d::insertKnotFull(int direction) {
		auto& knotVector = (direction == 0) ? uKnot : vKnot;
		auto degree = (direction == 0) ? uDegree : vDegree;

		std::map<double, int> insertionTime;
		double beg = knotVector.front();
		double end = knotVector.back();
		double prevKnot = beg;
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
				insertKnot(direction, insertion.first);
		}
	}
	BsplineSurface2d BsplineSurface2d::create(int uDegree, int vDegree, const KnotVector& uKnot, const KnotVector& vKnot, const ControlPoints& cpts) {
		BsplineSurface2d surface;
		surface.uDegree = uDegree;
		surface.vDegree = vDegree;
		surface.uKnot = uKnot;
		surface.uDomain.set(uKnot.front(), uKnot.back());
		surface.vKnot = vKnot;
		surface.vDomain.set(vKnot.front(), vKnot.back());
		surface.cpts = cpts;
		surface.updatePatches();
		return surface;
	}
	BsplineSurface2d::Ptr BsplineSurface2d::createPtr(int uDegree, int vDegree, const KnotVector& uKnot, const KnotVector& vKnot, const ControlPoints& cpts) {
		return std::make_shared<BsplineSurface2d>(create(uDegree, vDegree, uKnot, vKnot, cpts));
	}
	void BsplineSurface2d::updatePatches() {
		// Insert knots in both directions
		insertKnotFull(0);
		insertKnotFull(1);

		patches.clear();

		int uPatchNum = 0;
		int vPatchNum = 0;
		std::vector<double> uniqueKnotsU = { uKnot[0] };
		std::vector<double> uniqueKnotsV = { vKnot[0] };

		double prevKnot = uKnot[0];
		for (auto knot : uKnot) {
			if (knot != prevKnot) {
				uPatchNum++;
				prevKnot = knot;
				uniqueKnotsU.push_back(knot);
			}
		}
		prevKnot = vKnot[0];
		for (auto knot : vKnot) {
			if (knot != prevKnot) {
				vPatchNum++;
				prevKnot = knot;
				uniqueKnotsV.push_back(knot);
			}
		}

		patches.reserve(uPatchNum * vPatchNum);

		int uIndexA = 0, uIndexB = uDegree;
		int vIndexA = 0, vIndexB = vDegree;

		for (int i = 0; i < uPatchNum; i++) {
			Domain uSubdomain = Domain::create(0, 1);
			uSubdomain.set(uniqueKnotsU[i], uniqueKnotsU[i + 1]);

			for (int j = 0; j < vPatchNum; j++) {
				Domain vSubdomain = Domain::create(0, 1);
				vSubdomain.set(uniqueKnotsV[j], uniqueKnotsV[j + 1]);

				ControlPoints bezCpts;
				bezCpts.resize(uDegree + 1);
				for (int m = uIndexA; m <= uIndexB; m++) {
					bezCpts[m - uIndexA].resize(vDegree + 1);
					for (int n = vIndexA; n <= vIndexB; n++)
						bezCpts[m - uIndexA][n - vIndexA] = cpts[m][n];
				}

				Patch patch;
				patch.subdomain.a = uSubdomain;
				patch.subdomain.b = vSubdomain;
				patch.patch = BezierSurface2d::createPtr(uDegree, vDegree, bezCpts);

				patches.push_back(patch);

				vIndexA += vDegree;
				vIndexB += vDegree;
			}
			uIndexA += uDegree;
			uIndexB += uDegree;
			vIndexA = 0;
			vIndexB = vDegree;
		}
	}
	Vec2 BsplineSurface2d::evaluate(Real u, Real v) const {
		for (const auto& patch : patches) {
			if (patch.domainHas(u, v)) {
				double nu = (u - patch.subdomain.a.beg()) / patch.subdomain.a.width();
				double nv = (v - patch.subdomain.b.beg()) / patch.subdomain.b.width();
				return patch.patch->evaluate(nu, nv);
			}
		}
		throw(std::runtime_error("Invalid parameter for Bspline surface evaluation"));
	}
	Vec2 BsplineSurface2d::differentiate(Real u, Real v, int uOrder, int vOrder) const {
		for (const auto& patch : patches) {
			if (patch.domainHas(u, v)) {
				double uWidth, vWidth;
				uWidth = patch.subdomain.a.width();
				vWidth = patch.subdomain.b.width();
				double nu = (u - patch.subdomain.a.beg()) / uWidth;
				double nv = (v - patch.subdomain.b.beg()) / vWidth;
				auto vec = patch.patch->differentiate(nu, nv, uOrder, vOrder);

				if (uOrder == 1 && vOrder == 0)
					vec /= uWidth;
				else if (uOrder == 0 && vOrder == 1)
					vec /= vWidth;
				else if (uOrder == 2 && vOrder == 0)
					vec /= SQ(uWidth);
				else if (uOrder == 1 && vOrder == 1)
					vec /= (uWidth * vWidth);
				else if (uOrder == 0 && vOrder == 2)
					vec /= SQ(vWidth);
				return vec;
			}
		}
		throw(std::runtime_error("Invalid parameter for Bspline surface differentiation"));
	}
}