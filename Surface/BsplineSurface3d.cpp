#include "BsplineSurface3d.h"
#include <map>

namespace MN {
	// BsplineSurface Patch
	bool BsplineSurface3d::Patch::domainHas(double u, double v) const noexcept {
		return uSubdomain.has(u) && vSubdomain.has(v);
	}
	bool BsplineSurface3d::Patch::domainMeet(Domain& uDomain, Domain& vDomain) const noexcept {
		bool uMeet = false, vMeet = false;
		// uDomain
		if (uSubdomain.beg() <= uDomain.beg() && uSubdomain.end() > uDomain.beg()) {
			uMeet = true;
			if (uSubdomain.end() <= uDomain.end())
				uDomain.set(uDomain.beg(), uSubdomain.end());
		}
		else if (uSubdomain.beg() < uDomain.end() && uSubdomain.end() >= uDomain.beg()) {
			uMeet = true;
			uDomain.set(uSubdomain.beg(), uDomain.end());
		}
		else if (uSubdomain.beg() > uDomain.beg() && uSubdomain.end() <= uDomain.end()) {
			uMeet = true;
			uDomain = uSubdomain;
		}
		// vDomain
		if (vSubdomain.beg() <= vDomain.beg() && vSubdomain.end() > vDomain.beg()) {
			vMeet = true;
			if (vSubdomain.end() <= vDomain.end())
				vDomain.set(vDomain.beg(), vSubdomain.end());
		}
		else if (vSubdomain.beg() < vDomain.end() && vSubdomain.end() >= vDomain.end()) {
			vMeet = true;
			vDomain.set(vSubdomain.beg(), vDomain.end());
		}
		else if (vSubdomain.beg() > vDomain.beg() && vSubdomain.end() <= vDomain.end()) {
			vMeet = true;
			vDomain = vSubdomain;
		}
		return (uMeet && vMeet);
	}
	// BsplineSurface
	void BsplineSurface3d::insertKnot(int direction, double knot) {
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
					Vec3 tmp0 = { 0.0, 0.0, 0.0 }, tmp1 = { 0.0, 0.0, 0.0 };
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
					Vec3 tmp0 = { 0.0, 0.0, 0.0 }, tmp1 = { 0.0, 0.0, 0.0 };
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
	void BsplineSurface3d::insertKnotFull(int direction) {
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
	BsplineSurface3d BsplineSurface3d::create(int uDegree, int vDegree, const KnotVector& uKnot, const KnotVector& vKnot, const ControlPoints& cpts) {
		BsplineSurface3d surface;
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
	BsplineSurface3d::Ptr BsplineSurface3d::createPtr(int uDegree, int vDegree, const KnotVector& uKnot, const KnotVector& vKnot, const ControlPoints& cpts) {
		BsplineSurface3d surface = create(uDegree, vDegree, uKnot, vKnot, cpts);
		return std::make_shared<BsplineSurface3d>(surface);
	}
	void BsplineSurface3d::updatePatches() {
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
			Domain uSubdomain = Domain::create(uniqueKnotsU[i], uniqueKnotsU[i + 1]);

			for (int j = 0; j < vPatchNum; j++) {
				Domain vSubdomain = Domain::create(uniqueKnotsV[j], uniqueKnotsV[j + 1]);

				ControlPoints bezCpts;
				bezCpts.resize(uDegree + 1);
				for (int m = uIndexA; m <= uIndexB; m++) {
					bezCpts[m - uIndexA].resize(vDegree + 1);
					for (int n = vIndexA; n <= vIndexB; n++)
						bezCpts[m - uIndexA][n - vIndexA] = cpts[m][n];
				}

				Patch patch;
				patch.uSubdomain = uSubdomain;
				patch.vSubdomain = vSubdomain;
				patch.patch = BezierSurface3d::createPtr(uDegree, vDegree, bezCpts);

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
	Vec3 BsplineSurface3d::evaluate(double u, double v) const {
		for (const auto& patch : patches) {
			if (patch.domainHas(u, v)) {
				double nu = (u - patch.uSubdomain.beg()) / patch.uSubdomain.width();
				double nv = (v - patch.vSubdomain.beg()) / patch.vSubdomain.width();
				return patch.patch->evaluate(nu, nv);
			}
		}
		throw(std::runtime_error("Invalid parameter for Bspline surface evaluation"));
	}
	Vec3 BsplineSurface3d::differentiate(double u, double v, int uOrder, int vOrder) const {
		for (const auto& patch : patches) {
			if (patch.domainHas(u, v)) {
				double uWidth, vWidth;
				uWidth = patch.uSubdomain.width();
				vWidth = patch.vSubdomain.width();
				double nu = (u - patch.uSubdomain.beg()) / uWidth;
				double nv = (v - patch.vSubdomain.beg()) / vWidth;
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