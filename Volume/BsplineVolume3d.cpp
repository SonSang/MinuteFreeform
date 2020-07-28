#include "BsplineVolume3d.h"
#include <map>

namespace MN {
	// BsplineVolume3d
	void BsplineVolume3d::insertKnot(int direction, Real knot) {
		auto& knotVector = (direction == 0 ? uKnot : (direction == 1 ? vKnot : wKnot));
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
			degree = (direction == 0 ? uDegree : (direction == 1 ? vDegree : wDegree)),
			uSize = (int)cpts.size(),
			vSize = (int)cpts[0].size(),
			wSize = (int)cpts[0][0].size(),
			nSize = (direction == 0 ? uSize + 1 : (direction == 1 ? vSize + 1 : wSize + 1));
		Real
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
		int
			nuSize = uSize,
			nvSize = vSize,
			nwSize = wSize;
		if (direction == 0)
			nuSize++;
		else if (direction == 1)
			nvSize++;
		else
			nwSize++;

		ControlPoints nCpts;
		nCpts.resize(nuSize);
		for (int i = 0; i < nuSize; i++) {
			nCpts[i].resize(nvSize);
			for (int j = 0; j < nvSize; j++) {
				nCpts[i][j].resize(nwSize);
			}
		}

		if (direction == 0) {
			for (int j = 0; j < nvSize; j++) {
				for (int k = 0; k < nwSize; k++) {
					for (int i = 0; i < nuSize; i++) {
						Vec3 tmp0 = { 0.0, 0.0, 0.0 }, tmp1 = { 0.0, 0.0, 0.0 };
						if (i < uSize)
							tmp0 = cpts[i][j][k] * alpha[i];
						if (i > 0)
							tmp1 = cpts[i - 1][j][k] * (1 - alpha[i]);
						nCpts[i][j][k] = tmp0 + tmp1;
					}
				}
			}
		}
		else if (direction == 1) {
			for (int i = 0; i < nuSize; i++) {
				for (int k = 0; k < nwSize; k++) {
					for (int j = 0; j < nvSize; j++) {
						Vec3 tmp0 = { 0.0, 0.0, 0.0 }, tmp1 = { 0.0, 0.0, 0.0 };
						if (j < vSize)
							tmp0 = cpts[i][j][k] * alpha[j];
						if (j > 0)
							tmp1 = cpts[i][j - 1][k] * (1 - alpha[j]);
						nCpts[i][j][k] = tmp0 + tmp1;
					}
				}
			}
		}
		else {
			for (int i = 0; i < nuSize; i++) {
				for (int j = 0; j < nvSize; j++) {
					for (int k = 0; k < nwSize; k++) {
						Vec3 tmp0 = { 0.0, 0.0, 0.0 }, tmp1 = { 0.0, 0.0, 0.0 };
						if (k < wSize)
							tmp0 = cpts[i][j][k] * alpha[k];
						if (k > 0)
							tmp1 = cpts[i][j][k - 1] * (1 - alpha[k]);
						nCpts[i][j][k] = tmp0 + tmp1;
					}
				}
			}
		}
		cpts = nCpts;
	}
	void BsplineVolume3d::insertKnotFull(int direction) {
		auto& knotVector = ((direction == 0) ? uKnot : (direction == 1 ? vKnot : wKnot));
		auto degree = ((direction == 0) ? uDegree : (direction == 1 ? vDegree : wDegree));

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
				insertKnot(direction, insertion.first);
		}
	}
	BsplineVolume3d BsplineVolume3d::create(int uDegree, int vDegree, int wDegree, const KnotVector& uKnot, const KnotVector& vKnot, const KnotVector& wKnot, const ControlPoints& cpts) {
		BsplineVolume3d volume;
		volume.uDegree = uDegree;
		volume.vDegree = vDegree;
		volume.wDegree = wDegree;

		volume.uKnot = uKnot;
		volume.uDomain.set(uKnot.front(), uKnot.back());
		volume.vKnot = vKnot;
		volume.vDomain.set(vKnot.front(), vKnot.back());
		volume.wKnot = wKnot;
		volume.wDomain.set(wKnot.front(), wKnot.back());

		volume.cpts = cpts;
		volume.updatePatches();
		return volume;
	}
	BsplineVolume3d::Ptr BsplineVolume3d::createPtr(int uDegree, int vDegree, int wDegree, const KnotVector& uKnot, const KnotVector& vKnot, const KnotVector& wKnot, const ControlPoints& cpts) {
		return std::make_shared<BsplineVolume3d>(create(uDegree, vDegree, wDegree, uKnot, vKnot, wKnot, cpts));
	}
	void BsplineVolume3d::updatePatches() {
		// Insert knots in all directions
		insertKnotFull(0);
		insertKnotFull(1);
		insertKnotFull(2);

		patches.clear();

		int uPatchNum = 0;
		int vPatchNum = 0;
		int wPatchNum = 0;
		std::vector<double> uniqueKnotsU = { uKnot[0] };
		std::vector<double> uniqueKnotsV = { vKnot[0] };
		std::vector<double> uniqueKnotsW = { wKnot[0] };

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
		prevKnot = wKnot[0];
		for (auto knot : wKnot) {
			if (knot != prevKnot) {
				wPatchNum++;
				prevKnot = knot;
				uniqueKnotsW.push_back(knot);
			}
		}

		patches.reserve(uPatchNum * vPatchNum * wPatchNum);

		int uIndexA = 0, uIndexB = uDegree;
		int vIndexA = 0, vIndexB = vDegree;
		int wIndexA = 0, wIndexB = wDegree;

		for (int i = 0; i < uPatchNum; i++) {
			Domain uSubdomain = Domain::create(uniqueKnotsU[i], uniqueKnotsU[i + 1]);

			for (int j = 0; j < vPatchNum; j++) {
				Domain vSubdomain = Domain::create(uniqueKnotsV[j], uniqueKnotsV[j + 1]);

				for (int k = 0; k < wPatchNum; k++) {
					Domain wSubdomain = Domain::create(uniqueKnotsW[k], uniqueKnotsW[k + 1]);

					ControlPoints bezCpts;
					bezCpts.resize(uDegree + 1);

					for (int p = uIndexA; p <= uIndexB; p++) {
						bezCpts[p - uIndexA].resize(vDegree + 1);
						for (int q = vIndexA; q <= vIndexB; q++) {
							bezCpts[p - uIndexA][q - vIndexA].resize(wDegree + 1);
							for (int r = wIndexA; r <= wIndexB; r++) {
								bezCpts[p - uIndexA][q - vIndexA][r - wIndexA] = cpts[p][q][r];
							}
						}
					}

					Patch patch;
					patch.uSubdomain = uSubdomain;
					patch.vSubdomain = vSubdomain;
					patch.wSubdomain = wSubdomain;
					patch.patch = BezierVolume3d::createPtr(uDegree, vDegree, wDegree, bezCpts);

					patches.push_back(patch);

					wIndexA += wDegree;
					wIndexB += wDegree;
				}
				vIndexA += vDegree;
				vIndexB += vDegree;
				wIndexA = 0;
				wIndexB = wDegree;
			}
			uIndexA += uDegree;
			uIndexB += uDegree;
			vIndexA = 0;
			vIndexB = vDegree;
			wIndexA = 0;
			wIndexB = wDegree;
		}
	}

	Vec3 BsplineVolume3d::evaluate(Real u, Real v, Real w) const {
		for (const auto& patch : patches) {
			if (patch.domainHas(u, v, w)) {
				Real nu = (u - patch.uSubdomain.beg()) / patch.uSubdomain.width();
				Real nv = (v - patch.vSubdomain.beg()) / patch.vSubdomain.width();
				Real nw = (w - patch.wSubdomain.beg()) / patch.wSubdomain.width();
				return patch.patch->evaluate(nu, nv, nw);
			}
		}
		throw(std::runtime_error("Invalid parameter for Bspline surface evaluation"));
	}
	Vec3 BsplineVolume3d::differentiate(Real u, Real v, Real w, int uOrder, int vOrder, int wOrder) const {
		for (const auto& patch : patches) {
			if (patch.domainHas(u, v, w)) {
				Real uWidth, vWidth, wWidth;
				uWidth = patch.uSubdomain.width();
				vWidth = patch.vSubdomain.width();
				wWidth = patch.wSubdomain.width();
				Real nu = (u - patch.uSubdomain.beg()) / uWidth;
				Real nv = (v - patch.vSubdomain.beg()) / vWidth;
				Real nw = (w - patch.wSubdomain.beg()) / wWidth;
				auto vec = patch.patch->differentiate(nu, nv, nw, uOrder, vOrder, wOrder);

				if (uOrder == 1 && vOrder == 0 && wOrder == 0)
					vec /= uWidth;
				else if (uOrder == 0 && vOrder == 1 && wOrder == 0)
					vec /= vWidth;
				else if (uOrder == 0 && vOrder == 0 && wOrder == 1)
					vec /= wWidth;

				else if (uOrder == 2 && vOrder == 0 && wOrder == 0)
					vec /= SQ(uWidth);
				else if (uOrder == 1 && vOrder == 1 && wOrder == 0)
					vec /= (uWidth * vWidth);
				else if (uOrder == 1 && vOrder == 0 && wOrder == 1)
					vec /= (uWidth * wWidth);
				else if (uOrder == 0 && vOrder == 2 && wOrder == 0)
					vec /= SQ(vWidth);
				else if (uOrder == 0 && vOrder == 1 && wOrder == 1)
					vec /= (vWidth * wWidth);
				else if (uOrder == 0 && vOrder == 0 && wOrder == 2)
					vec /= SQ(wWidth);

				else if (uOrder == 3 && vOrder == 0 && wOrder == 0)
					vec /= (uWidth * uWidth * uWidth);
				else if (uOrder == 2 && vOrder == 1 && wOrder == 0)
					vec /= (uWidth * uWidth * vWidth);
				else if (uOrder == 2 && vOrder == 0 && wOrder == 1)
					vec /= (uWidth * uWidth * wWidth);
				else if (uOrder == 1 && vOrder == 2 && wOrder == 0)
					vec /= (uWidth * vWidth * vWidth);
				else if (uOrder == 1 && vOrder == 1 && wOrder == 1)
					vec /= (uWidth * vWidth * wWidth);
				else if (uOrder == 1 && vOrder == 0 && wOrder == 2)
					vec /= (uWidth * wWidth * wWidth);
				else if (uOrder == 0 && vOrder == 3 && wOrder == 0)
					vec /= (vWidth * vWidth * vWidth);
				else if (uOrder == 0 && vOrder == 2 && wOrder == 1)
					vec /= (vWidth * vWidth * wWidth);
				else if (uOrder == 0 && vOrder == 1 && wOrder == 2)
					vec /= (vWidth * wWidth * wWidth);
				else if (uOrder == 0 && vOrder == 0 && wOrder == 3)
					vec /= (wWidth * wWidth * wWidth);

				return vec;
			}
		}
		throw(std::runtime_error("Invalid parameter for Bspline surface differentiation"));
	}
}