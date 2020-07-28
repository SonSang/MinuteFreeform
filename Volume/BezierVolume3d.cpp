#include "BezierVolume3d.h"

namespace MN {
	// BezierVolume3d
	static const BezierVolume3d empty = BezierVolume3d::create(0, 0, 0, {}, false);

	void BezierVolume3d::subdivideCpts(const std::vector<Vec3>& cpts, Real t, std::vector<Vec3>& lower, std::vector<Vec3>& upper) {
		// De Casteljou's algorithm
		int size = (int)cpts.size();
		lower.resize(size);
		upper.resize(size);

		Real t1 = 1 - t;

		std::vector<Vec3> copy = cpts;
		lower[0] = copy.front();
		upper[size - 1] = copy.back();

		int cnt = 1;
		for (int i = size - 1; i > 0; i--) {
			std::vector<Vec3> tmp;
			tmp.resize(i);
			for (int j = 0; j < i; j++)
				tmp[j] = copy[j] * t1 + copy[j + 1] * t;
			copy = tmp;

			lower[cnt] = copy.front();
			upper[size - cnt - 1] = copy.back();
			cnt++;
		}
	}
	void BezierVolume3d::uSubdivide(Real u, BezierVolume3d& lower, BezierVolume3d& upper, bool buildMat) const {
		ControlPoints lowerCpts = cpts;
		ControlPoints upperCpts = cpts;
		int uSize = (int)cpts.size();
		int vSize = (int)cpts[0].size();
		int wSize = (int)cpts[0][0].size();

		std::vector<Vec3> lowerCol;
		std::vector<Vec3> upperCol;
		for (int i = 0; i < vSize; i++) {
			for (int j = 0; j < wSize; j++) {
				std::vector<Vec3> col;
				col.resize(uSize);
				for (int k = 0; k < uSize; k++)
					col[k] = cpts[k][i][j];
				subdivideCpts(col, u, lowerCol, upperCol);

				for (int k = 0; k < uSize; k++) {
					lowerCpts[k][i][j] = lowerCol[k];
					upperCpts[k][i][j] = upperCol[k];
				}
			}
		}

		lower = BezierVolume3d::create(uDegree, vDegree, wDegree, lowerCpts, buildMat);
		upper = BezierVolume3d::create(uDegree, vDegree, wDegree, upperCpts, buildMat);
	}
	void BezierVolume3d::vSubdivide(Real v, BezierVolume3d& lower, BezierVolume3d& upper, bool buildMat) const {
		ControlPoints lowerCpts = cpts;
		ControlPoints upperCpts = cpts;
		int uSize = (int)cpts.size();
		int vSize = (int)cpts[0].size();
		int wSize = (int)cpts[0][0].size();

		std::vector<Vec3> lowerCol;
		std::vector<Vec3> upperCol;
		for (int i = 0; i < uSize; i++) {
			for (int j = 0; j < wSize; j++) {
				std::vector<Vec3> col;
				col.resize(vSize);
				for (int k = 0; k < vSize; k++)
					col[k] = cpts[i][k][j];
				subdivideCpts(col, v, lowerCol, upperCol);

				for (int k = 0; k < vSize; k++) {
					lowerCpts[i][k][j] = lowerCol[k];
					upperCpts[i][k][j] = upperCol[k];
				}
			}
		}

		lower = BezierVolume3d::create(uDegree, vDegree, wDegree, lowerCpts, buildMat);
		upper = BezierVolume3d::create(uDegree, vDegree, wDegree, upperCpts, buildMat);
	}
	void BezierVolume3d::wSubdivide(Real w, BezierVolume3d& lower, BezierVolume3d& upper, bool buildMat) const {
		ControlPoints lowerCpts = cpts;
		ControlPoints upperCpts = cpts;
		int uSize = (int)cpts.size();
		int vSize = (int)cpts[0].size();
		int wSize = (int)cpts[0][0].size();

		std::vector<Vec3> lowerCol;
		std::vector<Vec3> upperCol;
		for (int i = 0; i < uSize; i++) {
			for (int j = 0; j < vSize; j++) {
				std::vector<Vec3> col;
				col.resize(wSize);
				for (int k = 0; k < wSize; k++)
					col[k] = cpts[i][j][k];
				subdivideCpts(col, w, lowerCol, upperCol);

				for (int k = 0; k < wSize; k++) {
					lowerCpts[i][j][k] = lowerCol[k];
					upperCpts[i][j][k] = upperCol[k];
				}
			}
		}

		lower = BezierVolume3d::create(uDegree, vDegree, wDegree, lowerCpts, buildMat);
		upper = BezierVolume3d::create(uDegree, vDegree, wDegree, upperCpts, buildMat);
	}
	BezierVolume3d::Ptr BezierVolume3d::subdivide(const Domain& uSubdomain, const Domain& vSubdomain, const Domain& wSubdomain) const {
		// Divide in u direction
		BezierVolume3d lower, upper;
		uSubdivide(uSubdomain.beg(), lower, upper, false);		// upper
		Real t = (uSubdomain.end() - uSubdomain.beg()) / (1.0 - uSubdomain.beg());
		upper.uSubdivide(t, lower, upper, false);					// lower

		// Divide in v direction
		lower.vSubdivide(vSubdomain.beg(), lower, upper, false);	// upper
		t = (vSubdomain.end() - vSubdomain.beg()) / (1.0 - vSubdomain.beg());
		upper.vSubdivide(t, lower, upper, false);					// lower

		// Divide in w direction
		lower.wSubdivide(wSubdomain.beg(), lower, upper, false);	// upper
		t = (wSubdomain.end() - wSubdomain.beg()) / (1.0 - wSubdomain.beg());
		upper.wSubdivide(t, lower, upper, false);					// lower

		lower.updateDerivMat();
		return std::make_shared<BezierVolume3d>(lower);
	}
	BezierVolume3d BezierVolume3d::create() {
		return empty;
	}
	BezierVolume3d BezierVolume3d::create(int uDegree, int vDegree, int wDegree, const ControlPoints& cpts, bool buildMat) {
		BezierVolume3d volume;
		volume.uDegree = uDegree;
		volume.vDegree = vDegree;
		volume.wDegree = wDegree;
		volume.uDomain.set(0, 1);
		volume.vDomain.set(0, 1);
		volume.wDomain.set(0, 1);
		volume.cpts = cpts;
		if (buildMat)
			volume.updateDerivMat();
		return volume;
	}
	BezierVolume3d::Ptr BezierVolume3d::createPtr(int uDegree, int vDegree, int wDegree, const ControlPoints& cpts, bool buildMat) {
		return std::make_shared<BezierVolume3d>(create(uDegree, vDegree, wDegree, cpts, buildMat));
	}

	void BezierVolume3d::updateDerivMat() {
		int uNum = (int)cpts.size();
		int vNum = (int)cpts[0].size();
		int wNum = (int)cpts[0][0].size();
		// First
		if (uDegree > 0) {
			// derivMatU
			derivMatU.resize(uNum - 1);
			for (int i = 0; i < uNum - 1; i++) {
				derivMatU[i].resize(vNum);
				for (int j = 0; j < vNum; j++) {
					derivMatU[i][j].resize(wNum);
					for (int k = 0; k < wNum; k++) {
						derivMatU[i][j][k] = (cpts[i + 1][j][k] - cpts[i][j][k]) * uDegree;
					}
				}
			}
		}
		if (vDegree > 0) {
			// derivMatV
			derivMatV.resize(uNum);
			for (int i = 0; i < uNum; i++) {
				derivMatV[i].resize(vNum - 1);
				for (int j = 0; j < vNum - 1; j++) {
					derivMatV[i][j].resize(wNum);
					for (int k = 0; k < wNum; k++) {
						derivMatV[i][j][k] = (cpts[i][j + 1][k] - cpts[i][j][k]) * vDegree;
					}
				}
			}
		}
		if (wDegree > 0) {
			// derivMatW
			derivMatW.resize(uNum);
			for (int i = 0; i < uNum; i++) {
				derivMatW[i].resize(vNum);
				for (int j = 0; j < vNum; j++) {
					derivMatW[i][j].resize(wNum - 1);
					for (int k = 0; k < wNum - 1; k++) {
						derivMatW[i][j][k] = (cpts[i][j][k + 1] - cpts[i][j][k]) * wDegree;
					}
				}
			}
		}

		// Second
		if (uDegree > 1) {
			// derivMatUU
			derivMatUU.resize(uNum - 2);
			for (int i = 0; i < uNum - 2; i++) {
				derivMatUU[i].resize(vNum);
				for (int j = 0; j < vNum; j++) {
					derivMatUU[i][j].resize(wNum);
					for (int k = 0; k < wNum; k++) {
						derivMatUU[i][j][k] = (derivMatU[i + 1][j][k] - derivMatU[i][j][k]) * (uDegree - 1);
					}
				}
			}
		}
		if (uDegree > 0 && vDegree > 0) {
			// derivMatUV
			derivMatUV.resize(uNum - 1);
			for (int i = 0; i < uNum - 1; i++) {
				derivMatUV[i].resize(vNum - 1);
				for (int j = 0; j < vNum - 1; j++) {
					derivMatUV[i][j].resize(wNum);
					for (int k = 0; k < wNum; k++) {
						derivMatUV[i][j][k] = (derivMatU[i][j + 1][k] - derivMatU[i][j][k]) * (vDegree);
					}
				}
			}
		}
		if (uDegree > 0 && wDegree > 0) {
			// derivMatUW
			derivMatUW.resize(uNum - 1);
			for (int i = 0; i < uNum - 1; i++) {
				derivMatUW[i].resize(vNum);
				for (int j = 0; j < vNum; j++) {
					derivMatUW[i][j].resize(wNum - 1);
					for (int k = 0; k < wNum - 1; k++) {
						derivMatUW[i][j][k] = (derivMatU[i][j][k + 1] - derivMatU[i][j][k]) * (wDegree);
					}
				}
			}
		}
		if (vDegree > 1) {
			// derivMatVV
			derivMatVV.resize(uNum);
			for (int i = 0; i < uNum; i++) {
				derivMatVV[i].resize(vNum - 2);
				for (int j = 0; j < vNum - 2; j++) {
					derivMatVV[i][j].resize(wNum);
					for (int k = 0; k < wNum; k++) {
						derivMatVV[i][j][k] = (derivMatV[i][j + 1][k] - derivMatV[i][j][k]) * (vDegree - 1);
					}
				}
			}
		}
		if (vDegree > 0 && wDegree > 0) {
			// derivMatVW
			derivMatVW.resize(uNum);
			for (int i = 0; i < uNum; i++) {
				derivMatVW[i].resize(vNum - 1);
				for (int j = 0; j < vNum - 1; j++) {
					derivMatVW[i][j].resize(wNum - 1);
					for (int k = 0; k < wNum - 1; k++) {
						derivMatVW[i][j][k] = (derivMatV[i][j][k + 1] - derivMatV[i][j][k]) * (wDegree);
					}
				}
			}
		}
		if (wDegree > 1) {
			// derivMatWW
			derivMatWW.resize(uNum);
			for (int i = 0; i < uNum; i++) {
				derivMatWW[i].resize(vNum);
				for (int j = 0; j < vNum; j++) {
					derivMatWW[i][j].resize(wNum - 2);
					for (int k = 0; k < wNum - 2; k++) {
						derivMatWW[i][j][k] = (derivMatW[i][j][k + 1] - derivMatW[i][j][k]) * (wDegree - 1);
					}
				}
			}
		}

		// Third
		if (uDegree > 2) {
			// derivMatUUU
			derivMatUUU.resize(uNum - 3);
			for (int i = 0; i < uNum - 3; i++) {
				derivMatUUU[i].resize(vNum);
				for (int j = 0; j < vNum; j++) {
					derivMatUUU[i][j].resize(wNum);
					for (int k = 0; k < wNum; k++) {
						derivMatUUU[i][j][k] = (derivMatUU[i + 1][j][k] - derivMatUU[i][j][k]) * (uDegree - 2);
					}
				}
			}
		}
		if (uDegree > 1 && vDegree > 0) {
			// derivMatUUV
			derivMatUUV.resize(uNum - 2);
			for (int i = 0; i < uNum - 2; i++) {
				derivMatUUV[i].resize(vNum - 1);
				for (int j = 0; j < vNum - 1; j++) {
					derivMatUUV[i][j].resize(wNum);
					for (int k = 0; k < wNum; k++) {
						derivMatUUV[i][j][k] = (derivMatUU[i][j + 1][k] - derivMatUU[i][j][k]) * (vDegree);
					}
				}
			}
		}
		if (uDegree > 1 && wDegree > 0) {
			// derivMatUUW
			derivMatUUW.resize(uNum - 2);
			for (int i = 0; i < uNum - 2; i++) {
				derivMatUUW[i].resize(vNum);
				for (int j = 0; j < vNum; j++) {
					derivMatUUW[i][j].resize(wNum - 1);
					for (int k = 0; k < wNum - 1; k++) {
						derivMatUUW[i][j][k] = (derivMatUU[i][j][k + 1] - derivMatUU[i][j][k]) * (wDegree);
					}
				}
			}
		}
		if (uDegree > 0 && vDegree > 1) {
			// derivMatUVV
			derivMatUVV.resize(uNum - 1);
			for (int i = 0; i < uNum - 1; i++) {
				derivMatUVV[i].resize(vNum - 2);
				for (int j = 0; j < vNum - 2; j++) {
					derivMatUVV[i][j].resize(wNum);
					for (int k = 0; k < wNum; k++) {
						derivMatUVV[i][j][k] = (derivMatUV[i][j + 1][k] - derivMatUV[i][j][k]) * (vDegree - 1);
					}
				}
			}
		}
		if (uDegree > 0 && vDegree > 0 && wDegree > 0) {
			// derivMatUVW
			derivMatUVW.resize(uNum - 1);
			for (int i = 0; i < uNum - 1; i++) {
				derivMatUVW[i].resize(vNum - 1);
				for (int j = 0; j < vNum - 1; j++) {
					derivMatUVW[i][j].resize(wNum - 1);
					for (int k = 0; k < wNum - 1; k++) {
						derivMatUVW[i][j][k] = (derivMatUV[i][j][k + 1] - derivMatUV[i][j][k]) * (wDegree);
					}
				}
			}
		}
		if (uDegree > 0 && wDegree > 1) {
			// derivMatUWW
			derivMatUWW.resize(uNum - 1);
			for (int i = 0; i < uNum - 1; i++) {
				derivMatUWW[i].resize(vNum);
				for (int j = 0; j < vNum; j++) {
					derivMatUWW[i][j].resize(wNum - 2);
					for (int k = 0; k < wNum - 2; k++) {
						derivMatUWW[i][j][k] = (derivMatUW[i][j][k + 1] - derivMatUW[i][j][k + 1]) * (wDegree - 1);
					}
				}
			}
		}
		if (vDegree > 2) {
			// derivMatVVV
			derivMatVVV.resize(uNum);
			for (int i = 0; i < uNum; i++) {
				derivMatVVV[i].resize(vNum - 3);
				for (int j = 0; j < vNum - 3; j++) {
					derivMatVVV[i][j].resize(wNum);
					for (int k = 0; k < wNum; k++) {
						derivMatVVV[i][j][k] = (derivMatVV[i][j + 1][k] - derivMatVV[i][j][k]) * (vDegree - 2);
					}
				}
			}
		}
		if (vDegree > 1 && wDegree > 0) {
			// derivMatVVW
			derivMatVVW.resize(uNum);
			for (int i = 0; i < uNum; i++) {
				derivMatVVW[i].resize(vNum - 2);
				for (int j = 0; j < vNum - 2; j++) {
					derivMatVVW[i][j].resize(wNum - 1);
					for (int k = 0; k < wNum - 1; k++) {
						derivMatVVW[i][j][k] = (derivMatVV[i][j][k + 1] - derivMatVV[i][j][k]) * (wDegree);
					}
				}
			}
		}
		if (vDegree > 0 && wDegree > 1) {
			// derivMatVWW
			derivMatVWW.resize(uNum);
			for (int i = 0; i < uNum; i++) {
				derivMatVWW[i].resize(vNum - 1);
				for (int j = 0; j < vNum - 1; j++) {
					derivMatVWW[i][j].resize(wNum - 2);
					for (int k = 0; k < wNum - 2; k++) {
						derivMatVWW[i][j][k] = (derivMatVW[i][j][k + 1] - derivMatVW[i][j][k]) * (wDegree - 1);
					}
				}
			}
		}
		if (wDegree > 2) {
			// derivMatWWW
			derivMatWWW.resize(uNum);
			for (int i = 0; i < uNum; i++) {
				derivMatWWW[i].resize(vNum);
				for (int j = 0; j < vNum; j++) {
					derivMatWWW[i][j].resize(wNum - 3);
					for (int k = 0; k < wNum - 3; k++) {
						derivMatWWW[i][j][k] = (derivMatWW[i][j][k + 1] - derivMatWW[i][j][k]) * (wDegree - 2);
					}
				}
			}
		}
	}
	Vec3 BezierVolume3d::evaluate(Real u, Real v, Real w) const {
		BasisVector uBasis, vBasis, wBasis;
		Bezier::calBasisVector(u, uDegree, uBasis);
		Bezier::calBasisVector(v, vDegree, vBasis);
		Bezier::calBasisVector(w, wDegree, wBasis);
		return tensorProduct(uBasis, vBasis, wBasis, cpts);
	}
	Vec3 BezierVolume3d::differentiate(Real u, Real v, Real w, int uOrder, int vOrder, int wOrder) const {
		BasisVector uBasis, vBasis, wBasis;
		if (uOrder == 0 && vOrder == 0 && wOrder == 0)
			return evaluate(u, v, w);
		else if (uOrder == 1 && vOrder == 0 && wOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			Bezier::calBasisVector(w, wDegree, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatU);
		}
		else if (uOrder == 0 && vOrder == 1 && wOrder == 0) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			Bezier::calBasisVector(w, wDegree, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatV);
		}
		else if (uOrder == 0 && vOrder == 0 && wOrder == 1) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			Bezier::calBasisVector(w, wDegree - 1, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatW);
		}
		else if (uOrder == 2 && vOrder == 0 && wOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 2, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			Bezier::calBasisVector(w, wDegree, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatUU);
		}
		else if (uOrder == 1 && vOrder == 1 && wOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			Bezier::calBasisVector(w, wDegree, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatUV);
		}
		else if (uOrder == 1 && vOrder == 0 && wOrder == 1) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			Bezier::calBasisVector(w, wDegree - 1, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatUW);
		}
		else if (uOrder == 0 && vOrder == 2 && wOrder == 0) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 2, vBasis);
			Bezier::calBasisVector(w, wDegree, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatVV);
		}
		else if (uOrder == 0 && vOrder == 1 && wOrder == 1) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			Bezier::calBasisVector(w, wDegree - 1, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatVW);
		}
		else if (uOrder == 0 && vOrder == 0 && wOrder == 2) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			Bezier::calBasisVector(w, wDegree - 2, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatWW);
		}
		else if (uOrder == 3 && vOrder == 0 && wOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 3, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			Bezier::calBasisVector(w, wDegree, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatUUU);
		}
		else if (uOrder == 2 && vOrder == 1 && wOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 2, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			Bezier::calBasisVector(w, wDegree, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatUUV);
		}
		else if (uOrder == 2 && vOrder == 0 && wOrder == 1) {
			Bezier::calBasisVector(u, uDegree - 2, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			Bezier::calBasisVector(w, wDegree - 1, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatUUW);
		}
		else if (uOrder == 1 && vOrder == 2 && wOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree - 2, vBasis);
			Bezier::calBasisVector(w, wDegree, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatUVV);
		}
		else if (uOrder == 1 && vOrder == 1 && wOrder == 1) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			Bezier::calBasisVector(w, wDegree - 1, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatUVW);
		}
		else if (uOrder == 1 && vOrder == 0 && wOrder == 2) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			Bezier::calBasisVector(w, wDegree - 2, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatUWW);
		}
		else if (uOrder == 0 && vOrder == 3 && wOrder == 0) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 3, vBasis);
			Bezier::calBasisVector(w, wDegree, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatVVV);
		}
		else if (uOrder == 0 && vOrder == 2 && wOrder == 1) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 2, vBasis);
			Bezier::calBasisVector(w, wDegree - 1, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatVVW);
		}
		else if (uOrder == 0 && vOrder == 1 && wOrder == 2) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			Bezier::calBasisVector(w, wDegree - 2, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatVWW);
		}
		else if (uOrder == 0 && vOrder == 0 && wOrder == 3) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			Bezier::calBasisVector(w, wDegree - 3, wBasis);
			return tensorProduct(uBasis, vBasis, wBasis, derivMatWWW);
		}
		else
			throw(std::runtime_error("Bezier volume differentiation is only allowed up to 3rd derivatives"));
	}
}