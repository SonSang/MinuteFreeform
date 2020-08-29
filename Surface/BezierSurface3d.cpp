/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#include "BezierSurface3d.h"

namespace MN {
	void BezierSurface3d::subdivideCpts(const std::vector<Vec3>& cpts, Real t, std::vector<Vec3>& lower, std::vector<Vec3>& upper) {
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
	void BezierSurface3d::uSubdivide(Real u, BezierSurface3d& lower, BezierSurface3d& upper) const {
		ControlPoints lowerCpts = cpts;
		ControlPoints upperCpts = cpts;
		int rowNum = (int)cpts.size();
		int colNum = (int)cpts[0].size();
		std::vector<Vec3> lowerCol;
		std::vector<Vec3> upperCol;
		for (int i = 0; i < colNum; i++) {
			std::vector<Vec3> col;
			col.resize(rowNum);
			for (int j = 0; j < rowNum; j++)
				col[j] = cpts[j][i];
			subdivideCpts(col, u, lowerCol, upperCol);

			for (int j = 0; j < rowNum; j++) {
				lowerCpts[j][i] = lowerCol[j];
				upperCpts[j][i] = upperCol[j];
			}
		}
		lower = BezierSurface3d::create(uDegree, vDegree, lowerCpts, false);
		upper = BezierSurface3d::create(uDegree, vDegree, upperCpts, false);
	}
	void BezierSurface3d::vSubdivide(Real v, BezierSurface3d& lower, BezierSurface3d& upper) const {
		ControlPoints lowerCpts = cpts;
		ControlPoints upperCpts = cpts;
		int rowNum = (int)cpts.size();
		int colNum = (int)cpts[0].size();
		std::vector<Vec3> lowerRow;
		std::vector<Vec3> upperRow;
		for (int i = 0; i < rowNum; i++) {
			std::vector<Vec3> row;
			row.resize(colNum);
			for (int j = 0; j < colNum; j++)
				row[j] = cpts[i][j];
			subdivideCpts(row, v, lowerRow, upperRow);

			for (int j = 0; j < colNum; j++) {
				lowerCpts[i][j] = lowerRow[j];
				upperCpts[i][j] = upperRow[j];
			}
		}
		lower = BezierSurface3d::create(uDegree, vDegree, lowerCpts, false);
		upper = BezierSurface3d::create(uDegree, vDegree, upperCpts, false);
	}
	BezierSurface3d::Ptr BezierSurface3d::subdivide(const Domain& uSubdomain, const Domain& vSubdomain) const {
		// Divide in u direction
		BezierSurface3d lower, upper;
		uSubdivide(uSubdomain.beg(), lower, upper);	// upper
		Real t = (uSubdomain.end() - uSubdomain.beg()) / (1.0 - uSubdomain.beg());
		upper.uSubdivide(t, lower, upper);			// lower

		// Divide in v direction
		lower.vSubdivide(vSubdomain.beg(), lower, upper);	// upper
		t = (vSubdomain.end() - vSubdomain.beg()) / (1.0 - vSubdomain.beg());
		upper.vSubdivide(t, lower, upper);			// lower

		lower.updateDerivMat();
		return std::make_shared<BezierSurface3d>(lower);
	}
	BezierSurface3d BezierSurface3d::create(int uDegree, int vDegree, const ControlPoints& cpts, bool buildMat) {
		BezierSurface3d surface;
		surface.setDomain(0, Domain::create(0, 1));
		surface.setDomain(1, Domain::create(0, 1));
		surface.setDegree(0, uDegree);
		surface.setDegree(1, vDegree);
		surface.setCpts(cpts);
		if (buildMat)
			surface.updateDerivMat();
		return surface;
	}
	BezierSurface3d::Ptr BezierSurface3d::createPtr(int uDegree, int vDegree, const ControlPoints& cpts, bool buildMat) {
		BezierSurface3d surface = create(uDegree, vDegree, cpts, buildMat);
		return std::make_shared<BezierSurface3d>(surface);
	}
	void BezierSurface3d::updateDerivMat() {
		int rowNum = (int)cpts.size();
		int colNum = (int)cpts[0].size();
		if (uDegree > 0) {
			// derivMatU
			derivMatU.resize(rowNum - 1);
			for (int i = 0; i < rowNum - 1; i++) {
				derivMatU[i].resize(colNum);
				for (int j = 0; j < colNum; j++) {
					derivMatU[i][j] = (cpts[i + 1][j] - cpts[i][j]) * uDegree;
				}
			}
		}
		if (vDegree > 0) {
			// derivMatV
			derivMatV.resize(rowNum);
			for (int i = 0; i < rowNum; i++) {
				derivMatV[i].resize(colNum - 1);
				for (int j = 0; j < colNum - 1; j++) {
					derivMatV[i][j] = (cpts[i][j + 1] - cpts[i][j]) * vDegree;
				}
			}
		}
		if (uDegree > 1) {
			// derivMatUU
			derivMatUU.resize(rowNum - 2);
			for (int i = 0; i < rowNum - 2; i++) {
				derivMatUU[i].resize(colNum);
				for (int j = 0; j < colNum; j++) {
					derivMatUU[i][j] = (derivMatU[i + 1][j] - derivMatU[i][j]) * (uDegree - 1);
				}
			}
		}
		if (uDegree > 0 && vDegree > 0) {
			// derivMatUV
			derivMatUV.resize(rowNum - 1);
			for (int i = 0; i < rowNum - 1; i++) {
				derivMatUV[i].resize(colNum - 1);
				for (int j = 0; j < colNum - 1; j++) {
					derivMatUV[i][j] = (derivMatV[i + 1][j] - derivMatV[i][j]) * uDegree;
				}
			}
		}
		if (vDegree > 1) {
			// derivMatVV
			derivMatVV.resize(rowNum);
			for (int i = 0; i < rowNum; i++) {
				derivMatVV[i].resize(colNum - 2);
				for (int j = 0; j < colNum - 2; j++) {
					derivMatVV[i][j] = (derivMatV[i][j + 1] - derivMatV[i][j]) * (vDegree - 1);
				}
			}
		}
		if (uDegree > 2) {
			// derivMatUUU
			derivMatUUU.resize(rowNum - 3);
			for (int i = 0; i < rowNum - 3; i++) {
				derivMatUUU[i].resize(colNum);
				for (int j = 0; j < colNum; j++) {
					derivMatUUU[i][j] = (derivMatUU[i + 1][j] - derivMatUU[i][j]) * (uDegree - 2);
				}
			}
		}
		if (uDegree > 1 && vDegree > 0) {
			// derivMatUUV
			derivMatUUV.resize(rowNum - 2);
			for (int i = 0; i < rowNum - 2; i++) {
				derivMatUUV[i].resize(colNum - 1);
				for (int j = 0; j < colNum - 1; j++) {
					derivMatUUV[i][j] = (derivMatUU[i][j + 1] - derivMatUU[i][j]) * vDegree;
				}
			}
		}
		if (uDegree > 0 && vDegree > 1) {
			// derivMatUVV
			derivMatUVV.resize(rowNum - 1);
			for (int i = 0; i < rowNum - 1; i++) {
				derivMatUVV[i].resize(colNum - 2);
				for (int j = 0; j < colNum - 2; j++) {
					derivMatUVV[i][j] = (derivMatVV[i + 1][j] - derivMatVV[i][j]) * uDegree;
				}
			}
		}
		if (vDegree > 2) {
			// derivMatVVV
			derivMatVVV.resize(rowNum);
			for (int i = 0; i < rowNum; i++) {
				derivMatVVV[i].resize(colNum - 3);
				for (int j = 0; j < colNum - 3; j++) {
					derivMatVVV[i][j] = (derivMatVV[i][j + 1] - derivMatVV[i][j]) * (vDegree - 2);
				}
			}
		}
	}
	Vec3 BezierSurface3d::evaluate(Real u, Real v) const {
		BasisVector uBasis, vBasis;
		Bezier::calBasisVector(u, uDegree, uBasis);
		Bezier::calBasisVector(v, vDegree, vBasis);
		return tensorProduct(uBasis, cpts, vBasis);
	}
	Vec3 BezierSurface3d::differentiate(Real u, Real v, int uOrder, int vOrder) const {
		BasisVector uBasis, vBasis;
		if (uOrder == 0 && vOrder == 0)
			return evaluate(u, v);
		else if (uOrder == 1 && vOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			return tensorProduct(uBasis, derivMatU, vBasis);
		}
		else if (uOrder == 0 && vOrder == 1) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			return tensorProduct(uBasis, derivMatV, vBasis);
		}
		else if (uOrder == 2 && vOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 2, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			return tensorProduct(uBasis, derivMatUU, vBasis);
		}
		else if (uOrder == 1 && vOrder == 1) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			return tensorProduct(uBasis, derivMatUV, vBasis);
		}
		else if (uOrder == 0 && vOrder == 2) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 2, vBasis);
			return tensorProduct(uBasis, derivMatVV, vBasis);
		}
		else if (uOrder == 3 && vOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 3, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			return tensorProduct(uBasis, derivMatUUU, vBasis);
		}
		else if (uOrder == 2 && vOrder == 1) {
			Bezier::calBasisVector(u, uDegree - 2, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			return tensorProduct(uBasis, derivMatUUV, vBasis);
		}
		else if (uOrder == 1 && vOrder == 2) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree - 2, vBasis);
			return tensorProduct(uBasis, derivMatUVV, vBasis);
		}
		else if (uOrder == 0 && vOrder == 3) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 3, vBasis);
			return tensorProduct(uBasis, derivMatVVV, vBasis);
		}
		else
			throw(std::runtime_error("Bezier surface differentiation is only allowed up to 2nd derivatives"));
	}
}