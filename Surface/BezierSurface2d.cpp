#include "BezierSurface2d.h"

namespace MN {
	// BezierSurface2d
	static const BezierSurface2d empty2d = BezierSurface2d::create(0, 0, {}, false);
	void BezierSurface2d::subdivideCpts(const std::vector<Vec2>& cpts, Real t, std::vector<Vec2>& lower, std::vector<Vec2>& upper) {
		// De Casteljou's algorithm
		int size = (int)cpts.size();
		lower.resize(size);
		upper.resize(size);

		Real t1 = 1 - t;

		std::vector<Vec2> copy = cpts;
		lower[0] = copy.front();
		upper[size - 1] = copy.back();

		int cnt = 1;
		for (int i = size - 1; i > 0; i--) {
			std::vector<Vec2> tmp;
			tmp.resize(i);
			for (int j = 0; j < i; j++)
				tmp[j] = copy[j] * t1 + copy[j + 1] * t;
			copy = tmp;

			lower[cnt] = copy.front();
			upper[size - cnt - 1] = copy.back();
			cnt++;
		}
	}
	void BezierSurface2d::uSubdivide(Real u, BezierSurface2d& lower, BezierSurface2d& upper) const {
		ControlPoints lowerCpts = cpts;
		ControlPoints upperCpts = cpts;
		int rowNum = (int)cpts.size();
		int colNum = (int)cpts[0].size();
		std::vector<Vec2> lowerCol;
		std::vector<Vec2> upperCol;
		for (int i = 0; i < colNum; i++) {
			std::vector<Vec2> col;
			col.resize(rowNum);
			for (int j = 0; j < rowNum; j++)
				col[j] = cpts[j][i];
			subdivideCpts(col, u, lowerCol, upperCol);

			for (int j = 0; j < rowNum; j++) {
				lowerCpts[j][i] = lowerCol[j];
				upperCpts[j][i] = upperCol[j];
			}
		}
		lower = BezierSurface2d::create(uDegree, vDegree, lowerCpts, false);
		upper = BezierSurface2d::create(uDegree, vDegree, upperCpts, false);
	}
	void BezierSurface2d::vSubdivide(Real v, BezierSurface2d& lower, BezierSurface2d& upper) const {
		ControlPoints lowerCpts = cpts;
		ControlPoints upperCpts = cpts;
		int rowNum = (int)cpts.size();
		int colNum = (int)cpts[0].size();
		std::vector<Vec2> lowerRow;
		std::vector<Vec2> upperRow;
		for (int i = 0; i < rowNum; i++) {
			std::vector<Vec2> row;
			row.resize(colNum);
			for (int j = 0; j < colNum; j++)
				row[j] = cpts[i][j];
			subdivideCpts(row, v, lowerRow, upperRow);

			for (int j = 0; j < colNum; j++) {
				lowerCpts[i][j] = lowerRow[j];
				upperCpts[i][j] = upperRow[j];
			}
		}
		lower = BezierSurface2d::create(uDegree, vDegree, lowerCpts, false);
		upper = BezierSurface2d::create(uDegree, vDegree, upperCpts, false);
	}
	BezierSurface2d::Ptr BezierSurface2d::subdivide(const Domain& uSubdomain, const Domain& vSubdomain) const {
		// Divide in u direction
		BezierSurface2d lower, upper;
		uSubdivide(uSubdomain.beg(), lower, upper);	// upper
		Real t = (uSubdomain.end() - uSubdomain.beg()) / (1.0 - uSubdomain.beg());
		upper.uSubdivide(t, lower, upper);			// lower

		// Divide in v direction
		lower.vSubdivide(vSubdomain.beg(), lower, upper);	// upper
		t = (vSubdomain.end() - vSubdomain.beg()) / (1.0 - vSubdomain.beg());
		upper.vSubdivide(t, lower, upper);			// lower

		lower.updateDerivMat();
		return std::make_shared<BezierSurface2d>(lower);
	}
	BezierSurface2d BezierSurface2d::create() {
		return empty2d;
	}
	BezierSurface2d BezierSurface2d::create(int uDegree, int vDegree, const ControlPoints& cpts, bool buildMat) {
		BezierSurface2d surface;
		surface.setDomain(0, Domain::create(0, 1));
		surface.setDomain(1, Domain::create(0, 1));
		surface.setDegree(0, uDegree);
		surface.setDegree(1, vDegree);
		surface.setCpts(cpts);
		if (buildMat)
			surface.updateDerivMat();
		return surface;
	}
	BezierSurface2d::Ptr BezierSurface2d::createPtr(int uDegree, int vDegree, const ControlPoints& cpts, bool buildMat) {
		BezierSurface2d surface = create(uDegree, vDegree, cpts, buildMat);
		return std::make_shared<BezierSurface2d>(surface);
	}
	void BezierSurface2d::updateDerivMat() {
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
	Vec2 BezierSurface2d::evaluate(Real u, Real v) const {
		BasisVector uBasis, vBasis;
		Bezier::calBasisVector(u, uDegree, uBasis);
		Bezier::calBasisVector(v, vDegree, vBasis);
		return tensorProduct(uBasis, cpts, vBasis);
	}
	Vec2 BezierSurface2d::differentiate(Real u, Real v, int uOrder, int vOrder) const {
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
	void BezierSurface2d::uIsoCurve(Real u, BezierCurve2d& curve) const {
		BezierCurve2d::ControlPoints curveCpts;
		curveCpts.resize(cpts[0].size());
		for (auto& pt : curveCpts)
			pt = Vec2::zero();

		BasisVector uBasis;
		Bezier::calBasisVector(u, uDegree, uBasis);
		for (int i = 0; i < (int)cpts.size(); i++)
			for (int j = 0; j < (int)cpts[0].size(); j++)
				curveCpts[j] += cpts[i][j] * uBasis[i];

		curve.setCpts(curveCpts);
		curve.setDegree(vDegree);
		curve.setDomain(vDomain);
		curve.updateDerivMat();
	}
	void BezierSurface2d::vIsoCurve(Real v, BezierCurve2d& curve) const {
		BezierCurve2d::ControlPoints curveCpts;
		curveCpts.resize(cpts.size());
		for (auto& pt : curveCpts)
			pt = Vec2::zero();

		BasisVector vBasis;
		Bezier::calBasisVector(v, vDegree, vBasis);
		for (int i = 0; i < (int)cpts[0].size(); i++)
			for (int j = 0; j < (int)cpts.size(); j++)
				curveCpts[j] += cpts[j][i] * vBasis[i];

		curve.setCpts(curveCpts);
		curve.setDegree(uDegree);
		curve.setDomain(uDomain);
		curve.updateDerivMat();
	}
}