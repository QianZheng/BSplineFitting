#include "cubic_b_spline.h"
#include <fstream>

#include <ANN/ANN.h>




void CubicBSplineCurve::setNewControl( const vector<Vector2d>& controlPs)
{
	clear();
	controls_ = controlPs;

	for( unsigned int i = 0; i<nb_control(); i++)
	{
		for( double fj = 0; fj <=1.0f; fj+= interal_)
		{
			Parameter temp(i,fj); 
			Vector2d p = getPos( temp ) ;
			positions_.push_back(p);
		}
	}
}


//************************************
// Method:    getPos
// Returns:   Eigen::Vector2d
// Function:  公式B(t)的展开形式
// Time:      2014/08/05
// Author:    Qian
//************************************
Vector2d CubicBSplineCurve::getPos( const Parameter& para) const
{
	MatrixXd cm(4,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;

	double tf = para.second;
	int ki = para.first; 

	MatrixXd  tm(1,4);
	tm << tf*tf*tf, tf*tf, tf, 1;


	size_t n = nb_control();
	MatrixXd pm(4,2);
	for( int i = 0; i < 4; i++)
	{
		pm(i,0) = controls_[(ki+i)%n].x()/6.0;
		pm(i,1) = controls_[(ki+i)%n].y()/6.0;
	}
	MatrixXd rm = tm*cm*pm;

	return Vector2d( rm(0,0), rm(0,1));
}

Vector2d CubicBSplineCurve::getFirstDiff(const Parameter& para) const
{
	MatrixXd cm(4,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;

	double tf = para.second;
	int ki = para.first; 

	MatrixXd  tm(1,4);
	tm << 3*tf*tf,2*tf, 1, 0;

	size_t n = nb_control();
	MatrixXd pm(4,2);
	for( int i = 0; i < 4; i++)
	{
		pm(i,0) =  controls_[(ki+i)%n].x()/ 6.0f;
		pm(i,1) =  controls_[(ki+i)%n].y() / 6.0f;
	}


	MatrixXd rm = tm*cm*pm;

	return Vector2d( rm(0,0), rm(0,1));

}


Vector2d CubicBSplineCurve::getSecondDiff( const Parameter& para ) const
{
	MatrixXd cm(4,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;


	double tf = para.second;
	int ki = para.first; 
	MatrixXd  tm(1,4);
	tm << 6*tf,2, 0, 0;

	int n = nb_control();
	MatrixXd pm(4,2);
	for( int i = 0; i < 4; i++)
	{
		pm(i,0) = controls_[(ki+i)%n].x()/6.0;
		pm(i,1) =  controls_[(ki+i)%n].y()/6.0;
	}
	MatrixXd rm = tm*cm*pm;

	return Vector2d( rm(0,0), rm(0,1));

}


// Refer: http://en.wikipedia.org/wiki/Curvature
double  CubicBSplineCurve::getCurvature(const Parameter& para)  const
{
	Vector2d fp = getFirstDiff( para );
	Vector2d sp = getSecondDiff( para );

	double kappa = abs( fp.x()*sp.y() - sp.x()*fp.y() );
	kappa = kappa / sqrt( pow( ( fp.x()*fp.x()+fp.y()*fp.y()), 3) );

	return kappa;
}



Vector2d CubicBSplineCurve::getTangent( const Parameter &para ) const
{
	Vector2d p = getFirstDiff(para);
	return p.normalized();

}

Vector2d CubicBSplineCurve::getNormal( const Parameter &para ) const
{
	Vector2d v = getTangent( para );
	return Vector2d( -v.y(), v.x() );

}


Vector2d CubicBSplineCurve::getCurvCenter( const Parameter &para) const
{
	Vector2d p = getPos(para);

	Vector2d fd = getFirstDiff( para );
	Vector2d sd = getSecondDiff( para );

	double p1 = ( fd.x()*fd.x() + fd.y()*fd.y() ) * fd.y();
	double p2 =  sd.y()*fd.x()-sd.x()*fd.y() ;
	double alpha = p.x() - p1/p2;

	double p3 = ( fd.x()*fd.x() + fd.y()*fd.y() ) * fd.x();
	double beta = p.y() + p3/p2;


	return Vector2d(alpha ,beta);

}


double CubicBSplineCurve::findFootPrint(const vector<Vector2d>& givepoints, 
									   vector<Parameter>& footPrints) const
{
	footPrints.clear();
	footPrints.resize( givepoints.size(), Parameter(0,0.0) );

	int iKNei = 1;
	int iDim = 2;
	size_t iNPts = positions_.size();
	double eps = 0;

	ANNpointArray dataPts = annAllocPts(iNPts, iDim); // allocate data points; // data points
	ANNpoint queryPt = annAllocPt(iDim);  // allocate query point

	ANNidxArray nnIdx = new ANNidx[iKNei]; // allocate near neigh indices
	ANNdistArray dists = new ANNdist[iKNei]; // allocate near neighbor dists

	for( int i = 0; i!= iNPts; ++i) {
		dataPts[i][0] = positions_[i].x();
		dataPts[i][1] = positions_[i].y();
	}
	ANNkd_tree* kdTree = new ANNkd_tree( // build search structure
		dataPts, // the data points
		iNPts, // number of points
		iDim);

	double squareSum = 0.0;
	for( int i = 0 ;i!= (int)givepoints.size(); ++i) {	
		queryPt[0] = givepoints[i].x();
		queryPt[1] = givepoints[i].y();
		kdTree->annkSearch( // search
			queryPt, // query point
			iKNei, // number of near neighbors
			nnIdx, // nearest neighbors (returned)
			dists, // distance (returned)
			eps); // error bound
		squareSum += dists[0];
		footPrints[i] =  getPara(nnIdx[0]) ;
	}
	
	delete[] nnIdx;
	delete[] dists;
	delete kdTree;
	annClose(); // done with ANN

	return squareSum;
}


CubicBSplineCurve::Parameter CubicBSplineCurve::getPara( int index ) const
{
	int num = (int)( positions_.size()/controls_.size());
	int ki = index/num;
	double tf = interal_*( index - ki*num );
	return make_pair( ki, tf);
}



VectorXd CubicBSplineCurve::getCoffe( const Parameter& para) const
{
	int ki = para.first;
	double tf = para.second;

	Matrix4d cm(4,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;

	MatrixXd  tv(1,4);
	tv << tf*tf*tf, tf*tf, tf, 1;

	MatrixXd rv = tv*cm; 

	VectorXd newv(nb_control());
	newv.setZero();
	for( int i = 0; i < 4; i++)
	{
		newv[ (ki+i)%nb_control() ] = rv(0,i)/6.0f;
	}
	return newv;
}


//temporary solution 
bool CubicBSplineCurve::checkSameSide(Vector2d p1, Vector2d p2 , Vector2d neip)
{
	Vector2d v1 = p2 - neip;
	Vector2d v2 = p1 - neip;
	bool b = true;

	if( v1.x()*v2.x() + v1.y()*v2.y() < 0)
	{
		b = false;
	}

	return  b; 
}



bool CubicBSplineCurve::checkInside(Vector2d p)
{
	int strip = 0.02/interal_;
	int    wn = 0;    // the winding number counter
	// loop through all edges of the polygon
	for (int i=0; i< (int)positions_.size(); i+=strip) 
	{   
		int j = (i+strip)/(int)positions_.size();
		// edge from V[i] to V[j]
		if (positions_[i].y() <= p.y() ) {       
			// start y <= P.y
			if ( positions_[j].y() > p.y() )      // an upward crossing
				if (isLeft( positions_[i], positions_[j], p) > 0)  // P left of edge
					++wn;            // have a valid up intersect
		}
		else {                       // start y > P.y (no test needed)
			if ( positions_[j].y() <= p.y())     // a downward crossing
				if ( isLeft( positions_[i], positions_[j], p) <0 ) // P right of edge
					--wn;            // have a valid down intersect
		}
	}
	if(wn == 0) 
		return false;
	else   
		return true; 

}

int CubicBSplineCurve::isLeft( Vector2d p0, Vector2d p1, Vector2d p2)
{
	return ( (p1.x() - p0.x()) * (p2.y() - p0.y())
		- (p2.x() - p0.x()) * (p1.y() - p0.y()) );
}



MatrixXd CubicBSplineCurve::getSIntegralSq()
{
	// compute P"(t)
	int controlNum  = nb_control();
	MatrixXd pm(2*controlNum, 2*controlNum);
	pm.setZero();

	Matrix2d tIntergrated;
	tIntergrated << 1/3.0, 1/2.0, 1/2.0, 1.0;
	Matrix2d tm;
	tm << 6, 0, 0, 2;


	MatrixXd cm(2,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0;
	cm = cm/6.0;


	Matrix4d coffm = cm.transpose()*tm.transpose()*tIntergrated*tm*cm; 
	for( int i = 0; i < controlNum; i++)
	{
		for( int j = 0; j < 4; j++)
		{
			for( int n = 0; n < 4; n++)
			{
				int kj = ( i+j ) % controlNum;
				int kn = (i+n) % controlNum;
				pm(kj,kn) += coffm(j,n);
				pm(controlNum+kj,controlNum+kn) += coffm(j,n);
			}
		}
	}
	return pm;
}


MatrixXd CubicBSplineCurve::getFIntegralSq()
{
	// compute P"(t)
	size_t controlNum  = nb_control();
	MatrixXd pm(2*controlNum, 2*controlNum);
	pm.setZero();


	Matrix3d tIntergrated;
	tIntergrated << 1/5.0 , 1/4.0, 1/3.0,
		1/4.0, 1/3.0, 1/2.0,
		1/3.0, 1/2.0, 1/1.0;

	MatrixXd cm(3,4);
	cm << -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0;
	cm = cm/6.0;

	Matrix3d tm;
	tm << 3, 0, 0, 0,2, 0, 0, 0,1;

	Matrix4d coffm = cm.transpose()*tm.transpose()*tIntergrated*tm*cm; 
	for( int i = 0; i < controlNum; i++)
	{
		for( int j = 0; j < 4; j++)
		{
			for( int n = 0; n < 4; n++)
			{
				int kj = ( i+j ) % controlNum;
				int kn = (i+n) % controlNum;
				pm(kj,kn) += coffm(j,n);
				pm(controlNum+kj,controlNum+kn) += coffm(j,n);
			}
		}
	}
	return pm;

}
