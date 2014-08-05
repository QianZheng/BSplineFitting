#include "read_write_asc.h"

#include <fstream>
#include <iostream>

bool CReadWriteAsc::readAsc( const string& filename, vector<Vector2d>& points )
{

	std::ifstream fin(filename);
	std::vector<double> coefficients;

	std::string line;
	int nRow = 0;
	int nCol = 0;
	if( fin.fail() )
		return false;


	while( std::getline( fin, line ) )
	{
		std::stringstream ss(line);
		double d = 0;
		nCol = 0;
		while( ss >> d)
		{
			coefficients.push_back( d);
			++nCol;
		}
		++nRow;
	}
	fin.close();

	if( nCol < 2)
		return false;



	points.resize( nRow);
	int idx = 0;
	for( int i = 0; i!= nRow; ++i )
	{
		points[i] = Vector2d( coefficients[idx+0], coefficients[idx+1]);
		idx += nCol;
	}
	return true;
}



bool CReadWriteAsc::writeAsc( const string& filename, const vector<Vector2d>& points )
{
	ofstream fout( filename.c_str());
	if( fout.fail() )
		return false;

	for( int i = 0; i!= points.size(); ++i )
	{
		fout << points[i].x() << " " << points[i].y()  << " 0" << endl;
	}
	fout.close();

	return true;

}
