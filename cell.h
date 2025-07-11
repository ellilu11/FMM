#include <vector>
#include <complex>

using namespace std;

class Cell
{
public:
	Cell(vector<complex<double>> &,
		 const double,
		 vector<Cell> &)

private:
	const std::complex center;
	const double s;
	vector<double> qs;
	vector<complex<double>> pos;
	vector<Cell> leaves;

};