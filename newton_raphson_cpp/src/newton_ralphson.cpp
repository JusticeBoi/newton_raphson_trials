#include <eigen3/Eigen/Dense>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <chrono>
#include <limits>

//#define DEBUG
using Eigen::Matrix2d;
using Eigen::Vector2d;

int main(int argc, char* argv[])
{
	assert(argc == 4 && "need exactly 3 arguments");
	int max_iteration = std::atoi(argv[3]);
	double x1 = std::atof(argv[1]);
	double x2 = std::atof(argv[2]);
	std::cout << "your initial guess for x1: " <<x1 <<" x2: "<<x2<<" and max iteration is : "<<max_iteration<< std::endl;


	Vector2d F,dFdx2,dFdx1;
	Matrix2d M ;


	auto start_read = std::chrono::steady_clock::now();
	for(int number_of_iterations = 0 ; number_of_iterations < max_iteration ;number_of_iterations++)
	{
		F << sin(x1)-cos(x2),cos(x1)-sin(x2);
		dFdx1<<cos(x1),-sin(x1);

		dFdx2<<sin(x2),-cos(x2);
		M << dFdx1, dFdx2 ;
		if(F(0) < std::numeric_limits<double>::epsilon() && F(1) < std::numeric_limits<double>::epsilon() ) break;
		Vector2d sol = M.fullPivLu().solve(-F) ;
		x1+=sol(0);
		x2+=sol(1);
				#ifdef DEBUG
				std::cout <<"F: \n"<< F << std::endl;
				std::cout <<"dFdx1: \n"<< dFdx1 << std::endl;
				std::cout <<"dFdx2: \n"<< dFdx2 << std::endl;
				std::cout <<"solved M x = -F : \n" << M <<"\n\n"<< -F <<"\n"<<"with solution : \n"<< sol<<"\n"<<std::endl;
				std::cout <<"x1 : "<<x1 << " x2 : " << x2 << std::endl;
				#endif

	}
	auto end_read = std::chrono::steady_clock::now();
	auto diff_read = end_read - start_read;
	std::cout <<"duration of newton_ralphson "<< std::chrono::duration <double, std::milli> (diff_read).count() << " ms" << std::endl;



	std::cout <<std::setprecision(30)<< "x1 : " <<x1 << " x2 : " << x2 <<" with error: \n\n"<<F<< std::endl;

	return 0;
}

