#include <iostream>
#include <cmath>

using std::cout;
using std::cin;

double powerf(double, int); // declaration

double powerf(double base, int exponent) // declaraing and defining
{
	double result =1;
	for(int i =0;i<exponent; i++)
	{
		result=result*base;
	}
	return result;
}

void print_pow(double base, int exponent)
{
	double myPower = powerf(base,exponent);
	cout<<base<<" rasied to the "<<exponent<<" power is "<<myPower<<std::endl;

}

int main() // main function
{ 
/*	std::cout << "Hello World\n"; */
	//int slices=5; // Datatype, indetifier, opertaor, value, semicolon | declaration- int slices; initialization- slices =5
	//cin >> slices;
	/* cout<<"Hello\n";*/
	//cout<<slices<<std::endl;
	//printf("%i\n",slices);
	//double power = pow(10,2);
	//cout<<power;
	/////
	double base; int exponent;
	cin >> base;
	cin >> exponent;
	//double myPower = powerf(base,exponent);
	//cout<<myPower<<std::endl;
	print_pow(base,exponent);
}





