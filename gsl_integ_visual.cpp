#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

#include <gsl/gsl_integration.h>

double com_dis(double x, void* Om, void* Or, void* OL) {
    double m = *(double*)Om;
    double r = *(double*)Or;
    double L = *(double*)OL;
    double s = m * (pow(1 + x, 3)) + r * (pow(1 + x, 2)) + L;
    //second term is 0 because \Omega_r = 0;
    double com_distance = 3000 / (0.7 * pow(s, 0.5));
    // c/H_0 = 3000 MPc
    return com_distance;
}

int
main(void)
{
    gsl_integration_workspace* work_ptr
        = gsl_integration_workspace_alloc(1000);

    double lower_limit = 0;	/* lower limit a */
    double upper_limit = 1;	/* upper limit b */
    double abs_error = 1.0e-8;	/* to avoid round-off problems */
    double rel_error = 1.0e-8;	/* the result will usually be much better */
    double result;		/* the result from the integration */
    double error;			/* the estimated error from the integration */

    double m = 0.3;
    double r = 0;
    double L = 0.7	// parameter in integrand
        double expected = 2314.28;	// exact answer

    gsl_function My_function;
    void* Om_ptr = &m;
    void* Ob_ptr = &r;
    void* OL_ptr = &L;

    My_function.function = &com_dis;
    My_function.Om = Om_ptr;
    My_function.Ob = Ob_ptr;
    My_function.OL = OL_ptr;

    gsl_integration_qags(&My_function, lower_limit, upper_limit,
        abs_error, rel_error, 1000, work_ptr, &result,
        &error);

    cout.setf(ios::fixed, ios::floatfield);	// output in fixed format
    cout.precision(18);		// 18 digits in doubles

    int width = 20;  // setw width for output
    cout << "result          = " << setw(width) << result << endl;
    cout << "exact result    = " << setw(width) << expected << endl;
    cout << "estimated error = " << setw(width) << error << endl;
    cout << "actual error    = " << setw(width) << result - expected << endl;
    cout << "intervals =  " << work_ptr->size << endl;

    return 0;
}
