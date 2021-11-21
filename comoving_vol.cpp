#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_integration.h>
using namespace std;

struct f_params {
    double Om;
    double Or;
    double OL;
};

struct g_params {
    double theta1;
    double theta2;
    double phi1;
    double phi2;
    double z1;
    double z2;
};

double com_dis(double x, void* par) {
    struct f_params* params = (struct f_params*)par;
    double m = params->Om;
    double r = params->Or;
    double L = params->OL;
    double s = m * (pow(1 + x, 3)) + r * (pow(1 + x, 4)) + L;
    //second term is 0 because \Omega_r = 0;
    double com_distance = 3000 / (0.7 * pow(s, 0.5));
    // c/H_0 = 3000 MPc
    return com_distance;
}

double Co_vol(double theta1, double theta2, double phi1, double phi2, double z1, double z2) {
    //struct g_params* params = (struct g_params*)p;
    double theta_v = cos(theta1) - cos(theta2);
    double phi_v = phi1 - phi2;

    //calculating r1 and r2 using GSL routine
    gsl_integration_workspace* work_ptr
        = gsl_integration_workspace_alloc(1000);

    double low = 0;	// lower limit a 
    double up1 = z1;
    double up2 = z2;    // upper limit a 
    double abs_error = 1.0e-8;	// to avoid round-off problems 
    double rel_error = 1.0e-8;	// the result will usually be much better 
    double result;		// the result from the integration //
    double error;

    double m = 0.3;
    double r = 0;
    double L = 0.7;

    gsl_function F;
    struct f_params  params = { m,r,L };

    F.function = &com_dis;
    F.params = &params;
    double r1 = gsl_integration_qags(&F, low, up1,
        abs_error, rel_error, 1000, work_ptr, &result,
        &error);
    double r2 = gsl_integration_qags(&F, low, up2,
        abs_error, rel_error, 1000, work_ptr, &result,
        &error);
    double r_v = (pow(r2,3) - pow(r1,3)) / 3;
    return theta_v * phi_v * r_v;
}

int main(void) {
    double RA1h;
    double RA1m;
    double RA1s;

    double RA2h;
    double RA2m;
    double RA2s;

    cout << "Enter the range of RA i.e. at first, RA(i)(Hours, minutes and seconds one by one)" << endl;
    cout << "Hours: " << endl;
    cin >> RA1h;
    cout << "Minutes: " << endl;
    cin >> RA1m;
    cout << "Seconds: " << endl;
    cin >> RA1s;

    cout << "Enter the range of RA i.e. RA(f)(Hours, minutes and seconds one by one)" << endl;
    cout << "Hours: " << endl;
    cin >> RA2h;
    cout << "Minutes: " << endl;
    cin >> RA2m;
    cout << "Seconds: " << endl;
    cin >> RA2s;
    double DEC1;
    double DEC2;
    cout << "Enter the range of DEC i.e. DEC(i) and DEC(f)" << endl;
    cin >> DEC1;
    cin >> DEC2;
    double Z1;
    double Z2;
    cout << "Enter the range of Z i.e. Z(i) and Z(f)" << endl;
    cin >> Z1;
    cin >> Z2;
    double theta1 = 90.0 - DEC1;
    double theta2 = 90.0 - DEC2;
    double phi1 = (RA1h + RA1m/60.0 + RA1s/3600.0)*15*0.0174533;
    double phi2 = (RA2h + RA2m / 60.0 + RA2s / 3600.0) * 15 * 0.0174533;
    cout.precision(9);		// 9 digits in doubles
    cout << "Here is the comoving volume for the given patch" << Co_vol(theta1, theta2, phi1, phi2, Z1, Z2) << endl;
    return 0;
}