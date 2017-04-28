#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double deltaL_prime;                                                               //
double K_L;                                                                        // Variables for L group
double S_L;                                                                        //
double L_macron_prime;                                                             //

double cal_deltaL_prime(double L_1, double L_2);                                   //
double cal_S_L(double L_macron_prime);                                             // Calculation functions for L group
double cal_L_macron_prime(double L_1, double L_2);                                 //

double deltaC_prime;                                                               //
double K_C;                                                                        //
double S_C;                                                                        //
double C_1_prime;                                                                  //
double C_2_prime;                                                                  // Variables for C group
double C_macron_prime;                                                             //
double C_macron;                                                                   //
double C_1;                                                                        //
double C_2;                                                                        //

double cal_deltaC_prime(double C_1_prime, double C_2_prime);                       //
double cal_C_prime(double a_prime, double b);                                      //
double cal_S_C(double C_macron_prime);                                             // Calculation functions for C group
double cal_C_macron_prime(double C_1_prime, double C_2_prime);                     //
double cal_C_macron(double C_1, double C_2);                                       //
double cal_C(double a, double b);                                                  //

double deltaH_prime;                                                               //
double K_H;                                                                        //
double S_H;                                                                        //
double delta_h_prime;                                                              // Variables for H group
double h_1_prime;                                                                  //
double h_2_prime;                                                                  //
double H_macron_prime;                                                             //

double cal_deltaH_prime(double C_1_prime, double C_2_prime, double delta_h_prime); //
double cal_S_H(double C_macron_prime, double T);                                   //
double cal_delta_h_prime(double h_1_prime, double h_2_prime);                      // Calculation functions for H group
double cal_h_prime(double a_prime, double b);                                      //
double cal_H_macron_prime(double h_1_prime, double h_2_prime);                     //

double a_1_prime;                                                                  //
double a_2_prime;                                                                  //
double G;                                                                          //
double T;                                                                          // Misc variables
double R_T;                                                                        //
double R_C;                                                                        //
double delta_theta;                                                                //

double cal_a_prime(double a, double G);                                            //
double cal_G(double C_macron);                                                     //
double cal_T(double H_macron_prime);                                               // Misc calculation functions
double cal_R_T(double R_C, double delta_theta);                                    //
double cal_R_C(double C_macron_prime);                                             //
double cal_delta_theta(double H_macron_prime);

double deltaE_L_group;                                                             //
double deltaE_C_group;                                                             // Helper variables to buld Delta E
double deltaE_H_group;                                                             //

double cal_deltaE_L_group(double deltaL_prime, double K_L, double S_L);            //
double cal_deltaE_C_group(double deltaC_prime, double K_C, double S_C);            // Helper groups to build Delta E
double cal_deltaE_H_group(double deltaH_prime, double K_H, double S_H);            //

double deltaE;                                                                     //
double cal_deltaE(double L_group, double C_group, double H_group, double R_T);     // Final build function

double run_deltaE(double L_1, double a_1, double b_1, double L_2, double a_2, double b_2);

int main (int argc, char* argv[]){
    if (argc != 7){
        printf("usage: %s L_1 a_1 b_1 L_2 a_2 b_2.\n", argv[0]);
        return 1;
    }
    double L_1, L_2, a_1, a_2, b_1, b_2;
    L_1 = atof(argv[1]);
    a_1 = atof(argv[2]);
    b_1 = atof(argv[3]);
    L_2 = atof(argv[4]);
    a_2 = atof(argv[5]);
    b_2 = atof(argv[6]);

    K_L = 1.0;
    K_C = 1.0;
    K_H = 1.0;

    run_deltaE( L_1, a_1, b_1, L_2, a_2, b_2 );

}

double cal_deltaL_prime(double L_1, double L_2){
    return L_2 - L_1;
}

double cal_L_macron_prime(double L_1, double L_2){
    return ( L_1 + L_2 ) / 2.0;
}

double cal_S_L(double L_macron_prime){
    return 1 + ( 0.015 * pow( L_macron_prime - 50, 2 ) ) \
            / sqrt( 20 + pow( L_macron_prime - 50, 2 ) );
}

double cal_C(double a, double b){
    return sqrt( a * a + b * b );
}

double cal_C_macron(double C_1, double C_2){
    return ( C_1 + C_2 ) / 2.0;
}

double cal_G(double C_macron){
    return 0.5 * ( 1 - sqrt( pow( C_macron, 7 ) / \
        ( pow( C_macron, 7 ) + pow( 25, 7 ) ) ) );
}

double cal_a_prime(double a, double G){
    return a * ( 1.0 + G );
}

double cal_C_prime(double a_prime, double b){
    return sqrt( a_prime * a_prime +  b * b );
}

double cal_deltaC_prime(double C_1_prime, double C_2_prime){
    return C_2_prime - C_1_prime;
}

double cal_C_macron_prime(double C_1_prime, double C_2_prime){
    return ( C_1_prime + C_2_prime ) / 2.0;
}

double cal_S_C(double C_macron_prime){
    return 1 + 0.045 * C_macron_prime;
}


double cal_h_prime(double a_prime, double b){
    double value;
    value = atan2( b, a_prime );
    if ( value >= 0 ){
        return value;
    } else {
        return value + 2.0 * M_PI;
    }
}

double cal_delta_h_prime(double h_1_prime, double h_2_prime){
    double value;
    value = h_2_prime - h_1_prime;
    if ( fabs( value ) <= M_PI ) {
        return value;
    } else if ( ( fabs ( value ) > M_PI ) && ( h_2_prime <= h_1_prime ) ) {
        return value + 2.0 * M_PI;
    } else {
        return value - 2.0 * M_PI;
    }
}

double cal_H_macron_prime(double h_1_prime, double h_2_prime){
    double value;
    value =  h_1_prime + h_2_prime;
    if( fabs ( h_1_prime - h_2_prime ) <= M_PI ){
        return value / 2.0;
    } else if ( value < 2 * M_PI ){
        return ( value / 2 ) + M_PI ;
    } else {
        return ( value / 2 ) - M_PI ;
    }
}

double cal_deltaH_prime(double C_1_prime, double C_2_prime, double delta_h_prime){
    return 2.0 * sqrt( C_1_prime * C_2_prime ) * sin( delta_h_prime / 2.0 );
}

double cal_T(double H_macron_prime){
   return 1 - 0.17 * cos( 1 * H_macron_prime - (     M_PI / 6  ) )  \
            + 0.24 * cos( 2 * H_macron_prime                     )  \
            + 0.32 * cos( 3 * H_macron_prime + (     M_PI / 30 ) )  \
            - 0.20 * cos( 4 * H_macron_prime - ( 7 * M_PI / 20 ) );
}

double cal_S_H(double C_macron_prime, double T){
    return 1 + 0.015 * C_macron_prime * T;
}

double cal_R_C(double C_macron_prime){
    return 2.0 * sqrt( pow( C_macron_prime, 7 ) / \
        ( pow( C_macron_prime, 7 ) + pow( 25, 7 ) ) );
}

double cal_delta_theta(double H_macron_prime){
    return ( M_PI / 6 ) * exp( -1 * pow( ( H_macron_prime - ( 55 * M_PI / 36 ) ) / 25, 2.0 ) );
}

double cal_R_T(double R_C, double delta_theta){
    return -1 * R_C * sin( 2.0 * delta_theta );
}

double cal_deltaE_L_group(double deltaL_prime, double K_L, double S_L){
    return deltaL_prime / ( K_L * S_L );
}
    
double cal_deltaE_C_group(double deltaC_prime, double K_C, double S_C){
    return deltaC_prime / ( K_C * S_C );
}

double cal_deltaE_H_group(double deltaH_prime, double K_H, double S_H){
    return deltaH_prime / ( K_H * S_H );
}

double cal_deltaE(double L_group, double C_group, double H_group, double R_T){
    return sqrt( pow( L_group, 2.0 ) + pow( C_group, 2.0 ) + pow( H_group, 2.0 ) + R_T * C_group * H_group );
}

double run_deltaE(double L_1, double a_1, double b_1, double L_2, double a_2, double b_2){
    deltaL_prime = cal_deltaL_prime( L_1, L_2 );
    L_macron_prime = cal_L_macron_prime( L_1, L_2);
    S_L = cal_S_L( L_macron_prime);

    C_1 = cal_C( a_1, b_1 );
    C_2 = cal_C( a_2, b_2 );
    C_macron = cal_C_macron( C_1, C_2 );
    G = cal_G( C_macron );
    a_1_prime = cal_a_prime( a_1, G );
    a_2_prime = cal_a_prime( a_2, G );
    C_1_prime = cal_C_prime( a_1_prime, b_1 );
    C_2_prime = cal_C_prime( a_2_prime, b_2 );
    deltaC_prime = cal_deltaC_prime( C_1_prime, C_2_prime );
    C_macron_prime = cal_C_macron_prime( C_1_prime, C_2_prime );
    S_C = cal_S_C( C_macron_prime );

    h_1_prime = cal_h_prime( a_1_prime, b_1 );
    h_2_prime = cal_h_prime( a_2_prime, b_2 );
    delta_h_prime = cal_delta_h_prime( h_1_prime, h_2_prime );
    H_macron_prime = cal_H_macron_prime( h_1_prime, h_2_prime );
    deltaH_prime = cal_deltaH_prime( C_1_prime, C_2_prime, delta_h_prime );
    T = cal_T( H_macron_prime );
    S_H = cal_S_H( C_macron_prime, T );
    R_C = cal_R_C( C_macron_prime );
    delta_theta = cal_delta_theta( H_macron_prime );
    R_T = cal_R_T( R_C, delta_theta );

    deltaE_L_group = cal_deltaE_L_group( deltaL_prime, K_L, S_L );
    deltaE_C_group = cal_deltaE_C_group( deltaC_prime, K_C, S_C );
    deltaE_H_group = cal_deltaE_H_group( deltaH_prime, K_H, S_H );

    deltaE = cal_deltaE( deltaE_L_group, deltaE_C_group, deltaE_H_group, R_T);
    
    printf("The delta E is: %f \n", deltaE);
    
    return deltaE;
}
