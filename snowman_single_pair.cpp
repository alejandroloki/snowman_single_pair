/*****************************/
/*a basic 3D MPC solvent code*/
/*****************************/


#include <stdlib.h>
#include <math.h>
#include <time.h>
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <malloc.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include <limits.h>
#define nxcell 32
#define nycell 32
#define nzcell 20
#define a 40000
#define k1 50000
#define ncolloid 2 //total colloid number
#define nparmpc 20480 //total mpc particle number
#define vavempc sqrt(1.5) //the root mean square velocity, kT=0.5, mmpc=1
#define eqtime 10000    //the number of simulation step
#define MIN_R2 1e-6      //minimum distance squared to prevent force explosion
#define MAX_FORCE 1e6    //maximum force value to prevent numerical overflow
#define ENABLE_NUMERICAL_CHECKS 0  // Set to 1 to enable numerical checks (slower), 0 for release mode (faster)
using namespace std;

double phi, theta;
double square_rcball1, square_rcball2;
double xlen = 1.0 * nxcell;   // system size
double ylen = 1.0 * nycell;
double zlen = 1.0 * nzcell;
double xhalf = xlen / 2;
double yhalf = ylen / 2;
double Pi = acos(-1.0);
static double rxmpc[nparmpc + 1];  // particle position after collision
static double rympc[nparmpc + 1];
static double rzmpc[nparmpc + 1];
static double vxmpc[nparmpc + 1];  // particle velocity of collision
static double vympc[nparmpc + 1];
static double vzmpc[nparmpc + 1];
double ang = 120 * Pi / 180;   // rotational angle
double rho = 10;
double dis = 2.5;
double sinang = sin(ang);
double cosang = cos(ang);
double dtmpc = 0.1; // MPC time step
double dtmd = 0.002;  //MD time step
static double fxmpc[nparmpc + 1]; // interaction in x direction
static double fympc[nparmpc + 1]; // interaction in y direction
static double fzmpc[nparmpc + 1]; // interaction in z direction
double sigma[ncolloid + 1];
double rc[ncolloid + 1];
double rn[ncolloid + 1];
double square_sigma1, square_sigma2;
double rd1 = 2.05;
double square_rd1 = rd1 * rd1;
double rbottem1 = 2.1;
double square_rbottem1 = rbottem1 * rbottem1;
double rc_att1 = 2.2;
double square_rc_att1 = rc_att1 * rc_att1;
double rd2 = 3.12;
double square_rd2 = rd2 * rd2;
double rbottem2 = 3.18;
double square_rbottem2 = rbottem2 * rbottem2;
double rc_att2 = 3.3;
double square_rc_att2 = rc_att2 * rc_att2;
double rw = 10.0;
double square_rw = rw * rw;
double rwall[ncolloid + 1];
double rcwall[ncolloid + 1];
double square_rwall1, square_rwall2;
double rball[ncolloid + 1];
double rcball[ncolloid + 1];
double square_rball1, square_rball2;
double square_rc[ncolloid + 1];
double square_rca[ncolloid + 1];
double square_rn[ncolloid + 1];
double rcx0[ncolloid + 1], rcy0[ncolloid + 1], rcz0[ncolloid + 1];
double rcx[ncolloid + 1], rcy[ncolloid + 1], rcz[ncolloid + 1];
double fcx0[ncolloid + 1], fcy0[ncolloid + 1], fcz0[ncolloid + 1];
double fcx[ncolloid + 1], fcy[ncolloid + 1], fcz[ncolloid + 1];
double vcx[ncolloid + 1], vcy[ncolloid + 1], vcz[ncolloid + 1];
double omegax[ncolloid + 1], omegay[ncolloid + 1], omegaz[ncolloid + 1];
double Momx0[ncolloid + 1], Momy0[ncolloid + 1], Momz0[ncolloid + 1];
double Momx[ncolloid + 1], Momy[ncolloid + 1], Momz[ncolloid + 1];
double nx[ncolloid + 1][3], ny[ncolloid + 1][3], nz[ncolloid + 1][3];
double nx_tor[ncolloid + 1][3], ny_tor[ncolloid + 1][3], nz_tor[ncolloid + 1][3];
double torquex[ncolloid + 1][3], torquey[ncolloid + 1][3], torquez[ncolloid + 1][3], torque[ncolloid + 1][3];
double costheta[ncolloid + 1][3], theta_st[ncolloid + 1][3];
double square_rea;
double square_bb[ncolloid + 1];
double Mass[ncolloid + 1];
double I[ncolloid + 1];
double Miu[ncolloid + 1];

/*random number generator*/

#define AMOL .23283064365386962890625e-9
#define NBIT 32
static int ir[256];
static int k;
double r250()
{
    int iran_ks;
    k = (k + 1) & 255;
    ir[k] = ir[(k - 250 & 255)] ^ ir[(k - 103 & 255)];
    iran_ks = ir[k];
    return iran_ks * AMOL + 0.5;
}
void wmup_ks(int nseed)
{
    int ibm, idum, i;
    double rdum;
    ibm = 2 * nseed + 1;
    for (k = 0; k < 256; k++)
    {
        idum = 0;
        for (i = 1; i <= NBIT; i++)
        {
            idum = idum * 2;
            ibm = ibm * 16807;
            if (ibm < 0) idum++;
        }
        ir[k] = idum;
    }
    k = 0;
    for (i = 0; i < 100000; i++)
        rdum = r250();
}

// Check if a double value is valid (not NaN, not Inf, within reasonable range)
int is_valid_double(double x) {
    return (x == x) && (x != INFINITY) && (x != -INFINITY) &&
        (x < 1e10) && (x > -1e10);
}

// Write error log with detailed information
void write_error_log(FILE* logfile, int step, int particle_id,
    const char* location, double val1, double val2, double val3) {
    fprintf(logfile, "ERROR at step %d, particle %d, location: %s\n",
        step, particle_id, location);
    fprintf(logfile, "  Values: %lf, %lf, %lf\n", val1, val2, val3);
    fprintf(logfile, "  Valid: %d, %d, %d\n",
        is_valid_double(val1), is_valid_double(val2), is_valid_double(val3));
    fflush(logfile);
}

// Optimized periodic boundary function (inline for performance)
inline void apply_periodic_boundary(double& dx, double& dy, double xhalf, double yhalf, double xlen, double ylen) {
    if (dx > xhalf) dx -= xlen;
    else if (dx < -xhalf) dx += xlen;
    if (dy > yhalf) dy -= ylen;
    else if (dy < -yhalf) dy += ylen;
}

// initiate velocity and position of MPC particles
void confmpc()
{
    int i, j;
    double phi, theta;
    double m1 = 0;
    double m2 = 0;
    double m3 = 0;
    double m4;
    double m5;
    double m6;
    double m7;

    // Initial positions for snowman colloid (ball 1 active, ball 2 passive)
    rcx[1] = 16.0;
    rcy[1] = 10.0;
    rcz[1] = 6.0;

    rcx[2] = 16.0;
    rcy[2] = 10.0 - dis;
    rcz[2] = 6.0;

    rcx0[1] = 16.0;
    rcy0[1] = 10.0;
    rcz0[1] = 6.0;

    rcx0[2] = 16.0;
    rcy0[2] = 10.0 - dis;
    rcz0[2] = 6.0;

    // Sizes and interaction ranges for colloids
    sigma[1] = 2.0;
    rc[1] = sigma[1] * pow(2.0, 1 / 12.0);
    rn[1] = rc_att1 + 0.9;
    square_sigma1 = sigma[1] * sigma[1];
    square_rc[1] = rc[1] * rc[1];
    square_rca[1] = square_rc_att1;
    square_rn[1] = rn[1] * rn[1];

    sigma[2] = 3.0;
    rc[2] = sigma[2] * pow(2.0, 1 / 12.0);
    rn[2] = rc_att2 + 0.9;
    square_sigma2 = sigma[2] * sigma[2];
    square_rc[2] = rc[2] * rc[2];
    square_rca[2] = square_rc_att2;
    square_rn[2] = rn[2] * rn[2];

    square_bb[1] = square_sigma1;
    square_bb[2] = square_sigma2;
    square_rea = square_rd1;

    for (j = 1; j <= ncolloid; j++)
    {
        Mass[j] = 4 * Pi * sigma[j] * sigma[j] * sigma[j] * rho / 3;
    }

    for (j = 1; j <= ncolloid; j++)
    {
        I[j] = 2 * Mass[j] * sigma[j] * sigma[j] / 5;
    }

    for (j = 1; j <= ncolloid; j++)
    {
        Miu[j] = Mass[j] / (Mass[j] + 1);
    }

    // initialize velocities and forces for colloids
    for (j = 1; j <= ncolloid; j++)
    {
        vcx[j] = 0.0;
        vcy[j] = 0.0;
        vcz[j] = 0.0;
        fcx0[j] = 0.0;
        fcy0[j] = 0.0;
        fcz0[j] = 0.0;
        fcx[j] = 0.0;
        fcy[j] = 0.0;
        fcz[j] = 0.0;
        omegax[j] = 0.0;
        omegay[j] = 0.0;
        omegaz[j] = 0.0;
        Momx0[j] = 0.0;
        Momy0[j] = 0.0;
        Momz0[j] = 0.0;
        Momx[j] = 0.0;
        Momy[j] = 0.0;
        Momz[j] = 0.0;
    }

    // Set initial orientation vectors for the snowman (two perpendicular directions)
    nx[1][1] = 0.0;
    ny[1][1] = 1.0;
    nz[1][1] = 0.0;

    nx[2][1] = 0.0;
    ny[2][1] = 1.0;
    nz[2][1] = 0.0;

    nx[1][2] = 0.0;
    ny[1][2] = 0.0;
    nz[1][2] = 1.0;

    nx[2][2] = 0.0;
    ny[2][2] = 0.0;
    nz[2][2] = 1.0;

    // Wall interaction radii (for top/bottom boundaries)
    rwall[1] = 2.27;
    rcwall[1] = rwall[1] * pow(2.0, 1 / 48.0);
    square_rwall1 = rwall[1] * rwall[1];

    rwall[2] = 3.4;
    rcwall[2] = rwall[2] * pow(2.0, 1 / 48.0);
    square_rwall2 = rwall[2] * rwall[2];

    // (No large third ball, so no rball needed beyond index 2)
    rball[1] = 9.5;
    rcball[1] = rball[1] * pow(2.0, 1 / 48.0);
    square_rball1 = rball[1] * rball[1];
    square_rcball1 = rcball[1] * rcball[1];

    rball[2] = 10.5;
    rcball[2] = rball[2] * pow(2.0, 1 / 48.0);
    square_rball2 = rball[2] * rball[2];
    square_rcball2 = rcball[2] * rcball[2];

    // Randomly place MPC fluid particles avoiding colloid exclusion volumes
    for (i = 1; i <= nparmpc; i++)
    {
    loop:
        rxmpc[i] = xlen * r250();
        rympc[i] = ylen * r250();
        rzmpc[i] = zlen * r250();

        for (j = 1; j <= ncolloid; j++)
        {
            m4 = rxmpc[i] - rcx0[j];
            m5 = rympc[i] - rcy0[j];
            m6 = rzmpc[i] - rcz0[j];

            if (m4 > xhalf) m4 = m4 - xlen;
            else if (m4 < -xhalf) m4 = m4 + xlen;
            if (m5 > yhalf) m5 = m5 - ylen;
            else if (m5 < -yhalf) m5 = m5 + ylen;
            m7 = m4 * m4 + m5 * m5 + m6 * m6;

            if (m7 <= square_rca[j])
            {
                goto loop;
            }
        }

        phi = 2 * Pi * r250();
        theta = Pi * r250();
        vxmpc[i] = vavempc * sin(theta) * cos(phi);
        vympc[i] = vavempc * sin(theta) * sin(phi);
        vzmpc[i] = vavempc * cos(theta);

        m1 = m1 + vxmpc[i];
        m2 = m2 + vympc[i];
        m3 = m3 + vzmpc[i];
    }

    // Zero net momentum of fluid
    m1 = m1 / nparmpc;
    m2 = m2 / nparmpc;
    m3 = m3 / nparmpc;
    for (i = 1; i <= nparmpc; i++)
    {
        vxmpc[i] = vxmpc[i] - m1;
        vympc[i] = vympc[i] - m2;
        vzmpc[i] = vzmpc[i] - m3;
    }
}

// Lennard-Jones repulsive force for colloid 1 (active) with fluid
double rforce1(double r2)
{
    // Prevent force explosion at very small distances
    if (r2 < MIN_R2) {
        r2 = MIN_R2;
    }

    double repforce1, ir2, ir4, ir8, ir12, ir24;
    ir2 = square_sigma1 / r2;
    ir4 = ir2 * ir2;
    ir8 = ir4 * ir4;
    ir12 = ir4 * ir8;
    ir24 = ir12 * ir12;
    repforce1 = (96 * ir24 - 48 * ir12) / r2;

    // Limit force magnitude to prevent numerical overflow
    if (repforce1 > MAX_FORCE) repforce1 = MAX_FORCE;
    if (repforce1 < -MAX_FORCE) repforce1 = -MAX_FORCE;

    return(repforce1);
}

// Lennard-Jones attractive force for colloid 1 with fluid
double aforce1(double r2)
{
    // Prevent force explosion at very small distances
    if (r2 < MIN_R2) {
        r2 = MIN_R2;
    }

    double attforce1, r1, r3, ir2, ir4, ir8, ir12, ir24;

    if (r2 <= square_rd1)
    {
        ir2 = square_sigma1 / r2;
        ir4 = ir2 * ir2;
        ir8 = ir4 * ir4;
        ir12 = ir4 * ir8;
        ir24 = ir12 * ir12;
        attforce1 = (96 * ir24 - 48 * ir12) / r2;

        // Limit force magnitude to prevent numerical overflow
        if (attforce1 > MAX_FORCE) attforce1 = MAX_FORCE;
        if (attforce1 < -MAX_FORCE) attforce1 = -MAX_FORCE;

        return(attforce1);
    }

    else if (r2 > square_rd1 && r2 <= square_rbottem1)
    {
        // Ensure r2 is not too small before taking sqrt
        if (r2 < MIN_R2) r2 = MIN_R2;
        r1 = sqrt(r2);
        attforce1 = -49212.5469405786 * r1 + 204062.456617842 - 211503.826889518 / r1;

        // Limit force magnitude to prevent numerical overflow
        if (attforce1 > MAX_FORCE) attforce1 = MAX_FORCE;
        if (attforce1 < -MAX_FORCE) attforce1 = -MAX_FORCE;
    }

    else if (r2 > square_rbottem1 && r2 <= square_rc_att1)
    {
        // Ensure r2 is not too small before taking sqrt
        if (r2 < MIN_R2) r2 = MIN_R2;
        r1 = sqrt(r2);
        attforce1 = 6000.0 * r1 - 25800.0 + 27720.0 / r1;

        // Limit force magnitude to prevent numerical overflow
        if (attforce1 > MAX_FORCE) attforce1 = MAX_FORCE;
        if (attforce1 < -MAX_FORCE) attforce1 = -MAX_FORCE;
    }
    return(attforce1);
}

// Lennard-Jones repulsive force for colloid 2 (passive) with fluid
double rforce2(double r2)
{
    // Prevent force explosion at very small distances
    if (r2 < MIN_R2) {
        r2 = MIN_R2;
    }

    double repforce2, ir2, ir4, ir8, ir12, ir24;
    ir2 = square_sigma2 / r2;
    ir4 = ir2 * ir2;
    ir8 = ir4 * ir4;
    ir12 = ir4 * ir8;
    ir24 = ir12 * ir12;
    repforce2 = (96 * ir24 - 48 * ir12) / r2;

    // Limit force magnitude to prevent numerical overflow
    if (repforce2 > MAX_FORCE) repforce2 = MAX_FORCE;
    if (repforce2 < -MAX_FORCE) repforce2 = -MAX_FORCE;

    return(repforce2);
}

// Lennard-Jones attractive force for colloid 2 with fluid
double aforce2(double r2)
{
    // Prevent force explosion at very small distances
    if (r2 < MIN_R2) {
        r2 = MIN_R2;
    }

    double attforce2, r1, r3, ir2, ir4, ir8, ir12, ir24;

    if (r2 <= square_rd2)
    {
        ir2 = square_sigma2 / r2;
        ir4 = ir2 * ir2;
        ir8 = ir4 * ir4;
        ir12 = ir4 * ir8;
        ir24 = ir12 * ir12;
        attforce2 = (96 * ir24 - 48 * ir12) / r2;

        // Limit force magnitude to prevent numerical overflow
        if (attforce2 > MAX_FORCE) attforce2 = MAX_FORCE;
        if (attforce2 < -MAX_FORCE) attforce2 = -MAX_FORCE;

        return(attforce2);
    }

    else if (r2 > square_rd2 && r2 <= square_rbottem2)
    {
        // Ensure r2 is not too small before taking sqrt
        if (r2 < MIN_R2) r2 = MIN_R2;
        r1 = sqrt(r2);
        attforce2 = -27507.25362186395 * r1 + 173255.7886128438 - 272789.0562631062 / r1;

        // Limit force magnitude to prevent numerical overflow
        if (attforce2 > MAX_FORCE) attforce2 = MAX_FORCE;
        if (attforce2 < -MAX_FORCE) attforce2 = -MAX_FORCE;
    }

    else if (r2 > square_rbottem2 && r2 <= square_rc_att2)
    {
        // Ensure r2 is not too small before taking sqrt
        if (r2 < MIN_R2) r2 = MIN_R2;
        r1 = sqrt(r2);
        attforce2 = 3472.222222115325 * r1 - 22500.0 + 36437.5 / r1;

        // Limit force magnitude to prevent numerical overflow
        if (attforce2 > MAX_FORCE) attforce2 = MAX_FORCE;
        if (attforce2 < -MAX_FORCE) attforce2 = -MAX_FORCE;
    }
    return(attforce2);
}

// Wall repulsive potential for colloid 1 (bottom/top walls)
double rpotentialwall1(double r2)
{
    double reppotentialwall1, ir2, ir4, ir8, ir12, ir24, ir48, ir96;
    ir2 = square_rwall1 / r2;
    ir4 = ir2 * ir2;
    ir8 = ir4 * ir4;
    ir12 = ir4 * ir8;
    ir24 = ir12 * ir12;
    ir48 = ir24 * ir24;
    ir96 = ir48 * ir48;
    reppotentialwall1 = 4 * (ir96 - ir48) + 1;
    return(reppotentialwall1);
}

// Wall repulsive potential for colloid 2
double rpotentialwall2(double r2)
{
    double reppotentialwall2, ir2, ir4, ir8, ir12, ir24, ir48, ir96;
    ir2 = square_rwall2 / r2;
    ir4 = ir2 * ir2;
    ir8 = ir4 * ir4;
    ir12 = ir4 * ir8;
    ir24 = ir12 * ir12;
    ir48 = ir24 * ir24;
    ir96 = ir48 * ir48;
    reppotentialwall2 = 4 * (ir96 - ir48) + 1;
    return(reppotentialwall2);
}

// Wall repulsive force for colloid 1
double rforcewall1(double r2)
{
    // 防止 r2 过小导致的除零或溢出
    if (r2 < 1e-12) {
        r2 = 1e-12;
    }

    double repforcewall1, ir2, ir4, ir8, ir12, ir24, ir48, ir96;
    ir2 = square_rwall1 / r2;
    ir4 = ir2 * ir2;
    ir8 = ir4 * ir4;
    ir12 = ir4 * ir8;
    ir24 = ir12 * ir12;
    ir48 = ir24 * ir24;
    ir96 = ir48 * ir48;
    repforcewall1 = (384 * ir96 - 192 * ir48) / r2;
    return(repforcewall1);
}

// Wall repulsive force for colloid 2
double rforcewall2(double r2)
{
    // 防止 r2 过小导致的除零或溢出
    if (r2 < 1e-12) {
        r2 = 1e-12;
    }

    double repforcewall2, ir2, ir4, ir8, ir12, ir24, ir48, ir96;
    ir2 = square_rwall2 / r2;
    ir4 = ir2 * ir2;
    ir8 = ir4 * ir4;
    ir12 = ir4 * ir8;
    ir24 = ir12 * ir12;
    ir48 = ir24 * ir24;
    ir96 = ir48 * ir48;
    repforcewall2 = (384 * ir96 - 192 * ir48) / r2;
    return(repforcewall2);
}

/*main function*/
int main()
{
    int nseed;                      // random number seed
    nseed = (unsigned)time(NULL);
    wmup_ks(nseed);
    FILE* fpw1;
    confmpc();
    int i, j, k, l, m;
    int inside[ncolloid + 1];
    double record1[ncolloid + 1];
    static int record2[nparmpc + 1];
    static int record3[nparmpc + 1];
    static int kind[nparmpc + 1];
    double bx = xlen / nxcell;
    double by = ylen / nycell;
    double bz = zlen / nzcell;
    double p1 = 1.0;
    double p2 = 1.0 / 1000.0;
    double cosphi;
    static int countcellmpc[nxcell * nycell * (nzcell + 2) + 1];  // number of particle in each collision box
    int me1, me2, me3, me4; // medium variable
    double med1, med2, med3, med4, med5, med6, med7, med8;
    static double vxmc[nxcell * nycell * (nzcell + 2) + 1];  // center of mass velocity
    static double vymc[nxcell * nycell * (nzcell + 2) + 1];
    static double vzmc[nxcell * nycell * (nzcell + 2) + 1];
    static double dirx[nxcell * nycell * (nzcell + 2) + 1];  // rotational axis
    static double diry[nxcell * nycell * (nzcell + 2) + 1];
    static double dirz[nxcell * nycell * (nzcell + 2) + 1];
    double vx, vy, vz; //relative velocity
    double vpx, vpy, vpz; //relative velocity parallel to random vector
    double vvx, vvy, vvz; //relative velocity vertical to random vector
    static int cellmpc[nparmpc + 1]; // the collision box that a MPC particle belongs to
    double xshift, yshift, zshift;
    static int npar[nparmpc + 1];
    static int mpar[ncolloid + 1][a + 1];
    double E, E1, E2, E3, E4, E5;
    static double fx[nparmpc + 1];
    static double fy[nparmpc + 1];
    static double fz[nparmpc + 1];
    double Ebef, Eaft;
    double Ebef1, Eaft1;
    double Ebef2, Eaft2;
    double Ebef3, Eaft3;
    double Ebef4, Eaft4;
    double Ebef5, Eaft5;
    double nx_st, ny_st, nz_st;
    double vrx, vry, vrz;
    double vrvx, vrvy, vrvz;
    double vrpx, vrpy, vrpz;
    double dt;

    // Time tracking variables for progress monitoring
    clock_t start_time, current_time, last_progress_time;
    double elapsed_seconds, estimated_total_seconds, remaining_seconds;
    double time_per_step;
    int progress_interval = 1000;  // Print progress every N steps

    fpw1 = fopen("a1.txt", "w");
    if (fpw1 == NULL)
    {
        printf("Cannot open this file!! ");
        exit(0);
    }

    // Write Tecplot format header
    fprintf(fpw1, "Variables = \"x\" \"y\" \"z\" \"Type\" \"Radius\"\n");
    fflush(fpw1);

    // Initialize error log file
    FILE* error_log = fopen("error_log.txt", "w");
    if (error_log == NULL)
    {
        printf("Cannot open error log file!!\n");
        exit(0);
    }

    // Print startup information and configuration
    printf("\n========================================\n");
    printf("MPC Snowman Simulation Started\n");
    printf("========================================\n");
    printf("Configuration:\n");
    printf("  Total particles (MPC): %d\n", nparmpc);
    printf("  Total simulation steps: %d\n", eqtime);
    printf("  Colloid number: %d\n", ncolloid);
    printf("  Output interval: every 20 steps\n");
    printf("  Progress update: every %d steps\n", progress_interval);

    // Calculate total operations and estimate time
    long long total_operations = (long long)eqtime * (long long)nparmpc;
    printf("  Total operations: %.2e\n", (double)total_operations);

    // Rough time estimation (assuming ~1e6 operations per second, adjust based on your system)
    double estimated_hours = (double)total_operations / 1e6 / 3600.0;
    if (estimated_hours > 1.0) {
        printf("  Estimated runtime: ~%.1f hours (rough estimate)\n", estimated_hours);
    }
    else {
        double estimated_minutes = estimated_hours * 60.0;
        printf("  Estimated runtime: ~%.1f minutes (rough estimate)\n", estimated_minutes);
    }

    // Provide parameter suggestions if runtime is too long
    if (estimated_hours > 1.0) {
        printf("\n*** Parameter Optimization Suggestion ***\n");
        printf("Estimated runtime is long. For quick testing, consider:\n");
        printf("  - Reduce steps: #define eqtime 10000  (for quick test)\n");
        printf("  - Reduce particles: #define nparmpc 20480  (10x faster)\n");
        printf("  - Combined quick test: eqtime=10000, nparmpc=20480\n");
        printf("  This would reduce runtime to ~%.1f minutes\n",
            (double)10000 * 20480 / 1e6 / 60.0);
    }
    printf("========================================\n\n");

    // initialize collision cell count
    for (i = 1; i <= nxcell * nycell * (nzcell + 2); i++)
    {
        countcellmpc[i] = 0;
    }

    // set all fluid particles initially as type 1 (reactant)
    for (j = 1; j <= nparmpc; j++)
    {
        kind[j] = 1;
    }

    // Record start time and begin simulation
    start_time = clock();
    last_progress_time = start_time;
    printf("Simulation started...\n");
    printf("Progress will be updated every %d steps.\n\n", progress_interval);
    fflush(stdout);

    // SRD simulation main loop
    for (i = 1; i <= eqtime; i++)
    {
        // reset forces and markers (optimized: only reset when needed)
        // Note: forces will be reset when used, so we can skip full reset
        // Only reset markers that are actually checked
        for (j = 1; j <= nparmpc; j++)
        {
            npar[j] = 0;
            record2[j] = 0;
            // record3 and forces will be set when needed, no need to reset here
        }

        // Reset colloid interaction lists
        inside[1] = 0;
        inside[2] = 0;
        // mpar arrays will be filled from scratch, no need to clear all

        // Assign fluid particles to colloid interaction lists (optimized: single loop for all colloids)
        for (j = 1; j <= nparmpc; j++)
        {
            // Check distance to colloid 1
            med1 = rxmpc[j] - rcx[1];
            med2 = rympc[j] - rcy[1];
            med3 = rzmpc[j] - rcz[1];
            apply_periodic_boundary(med1, med2, xhalf, yhalf, xlen, ylen);
            med4 = med1 * med1 + med2 * med2 + med3 * med3;

            if (med4 <= square_rn[1])
            {
                npar[j] = 1;
                inside[1] = inside[1] + 1;
                mpar[1][inside[1]] = j;
            }

            // Check distance to colloid 2
            med1 = rxmpc[j] - rcx[2];
            med2 = rympc[j] - rcy[2];
            med3 = rzmpc[j] - rcz[2];
            apply_periodic_boundary(med1, med2, xhalf, yhalf, xlen, ylen);
            med4 = med1 * med1 + med2 * med2 + med3 * med3;

            if (med4 <= square_rn[2])
            {
                npar[j] = 1;
                inside[2] = inside[2] + 1;
                mpar[2][inside[2]] = j;
            }
        }

        // === Check validity of all MPC particles before free streaming ===
#if ENABLE_NUMERICAL_CHECKS
        for (j = 1; j <= nparmpc; j++) {
            if (!is_valid_double(rxmpc[j]) || !is_valid_double(rympc[j]) ||
                !is_valid_double(rzmpc[j]) || !is_valid_double(vxmpc[j]) ||
                !is_valid_double(vympc[j]) || !is_valid_double(vzmpc[j])) {
                write_error_log(error_log, i, j, "Before free streaming",
                    rxmpc[j], rympc[j], rzmpc[j]);
                fprintf(error_log, "  Velocities: %lf, %lf, %lf\n",
                    vxmpc[j], vympc[j], vzmpc[j]);
                fprintf(error_log, "  CRITICAL: Invalid particle state detected!\n");
                fflush(error_log);
                fclose(error_log);
                exit(1);
            }
        }
#endif

        // Free streaming of MPC particles (no collision with colloids)
        for (j = 1; j <= nparmpc; j++)
        {
            if (npar[j] == 0)
            {
                if (rzmpc[j] >= 0.9 && rzmpc[j] <= (zlen - 0.9))
                {
                    rxmpc[j] = rxmpc[j] + vxmpc[j] * dtmpc;
                    rympc[j] = rympc[j] + vympc[j] * dtmpc;
                    rzmpc[j] = rzmpc[j] + vzmpc[j] * dtmpc;
                }
                else if (rzmpc[j] < 0.9)
                {
                    if (vzmpc[j] >= 0.0)
                    {
                        rxmpc[j] = rxmpc[j] + vxmpc[j] * dtmpc;
                        rympc[j] = rympc[j] + vympc[j] * dtmpc;
                        rzmpc[j] = rzmpc[j] + vzmpc[j] * dtmpc;
                    }
                    else if (vzmpc[j] < 0.0)
                    {
                        // Prevent divide-by-zero: check if vzmpc is too small
                        if (fabs(vzmpc[j]) < 1e-10) {
                            // Use default small time step if velocity is too small
                            dt = dtmpc;
                        }
                        else {
                            dt = -rzmpc[j] / vzmpc[j];
                        }
                        if (dt > dtmpc)
                        {
                            rxmpc[j] = rxmpc[j] + vxmpc[j] * dtmpc;
                            rympc[j] = rympc[j] + vympc[j] * dtmpc;
                            rzmpc[j] = rzmpc[j] + vzmpc[j] * dtmpc;
                        }
                        else if (dt <= dtmpc)
                        {
                            rxmpc[j] = rxmpc[j] + vxmpc[j] * dt;
                            rympc[j] = rympc[j] + vympc[j] * dt;
                            rzmpc[j] = rzmpc[j] + vzmpc[j] * dt;
                            vxmpc[j] = -vxmpc[j];
                            vympc[j] = -vympc[j];
                            vzmpc[j] = -vzmpc[j];
                            rxmpc[j] = rxmpc[j] + vxmpc[j] * (dtmpc - dt);
                            rympc[j] = rympc[j] + vympc[j] * (dtmpc - dt);
                            rzmpc[j] = rzmpc[j] + vzmpc[j] * (dtmpc - dt);
                        }
                    }
                }
                else if (rzmpc[j] > (zlen - 0.9))
                {
                    if (vzmpc[j] <= 0.0)
                    {
                        rxmpc[j] = rxmpc[j] + vxmpc[j] * dtmpc;
                        rympc[j] = rympc[j] + vympc[j] * dtmpc;
                        rzmpc[j] = rzmpc[j] + vzmpc[j] * dtmpc;
                    }
                    else if (vzmpc[j] > 0.0)
                    {
                        // Prevent divide-by-zero: check if vzmpc is too small
                        if (fabs(vzmpc[j]) < 1e-10) {
                            // Use default small time step if velocity is too small
                            dt = dtmpc;
                        }
                        else {
                            dt = (zlen - rzmpc[j]) / vzmpc[j];
                        }
                        if (dt > dtmpc)
                        {
                            rxmpc[j] = rxmpc[j] + vxmpc[j] * dtmpc;
                            rympc[j] = rympc[j] + vympc[j] * dtmpc;
                            rzmpc[j] = rzmpc[j] + vzmpc[j] * dtmpc;
                        }
                        else if (dt <= dtmpc)
                        {
                            rxmpc[j] = rxmpc[j] + vxmpc[j] * dt;
                            rympc[j] = rympc[j] + vympc[j] * dt;
                            rzmpc[j] = rzmpc[j] + vzmpc[j] * dt;
                            vxmpc[j] = -vxmpc[j];
                            vympc[j] = -vympc[j];
                            vzmpc[j] = -vzmpc[j];
                            rxmpc[j] = rxmpc[j] + vxmpc[j] * (dtmpc - dt);
                            rympc[j] = rympc[j] + vympc[j] * (dtmpc - dt);
                            rzmpc[j] = rzmpc[j] + vzmpc[j] * (dtmpc - dt);
                        }
                    }
                }
            }
        }

        // === 在这里补充 x, y 的周期边界 ===
        for (j = 1; j <= nparmpc; j++)
        {
            // x 方向周期
            while (rxmpc[j] >= xlen) rxmpc[j] -= xlen;
            while (rxmpc[j] < 0.0)  rxmpc[j] += xlen;

            // y 方向周期
            while (rympc[j] >= ylen) rympc[j] -= ylen;
            while (rympc[j] < 0.0)  rympc[j] += ylen;

#if ENABLE_NUMERICAL_CHECKS
            // Check validity after periodic boundary treatment
            if (!is_valid_double(rxmpc[j]) || !is_valid_double(rympc[j]) ||
                !is_valid_double(rzmpc[j])) {
                write_error_log(error_log, i, j, "After periodic boundary",
                    rxmpc[j], rympc[j], rzmpc[j]);
            }
#endif
        }


        // SRD collision step: assign random shifts to grid
        xshift = r250() - 0.5;
        yshift = r250() - 0.5;
        zshift = r250() - 0.5;

        for (j = 1; j <= nxcell * nycell * (nzcell + 2); j++)
        {
            countcellmpc[j] = 0;
        }

        for (j = 1; j <= nparmpc; j++)
        {
            double raw_x = rxmpc[j] + xshift;
            double raw_y = rympc[j] + yshift;
            double raw_z = rzmpc[j] + zshift;

#if ENABLE_NUMERICAL_CHECKS
            // Check validity of position and shift values before calculation
            if (!is_valid_double(raw_x) || !is_valid_double(raw_y) ||
                !is_valid_double(raw_z)) {
                write_error_log(error_log, i, j, "Before cell assignment",
                    raw_x, raw_y, raw_z);
                fprintf(error_log, "  Original pos: %lf, %lf, %lf\n",
                    rxmpc[j], rympc[j], rzmpc[j]);
                fprintf(error_log, "  Shift: %lf, %lf, %lf\n",
                    xshift, yshift, zshift);
                fflush(error_log);
                fclose(error_log);
                exit(1);
            }
#endif

            // Improved periodic boundary handling using modulo arithmetic
            me1 = int(raw_x + 1);
            if (me1 > nxcell) {
                me1 = ((me1 - 1) % nxcell) + 1;
            }
            if (me1 <= 0) {
                me1 = ((me1 - 1) % nxcell + nxcell) % nxcell + 1;
            }

            me2 = int(raw_y + 1);
            if (me2 > nycell) {
                me2 = ((me2 - 1) % nycell) + 1;
            }
            if (me2 <= 0) {
                me2 = ((me2 - 1) % nycell + nycell) % nycell + 1;
            }

            me3 = int(raw_z + 1);
            if (me3 > nzcell + 1) {
                me3 = ((me3 - 1) % (nzcell + 2)) + 1;
            }
            if (me3 <= 0) {
                me3 = ((me3 - 1) % (nzcell + 2) + (nzcell + 2)) % (nzcell + 2) + 1;
            }

            // Validate me1, me2, me3 are in valid ranges
            if (me1 < 1 || me1 > nxcell || me2 < 1 || me2 > nycell ||
                me3 < 1 || me3 >(nzcell + 2)) {
                fprintf(error_log, "ERROR at step %d, particle %d: Invalid cell indices\n",
                    i, j);
                fprintf(error_log, "  me1=%d (should be 1-%d), me2=%d (should be 1-%d), me3=%d (should be 1-%d)\n",
                    me1, nxcell, me2, nycell, me3, nzcell + 2);
                fprintf(error_log, "  raw_x=%lf, raw_y=%lf, raw_z=%lf\n",
                    raw_x, raw_y, raw_z);
                fprintf(error_log, "  rxmpc=%lf, rympc=%lf, rzmpc=%lf\n",
                    rxmpc[j], rympc[j], rzmpc[j]);
                fflush(error_log);
                fclose(error_log);
                exit(1);
            }

            me4 = (me3 - 1) * nxcell * nycell + (me2 - 1) * nxcell + me1;

            // Validate me4 is in valid range
            int max_cell = nxcell * nycell * (nzcell + 2);
            if (me4 < 1 || me4 > max_cell) {
                fprintf(error_log, "ERROR at step %d, particle %d: Invalid me4 index\n",
                    i, j);
                fprintf(error_log, "  me4=%d (should be 1-%d)\n", me4, max_cell);
                fprintf(error_log, "  me1=%d, me2=%d, me3=%d\n", me1, me2, me3);
                fprintf(error_log, "  Calculation: (%d-1)*%d*%d + (%d-1)*%d + %d = %d\n",
                    me3, nxcell, nycell, me2, nxcell, me1, me4);
                fflush(error_log);
                fclose(error_log);
                exit(1);
            }

            cellmpc[j] = me4;
            countcellmpc[me4] = countcellmpc[me4] + 1;
        }

        for (j = 1; j <= nxcell * nycell * (nzcell + 2); j++)
        {
            if (countcellmpc[j] > 0)
            {
                vxmc[j] = 0.0;
                vymc[j] = 0.0;
                vzmc[j] = 0.0;
            }
        }

        for (j = 1; j <= nparmpc; j++)
        {
            vxmc[cellmpc[j]] = vxmc[cellmpc[j]] + vxmpc[j];
            vymc[cellmpc[j]] = vymc[cellmpc[j]] + vympc[j];
            vzmc[cellmpc[j]] = vzmc[cellmpc[j]] + vzmpc[j];
        }

        for (j = 1; j <= nxcell * nycell * (nzcell + 2); j++)
        {
            if (countcellmpc[j] > 0)
            {
                vxmc[j] = vxmc[j] / countcellmpc[j];
                vymc[j] = vymc[j] / countcellmpc[j];
                vzmc[j] = vzmc[j] / countcellmpc[j];
                // choose random rotation axis
                phi = 2 * Pi * r250();
                theta = ang;
                dirx[j] = sin(theta) * cos(phi);
                diry[j] = sin(theta) * sin(phi);
                dirz[j] = cos(theta);
            }
        }

        for (j = 1; j <= nparmpc; j++)
        {
            // relative velocity to cell center of mass
            vx = vxmpc[j] - vxmc[cellmpc[j]];
            vy = vympc[j] - vymc[cellmpc[j]];
            vz = vzmpc[j] - vzmc[cellmpc[j]];
            // parallel component of velocity
            vpx = vx * dirx[cellmpc[j]] + vy * diry[cellmpc[j]] + vz * dirz[cellmpc[j]];
            // perpendicular component of velocity
            vvx = vx - vpx * dirx[cellmpc[j]];
            vvy = vy - vpx * diry[cellmpc[j]];
            vvz = vz - vpx * dirz[cellmpc[j]];
            // rotate perpendicular component
            vxmpc[j] = vxmc[cellmpc[j]] + vpx * dirx[cellmpc[j]] + vvy * cosang + vvz * sinang;
            vympc[j] = vymc[cellmpc[j]] + vpx * diry[cellmpc[j]] + vvz * cosang - vvy * sinang;
            vzmpc[j] = vzmc[cellmpc[j]] + vpx * dirz[cellmpc[j]] + vx * sinang + vy * cosang;
        }

        // Chemical reaction: ball 1 (active) can convert fluid particles of type 1 to type 2
        for (m = 1; m <= a; m++)
        {
            if (mpar[1][m] == 0) break;
            j = mpar[1][m];

            med1 = rxmpc[j] - rcx[1];
            med2 = rympc[j] - rcy[1];
            med3 = rzmpc[j] - rcz[1];

            if (med1 > xhalf) med1 = med1 - xlen;
            else if (med1 < -xhalf) med1 = med1 + xlen;
            if (med2 > yhalf) med2 = med2 - ylen;
            else if (med2 < -yhalf) med2 = med2 + ylen;
            med4 = med1 * med1 + med2 * med2 + med3 * med3;

            if (med4 <= square_rea)
            {
                if (kind[j] == 1)
                {
                    if (r250() <= p1)
                    {
                        record3[j] = 1;
                    }
                }
            }
        }

        for (m = 1; m <= a; m++)
        {
            if (mpar[1][m] == 0) break;
            j = mpar[1][m];

            med5 = rxmpc[j] - rcx[2];
            med6 = rympc[j] - rcy[2];
            med7 = rzmpc[j] - rcz[2];
            if (med5 > xhalf) med5 = med5 - xlen;
            else if (med5 < -xhalf) med5 = med5 + xlen;
            if (med6 > yhalf) med6 = med6 - ylen;
            else if (med6 < -yhalf) med6 = med6 + ylen;
            med8 = med5 * med5 + med6 * med6 + med7 * med7;

            // Reaction completion: if particle was activated and is far from ball 2, convert to product
            if (med8 > square_rc_att2)
            {
                if (record3[j] == 1)
                {
                    kind[j] = 2;
                    record3[j] = 0;
                }
            }
        }

        // Calculate new forces on fluid particles from colloids
        for (m = 1; m <= a; m++)
        {
            if (mpar[1][m] == 0) break;
            j = mpar[1][m];

            med1 = rxmpc[j] - rcx[1];
            med2 = rympc[j] - rcy[1];
            med3 = rzmpc[j] - rcz[1];

            if (med1 > xhalf) med1 = med1 - xlen;
            else if (med1 < -xhalf) med1 = med1 + xlen;
            if (med2 > yhalf) med2 = med2 - ylen;
            else if (med2 < -yhalf) med2 = med2 + ylen;
            med4 = med1 * med1 + med2 * med2 + med3 * med3;

            // Check if distance is too small and log warning
            if (med4 < MIN_R2) {
                fprintf(error_log, "WARNING at step %d, particle %d, colloid 1: Distance too small (r2=%le), using MIN_R2\n",
                    i, j, med4);
                fprintf(error_log, "  Particle position: %lf, %lf, %lf\n", rxmpc[j], rympc[j], rzmpc[j]);
                fprintf(error_log, "  Colloid position: %lf, %lf, %lf\n", rcx[1], rcy[1], rcz[1]);
                fprintf(error_log, "  Distance vector: %lf, %lf, %lf\n", med1, med2, med3);
                fflush(error_log);
                med4 = MIN_R2;  // Use minimum distance to prevent force explosion
            }

            if (kind[j] == 1)
            {
                if (med4 <= square_rc[1])
                {
                    med5 = rforce1(med4);
                    med6 = med1 * med5;
                    med7 = med2 * med5;
                    med8 = med3 * med5;
                    fxmpc[j] = fxmpc[j] + med6;
                    fympc[j] = fympc[j] + med7;
                    fzmpc[j] = fzmpc[j] + med8;
                    fcx[1] = fcx[1] - med6;
                    fcy[1] = fcy[1] - med7;
                    fcz[1] = fcz[1] - med8;
                }
            }
            else if (kind[j] == 2)
            {
                if (med4 <= square_rca[1])
                {
                    med5 = aforce1(med4);
                    med6 = med1 * med5;
                    med7 = med2 * med5;
                    med8 = med3 * med5;
                    fxmpc[j] = fxmpc[j] + med6;
                    fympc[j] = fympc[j] + med7;
                    fzmpc[j] = fzmpc[j] + med8;
                    fcx[1] = fcx[1] - med6;
                    fcy[1] = fcy[1] - med7;
                    fcz[1] = fcz[1] - med8;
                }
            }
        }

        for (m = 1; m <= a; m++)
        {
            if (mpar[2][m] == 0) break;
            j = mpar[2][m];

            med1 = rxmpc[j] - rcx[2];
            med2 = rympc[j] - rcy[2];
            med3 = rzmpc[j] - rcz[2];

            if (med1 > xhalf) med1 = med1 - xlen;
            else if (med1 < -xhalf) med1 = med1 + xlen;
            if (med2 > yhalf) med2 = med2 - ylen;
            else if (med2 < -yhalf) med2 = med2 + ylen;
            med4 = med1 * med1 + med2 * med2 + med3 * med3;

            // Check if distance is too small and log warning
            if (med4 < MIN_R2) {
                fprintf(error_log, "WARNING at step %d, particle %d, colloid 2: Distance too small (r2=%le), using MIN_R2\n",
                    i, j, med4);
                fprintf(error_log, "  Particle position: %lf, %lf, %lf\n", rxmpc[j], rympc[j], rzmpc[j]);
                fprintf(error_log, "  Colloid position: %lf, %lf, %lf\n", rcx[2], rcy[2], rcz[2]);
                fprintf(error_log, "  Distance vector: %lf, %lf, %lf\n", med1, med2, med3);
                fflush(error_log);
                med4 = MIN_R2;  // Use minimum distance to prevent force explosion
            }

            if (kind[j] == 1)
            {
                if (med4 <= square_rc[2])
                {
                    med5 = rforce2(med4);
                    med6 = med1 * med5;
                    med7 = med2 * med5;
                    med8 = med3 * med5;
                    fxmpc[j] = fxmpc[j] + med6;
                    fympc[j] = fympc[j] + med7;
                    fzmpc[j] = fzmpc[j] + med8;
                    fcx[2] = fcx[2] - med6;
                    fcy[2] = fcy[2] - med7;
                    fcz[2] = fcz[2] - med8;
                }
            }
            else if (kind[j] == 2)
            {
                if (med4 <= square_rca[2])
                {
                    med5 = aforce2(med4);
                    med6 = med1 * med5;
                    med7 = med2 * med5;
                    med8 = med3 * med5;
                    fxmpc[j] = fxmpc[j] + med6;
                    fympc[j] = fympc[j] + med7;
                    fzmpc[j] = fzmpc[j] + med8;
                    fcx[2] = fcx[2] - med6;
                    fcy[2] = fcy[2] - med7;
                    fcz[2] = fcz[2] - med8;
                }
            }
        }

        // === Check validity of forces after force calculation ===
#if ENABLE_NUMERICAL_CHECKS
        for (j = 1; j <= nparmpc; j++) {
            if (!is_valid_double(fxmpc[j]) || !is_valid_double(fympc[j]) ||
                !is_valid_double(fzmpc[j])) {
                write_error_log(error_log, i, j, "After force calculation",
                    fxmpc[j], fympc[j], fzmpc[j]);
            }
        }
#endif

        // Integrate fluid particle velocities (avoid double counting via record2)
        for (k = 1; k <= ncolloid; k++)
        {
            for (m = 1; m <= a; m++)
            {
                if (mpar[k][m] == 0) break;
                j = mpar[k][m];

                if (record2[j] == 1) continue;
                else if (record2[j] == 0)
                {
                    vxmpc[j] = vxmpc[j] + dtmd * (fxmpc[j] + fx[j]) / 2;
                    vympc[j] = vympc[j] + dtmd * (fympc[j] + fy[j]) / 2;
                    vzmpc[j] = vzmpc[j] + dtmd * (fzmpc[j] + fz[j]) / 2;
                    record2[j] = 1;
                }
            }
        }

        for (k = 1; k <= ncolloid; k++)
        {
            for (m = 1; m <= a; m++)
            {
                if (mpar[k][m] == 0) break;
                j = mpar[k][m];
                record2[j] = 0;
            }
        }

        // Rigid constraint between ball 1 and ball 2 (maintain snowman structure)
        med1 = rcx[1] - rcx[2];
        med2 = rcy[1] - rcy[2];
        med3 = rcz[1] - rcz[2];

        if (med1 > xhalf) med1 = med1 - xlen;
        else if (med1 < -xhalf) med1 = med1 + xlen;
        if (med2 > yhalf) med2 = med2 - ylen;
        else if (med2 < -yhalf) med2 = med2 + ylen;
        med4 = med1 * med1 + med2 * med2 + med3 * med3;
        med5 = sqrt(med4);

        // Prevent divide-by-zero: if med5 is too small, use default direction
        if (med5 < 1e-10) {
            fprintf(error_log, "WARNING at step %d: med5 too small (%le), using default direction\n",
                i, med5);
            fflush(error_log);
            med5 = 1.0;
            nx_st = 0.0;
            ny_st = 1.0;
            nz_st = 0.0;
        }
        else {
            // Apply spring force to maintain distance 'dis' between ball1 and ball2
            nx_st = med1 / med5;
            ny_st = med2 / med5;
            nz_st = med3 / med5;
        }

        fcx[1] = fcx[1] - med1 * k1 * (med5 - dis) / med5;
        fcy[1] = fcy[1] - med2 * k1 * (med5 - dis) / med5;
        fcz[1] = fcz[1] - med3 * k1 * (med5 - dis) / med5;
        fcx[2] = fcx[2] + med1 * k1 * (med5 - dis) / med5;
        fcy[2] = fcy[2] + med2 * k1 * (med5 - dis) / med5;
        //fcz[2] = fcx[2] + med3 * k1 * (med5 - dis) / med5;
        fcz[2] = fcz[2] + med3 * k1 * (med5 - dis) / med5;

        // Apply torque to enforce orientation constraints (keep orientation vectors aligned)
        for (k = 1; k <= 2; k++)
        {
            costheta[k][1] = nx[k][1] * nx_st + ny[k][1] * ny_st + nz[k][1] * nz_st;
            theta_st[k][1] = acos(costheta[k][1]);
            torque[k][1] = k1 * theta_st[k][1];

            nx_tor[k][1] = ny[k][1] * nz_st - nz[k][1] * ny_st;
            ny_tor[k][1] = nz[k][1] * nx_st - nx[k][1] * nz_st;
            nz_tor[k][1] = nx[k][1] * ny_st - ny[k][1] * nx_st;

            torquex[k][1] = torque[k][1] * nx_tor[k][1];
            torquey[k][1] = torque[k][1] * ny_tor[k][1];
            torquez[k][1] = torque[k][1] * nz_tor[k][1];

            omegax[k] = omegax[k] + dtmd * torquex[k][1] / (2 * I[k]);
            omegay[k] = omegay[k] + dtmd * torquey[k][1] / (2 * I[k]);
            omegaz[k] = omegaz[k] + dtmd * torquez[k][1] / (2 * I[k]);
        }

        costheta[1][2] = nx[1][2] * nx[2][2] + ny[1][2] * ny[2][2] + nz[1][2] * nz[2][2];
        theta_st[1][2] = acos(costheta[1][2]);
        torque[1][2] = k1 * theta_st[1][2];

        nx_tor[1][2] = ny[1][2] * nz[2][2] - nz[1][2] * ny[2][2];
        ny_tor[1][2] = nz[1][2] * nx[2][2] - nx[1][2] * nz[2][2];
        nz_tor[1][2] = nx[1][2] * ny[2][2] - ny[1][2] * nx[2][2];

        torquex[1][2] = torque[1][2] * nx_tor[1][2];
        torquey[1][2] = torque[1][2] * ny_tor[1][2];
        torquez[1][2] = torque[1][2] * nz_tor[1][2];
        torquex[2][2] = -torque[1][2] * nx_tor[1][2];
        torquey[2][2] = -torque[1][2] * ny_tor[1][2];
        torquez[2][2] = -torque[1][2] * nz_tor[1][2];

        for (k = 1; k <= 2; k++)
        {
            omegax[k] = omegax[k] + dtmd * torquex[k][2] / (2 * I[k]);
            omegay[k] = omegay[k] + dtmd * torquey[k][2] / (2 * I[k]);
            omegaz[k] = omegaz[k] + dtmd * torquez[k][2] / (2 * I[k]);
        }

        // Wall bounce-back conditions for colloids
        if (rcz[1] < rcwall[1])
        {
            fcz[1] = fcz[1] + rcz[1] * rforcewall1(rcz[1] * rcz[1]);
        }
        else if (rcz[1] >= (zlen - rcwall[1]))
        {
            fcz[1] = fcz[1] - (zlen - rcz[1]) * rforcewall1((zlen - rcz[1]) * (zlen - rcz[1]));
        }

        if (rcz[2] < rcwall[2])
        {
            fcz[2] = fcz[2] + rcz[2] * rforcewall2(rcz[2] * rcz[2]);
        }
        else if (rcz[2] >= (zlen - rcwall[2]))
        {
            fcz[2] = fcz[2] - (zlen - rcz[2]) * rforcewall2((zlen - rcz[2]) * (zlen - rcz[2]));
        }

        // Integrate colloid translational and rotational motion (velocity Verlet half-step)
        for (k = 1; k <= ncolloid; k++)
        {
            vcx[k] = vcx[k] + dtmd * (fcx[k] + fcx0[k]) / (2 * Mass[k]);
            vcy[k] = vcy[k] + dtmd * (fcy[k] + fcy0[k]) / (2 * Mass[k]);
            vcz[k] = vcz[k] + dtmd * (fcz[k] + fcz0[k]) / (2 * Mass[k]);
            omegax[k] = omegax[k] + dtmd * (Momx[k] + Momx0[k]) / (2 * I[k]);
            omegay[k] = omegay[k] + dtmd * (Momy[k] + Momy0[k]) / (2 * I[k]);
            omegaz[k] = omegaz[k] + dtmd * (Momz[k] + Momz0[k]) / (2 * I[k]);
        }

        // Bounce-back boundary condition for fluid on colloid 2 (no-slip on passive sphere)
        for (m = 1; m <= a; m++)
        {
            if (mpar[2][m] == 0) break;
            j = mpar[2][m];

            if (kind[j] == 1)
            {
                med1 = rxmpc[j] - rcx[2];
                med2 = rympc[j] - rcy[2];
                med3 = rzmpc[j] - rcz[2];

                if (med1 > xhalf) med1 = med1 - xlen;
                else if (med1 < -xhalf) med1 = med1 + xlen;
                if (med2 > yhalf) med2 = med2 - ylen;
                else if (med2 < -yhalf) med2 = med2 + ylen;
                med4 = med1 * med1 + med2 * med2 + med3 * med3;

                cosphi = med1 * (vxmpc[j] - vcx[2]) + med2 * (vympc[j] - vcy[2]) + med3 * (vzmpc[j] - vcz[2]);

                if (med4 <= square_bb[2] && cosphi <= 0)
                {
                    // Prevent divide-by-zero: check if med4 is too small
                    if (med4 < 1e-10) {
                        fprintf(error_log, "WARNING at step %d, particle %d: med4 too small (%le), skipping bounce-back\n",
                            i, j, med4);
                        fflush(error_log);
                        continue;
                    }

                    vrx = vxmpc[j] - vcx[2] - (omegay[2] * med3 - omegaz[2] * med2);
                    vry = vympc[j] - vcy[2] - (omegaz[2] * med1 - omegax[2] * med3);
                    vrz = vzmpc[j] - vcz[2] - (omegax[2] * med2 - omegay[2] * med1);

                    vrvx = (med1 * vrx + med2 * vry + med3 * vrz) * med1 / med4;
                    vrvy = (med1 * vrx + med2 * vry + med3 * vrz) * med2 / med4;
                    vrvz = (med1 * vrx + med2 * vry + med3 * vrz) * med3 / med4;

                    vrpx = vrx - vrvx;
                    vrpy = vry - vrvy;
                    vrpz = vrz - vrvz;

                    vxmpc[j] = vxmpc[j] - 2 * Miu[2] * vrvx - 2 * Miu[2] * I[2] * vrpx / (Miu[2] * med4 + I[2]);
                    vympc[j] = vympc[j] - 2 * Miu[2] * vrvy - 2 * Miu[2] * I[2] * vrpy / (Miu[2] * med4 + I[2]);
                    vzmpc[j] = vzmpc[j] - 2 * Miu[2] * vrvz - 2 * Miu[2] * I[2] * vrpz / (Miu[2] * med4 + I[2]);

                    vcx[2] = vcx[2] + 2 * Miu[2] * vrvx / Mass[2] + 2 * Miu[2] * I[2] * vrpx / ((Miu[2] * med4 + I[2]) * Mass[2]);
                    vcy[2] = vcy[2] + 2 * Miu[2] * vrvy / Mass[2] + 2 * Miu[2] * I[2] * vrpy / ((Miu[2] * med4 + I[2]) * Mass[2]);
                    vcz[2] = vcz[2] + 2 * Miu[2] * vrvz / Mass[2] + 2 * Miu[2] * I[2] * vrpz / ((Miu[2] * med4 + I[2]) * Mass[2]);

                    omegax[2] = omegax[2] + 2 * Miu[2] * (med2 * vrpz - med3 * vrpy) / (Miu[2] * med4 + I[2]);
                    omegay[2] = omegay[2] + 2 * Miu[2] * (med3 * vrpx - med1 * vrpz) / (Miu[2] * med4 + I[2]);
                    omegaz[2] = omegaz[2] + 2 * Miu[2] * (med1 * vrpy - med2 * vrpx) / (Miu[2] * med4 + I[2]);
                }
            }
        }

        // (Optional) Compute energies if needed – omitted or commented out for performance

        // Integrate colloid positions (full step)
        for (k = 1; k <= ncolloid; k++)
        {
            fcx0[k] = fcx[k];
            fcy0[k] = fcy[k];
            fcz0[k] = fcz[k];
            Momx0[k] = Momx[k];
            Momy0[k] = Momy[k];
            Momz0[k] = Momz[k];
            fcx[k] = 0.0;
            fcy[k] = 0.0;
            fcz[k] = 0.0;
            Momx[k] = 0.0;
            Momy[k] = 0.0;
            Momz[k] = 0.0;
            rcx[k] = rcx[k] + dtmd * vcx[k] + dtmd * dtmd * fcx0[k] / (2 * Mass[k]);
            rcy[k] = rcy[k] + dtmd * vcy[k] + dtmd * dtmd * fcy0[k] / (2 * Mass[k]);
            rcz[k] = rcz[k] + dtmd * vcz[k] + dtmd * dtmd * fcz0[k] / (2 * Mass[k]);

            if (rcx[k] <= 0) rcx[k] = xlen + rcx[k];
            else if (rcx[k] > xlen) rcx[k] = rcx[k] - xlen;
            if (rcy[k] <= 0) rcy[k] = ylen + rcy[k];
            else if (rcy[k] > ylen) rcy[k] = rcy[k] - ylen;

            // z 方向使用反射边界，防止雪人球心逃离计算区域
            double lower_z = rcwall[k];
            double upper_z = zlen - rcwall[k];

            if (rcz[k] < lower_z)
            {
                rcz[k] = lower_z;
                vcz[k] = fabs(vcz[k]);
            }
            else if (rcz[k] > upper_z)
            {
                rcz[k] = upper_z;
                vcz[k] = -fabs(vcz[k]);
            }
        }

        // === Check validity of colloid positions after update ===
#if ENABLE_NUMERICAL_CHECKS
        for (k = 1; k <= ncolloid; k++) {
            if (!is_valid_double(rcx[k]) || !is_valid_double(rcy[k]) ||
                !is_valid_double(rcz[k])) {
                fprintf(error_log, "ERROR at step %d, colloid %d: Invalid position\n",
                    i, k);
                fprintf(error_log, "  Position: %lf, %lf, %lf\n",
                    rcx[k], rcy[k], rcz[k]);
                fflush(error_log);
            }
        }
#endif

        // Update colloid orientation (rotation due to angular velocity)
        for (k = 1; k <= 2; k++)
        {
            for (j = 1; j <= 2; j++)
            {
                med1 = nx[k][j];
                med2 = ny[k][j];
                med3 = nz[k][j];

                med4 = omegay[k] * med3 - omegaz[k] * med2;
                med5 = omegaz[k] * med1 - omegax[k] * med3;
                med6 = omegax[k] * med2 - omegay[k] * med1;

                nx[k][j] = med1 + med4 * dtmd;
                ny[k][j] = med2 + med5 * dtmd;
                nz[k][j] = med3 + med6 * dtmd;

                med7 = nx[k][j] * nx[k][j] + ny[k][j] * ny[k][j] + nz[k][j] * nz[k][j];
                med8 = sqrt(med7);

                // Prevent divide-by-zero: check if med8 is too small
                if (med8 < 1e-10) {
                    fprintf(error_log, "WARNING at step %d, colloid %d, dir %d: med8 too small (%le), using default direction\n",
                        i, k, j, med8);
                    fflush(error_log);
                    nx[k][j] = 0.0;
                    ny[k][j] = 1.0;
                    nz[k][j] = 0.0;
                }
                else {
                    nx[k][j] = nx[k][j] / med8;
                    ny[k][j] = ny[k][j] / med8;
                    nz[k][j] = nz[k][j] / med8;
                }
            }
        }

        // Move fluid particles with updated velocities (second half-step of velocity Verlet)
        for (k = 1; k <= ncolloid; k++)
        {
            for (m = 1; m <= a; m++)
            {
                if (mpar[k][m] == 0) break;
                j = mpar[k][m];

                if (record2[j] == 1) continue;
                else if (record2[j] == 0)
                {
                    fx[j] = fxmpc[j];
                    fy[j] = fympc[j];
                    fz[j] = fzmpc[j];
                    fxmpc[j] = 0.0;
                    fympc[j] = 0.0;
                    fzmpc[j] = 0.0;

                    if (rzmpc[j] >= 0.05 && rzmpc[j] <= (zlen - 0.05))
                    {
                        rxmpc[j] = rxmpc[j] + dtmd * vxmpc[j] + dtmd * dtmd * fx[j] / 2;
                        rympc[j] = rympc[j] + dtmd * vympc[j] + dtmd * dtmd * fy[j] / 2;
                        rzmpc[j] = rzmpc[j] + dtmd * vzmpc[j] + dtmd * dtmd * fz[j] / 2;
                    }
                    else if (rzmpc[j] < 0.05)
                    {
                        if (vzmpc[j] >= 0.0)
                        {
                            rxmpc[j] = rxmpc[j] + dtmd * vxmpc[j] + dtmd * dtmd * fx[j] / 2;
                            rympc[j] = rympc[j] + dtmd * vympc[j] + dtmd * dtmd * fy[j] / 2;
                            rzmpc[j] = rzmpc[j] + dtmd * vzmpc[j] + dtmd * dtmd * fz[j] / 2;
                        }
                        else if (vzmpc[j] < 0.0)
                        {
                            // Prevent divide-by-zero: check if vzmpc is too small
                            if (fabs(vzmpc[j]) < 1e-10) {
                                // Use default small time step if velocity is too small
                                dt = dtmd;
                            }
                            else {
                                dt = -rzmpc[j] / vzmpc[j];
                            }
                            if (dt > dtmd)
                            {
                                rxmpc[j] = rxmpc[j] + dtmd * vxmpc[j] + dtmd * dtmd * fx[j] / 2;
                                rympc[j] = rympc[j] + dtmd * vympc[j] + dtmd * dtmd * fy[j] / 2;
                                rzmpc[j] = rzmpc[j] + dtmd * vzmpc[j] + dtmd * dtmd * fz[j] / 2;
                            }
                            else if (dt <= dtmd)
                            {
                                rxmpc[j] = rxmpc[j] + dt * vxmpc[j] + dt * dt * fx[j] / 2;
                                rympc[j] = rympc[j] + dt * vympc[j] + dt * dt * fy[j] / 2;
                                rzmpc[j] = rzmpc[j] + dt * vzmpc[j] + dt * dt * fz[j] / 2;
                                vxmpc[j] = -vxmpc[j];
                                vympc[j] = -vympc[j];
                                vzmpc[j] = -vzmpc[j];
                                rxmpc[j] = rxmpc[j] + (dtmd - dt) * vxmpc[j] + (dtmd - dt) * dtmd * fx[j] / 2;
                                rympc[j] = rympc[j] + (dtmd - dt) * vympc[j] + (dtmd - dt) * dtmd * fy[j] / 2;
                                rzmpc[j] = rzmpc[j] + (dtmd - dt) * vzmpc[j] + (dtmd - dt) * dtmd * fz[j] / 2;
                            }
                        }
                    }
                    else if (rzmpc[j] > (zlen - 0.05))
                    {
                        if (vzmpc[j] <= 0.0)
                        {
                            rxmpc[j] = rxmpc[j] + dtmd * vxmpc[j] + dtmd * dtmd * fx[j] / 2;
                            rympc[j] = rympc[j] + dtmd * vympc[j] + dtmd * dtmd * fy[j] / 2;
                            rzmpc[j] = rzmpc[j] + dtmd * vzmpc[j] + dtmd * dtmd * fz[j] / 2;
                        }
                        else if (vzmpc[j] > 0.0)
                        {
                            // Prevent divide-by-zero: check if vzmpc is too small
                            if (fabs(vzmpc[j]) < 1e-10) {
                                // Use default small time step if velocity is too small
                                dt = dtmd;
                            }
                            else {
                                dt = (zlen - rzmpc[j]) / vzmpc[j];
                            }
                            if (dt > dtmd)
                            {
                                rxmpc[j] = rxmpc[j] + dtmd * vxmpc[j] + dtmd * dtmd * fx[j] / 2;
                                rympc[j] = rympc[j] + dtmd * vympc[j] + dtmd * dtmd * fy[j] / 2;
                                rzmpc[j] = rzmpc[j] + dtmd * vzmpc[j] + dtmd * dtmd * fz[j] / 2;
                            }
                            else if (dt <= dtmd)
                            {
                                rxmpc[j] = rxmpc[j] + dt * vxmpc[j] + dt * dt * fx[j] / 2;
                                rympc[j] = rympc[j] + dt * vympc[j] + dt * dt * fy[j] / 2;
                                rzmpc[j] = rzmpc[j] + dt * vzmpc[j] + dt * dt * fz[j] / 2;
                                vxmpc[j] = -vxmpc[j];
                                vympc[j] = -vympc[j];
                                vzmpc[j] = -vzmpc[j];
                                rxmpc[j] = rxmpc[j] + (dtmd - dt) * vxmpc[j] + (dtmd - dt) * dtmd * fx[j] / 2;
                                rympc[j] = rympc[j] + (dtmd - dt) * vympc[j] + (dtmd - dt) * dtmd * fy[j] / 2;
                                rzmpc[j] = rzmpc[j] + (dtmd - dt) * vzmpc[j] + (dtmd - dt) * dtmd * fz[j] / 2;
                            }
                        }
                    }
                    record2[j] = 1;
                }
            }
        }

        // === 再做一次 x, y 周期边界 ===
        for (j = 1; j <= nparmpc; j++)
        {
            while (rxmpc[j] >= xlen) rxmpc[j] -= xlen;
            while (rxmpc[j] < 0.0)  rxmpc[j] += xlen;

            while (rympc[j] >= ylen) rympc[j] -= ylen;
            while (rympc[j] < 0.0)  rympc[j] += ylen;
        }

        for (k = 1; k <= ncolloid; k++)
        {
            for (m = 1; m <= a; m++)
            {
                if (mpar[k][m] == 0) break;
                j = mpar[k][m];
                record2[j] = 0;
            }
        }

        // (Optional) After-collision energy calculations omitted

        // Update previous forces for next iteration (already done above)

        // Print progress information
        if (i % progress_interval == 0 || i == 1)
        {
            current_time = clock();
            elapsed_seconds = ((double)(current_time - start_time)) / CLOCKS_PER_SEC;

            if (i > 1) {
                time_per_step = elapsed_seconds / i;
                estimated_total_seconds = time_per_step * eqtime;
                remaining_seconds = estimated_total_seconds - elapsed_seconds;
            }
            else {
                time_per_step = 0.0;
                estimated_total_seconds = 0.0;
                remaining_seconds = 0.0;
            }

            double progress_percent = 100.0 * i / eqtime;
            int elapsed_min = (int)(elapsed_seconds / 60);
            int elapsed_sec = (int)(elapsed_seconds) % 60;
            int remaining_min = (int)(remaining_seconds / 60);
            int remaining_sec = (int)(remaining_seconds) % 60;

            printf("Step %d/%d (%.1f%%) | Elapsed: %dm %ds",
                i, eqtime, progress_percent, elapsed_min, elapsed_sec);

            if (i > 1 && remaining_seconds > 0) {
                printf(" | ETA: %dm %ds", remaining_min, remaining_sec);
            }
            if (i > 1 && time_per_step > 0) {
                printf(" | %.3f sec/step", time_per_step);
            }
            printf("\n");
            fflush(stdout);
        }

        // Print positions of colloids at intervals in Tecplot format
        if (i % 20 == 0 || i == 1)
        {
            // Write Zone header
            fprintf(fpw1, " Zone T=\"%d\"\n", i);

            // Write radius value (using sigma[1] as reference, or average)
            fprintf(fpw1, "        %.8f\n", sigma[1]);

            // Write colloid data: x y z Type Radius
            for (k = 1; k <= ncolloid; k++)
            {
                fprintf(fpw1, "       %20.8f %20.8f %20.8f %10.2f %10.2f\n",
                    rcx[k], rcy[k], rcz[k], (double)k, sigma[k]);
            }
            fflush(fpw1);  // Force write to file immediately

            if (i == 1) {
                printf("First step completed! Data written to a1.txt\n");
                fflush(stdout);
            }
        }

        // Thermostat: rescale fluid velocities every 100 steps to maintain temperature
        if (i % 100 == 0)
        {
            double v2sum = 0.0;
            for (int j = 1; j <= nparmpc; j++)
            {
                v2sum += vxmpc[j] * vxmpc[j] + vympc[j] * vympc[j] + vzmpc[j] * vzmpc[j];
            }
            double currentTemp = v2sum / (3.0 * nparmpc);

            // Prevent divide-by-zero: check if currentTemp is too small
            if (currentTemp < 1e-10) {
                fprintf(error_log, "WARNING at step %d: currentTemp too small (%le), skipping temperature rescaling\n",
                    i, currentTemp);
                fflush(error_log);
                continue;
            }

            double scale = sqrt(0.5 / currentTemp);
            for (int j = 1; j <= nparmpc; j++)
            {
                vxmpc[j] = vxmpc[j] * scale;
                vympc[j] = vympc[j] * scale;
                vzmpc[j] = vzmpc[j] * scale;
            }
        }
    }

    // Calculate and print final statistics
    current_time = clock();
    elapsed_seconds = ((double)(current_time - start_time)) / CLOCKS_PER_SEC;
    time_per_step = elapsed_seconds / eqtime;

    int total_hours = (int)(elapsed_seconds / 3600);
    int total_min = (int)(elapsed_seconds / 60) % 60;
    int total_sec = (int)(elapsed_seconds) % 60;

    int total_output_points = (eqtime / 20) + 1;  // +1 for step 1

    printf("\n========================================\n");
    printf("Simulation Completed Successfully!\n");
    printf("========================================\n");
    printf("Statistics:\n");
    printf("  Total steps: %d\n", eqtime);
    printf("  Total runtime: ");
    if (total_hours > 0) {
        printf("%dh %dm %ds", total_hours, total_min, total_sec);
    }
    else if (total_min > 0) {
        printf("%dm %ds", total_min, total_sec);
    }
    else {
        printf("%ds", total_sec);
    }
    printf("\n");
    printf("  Average time per step: %.4f seconds\n", time_per_step);
    printf("  Total output data points: %d\n", total_output_points);
    printf("  Output file: a1.txt\n");
    printf("========================================\n\n");

    fclose(error_log);
    fclose(fpw1);
    return 0;
}
