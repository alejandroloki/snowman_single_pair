
/*****************************/
/*Modified: Single Snowman Model (Active-Inert Pair)*/
/*只保留一对活性-惰性胶体粒子雪人模型*/
/*****************************/

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <malloc.h>
#include <iostream>
#include <fstream>
#define	nxcell 32
#define	nycell 32
#define	nzcell 25
#define a 40000
#define k1 1000
#define k2 0.001
#define k3 0.003
#define ncolloid 2 //total colloid number (只保留一对雪人模型)
#define nparmpc 256000 //total mpc particle number
#define vavempc sqrt(1.5) //the root mean square velocity, kT=0.5, mmpc=1
#define eqtime 200000    //the number of simulation step
using namespace std;

double xlen=1.0*nxcell;   // system size
double ylen=1.0*nycell;
double zlen=1.0*nzcell;
double xhalf=xlen/2;
double yhalf=ylen/2;
double Pi=acos(-1.0);
static double rxmpc[nparmpc+1];  // particle position after collision
static double rympc[nparmpc+1];
static double rzmpc[nparmpc+1];
static double vxmpc[nparmpc+1];  // particle velocity of collision
static double vympc[nparmpc+1];
static double vzmpc[nparmpc+1];
double ang=120*Pi/180;   // rotational angle
double rho=10;
double dis=2.2;
double sinang=sin(ang);
double cosang=cos(ang);
double dtmpc=0.1; // MPC time step
double dtmd=0.002;  //MD time step
static double fxmpc[nparmpc+1]; // interaction in x direction
static double fympc[nparmpc+1]; // interaction in y direction
static double fzmpc[nparmpc+1]; // interaction in z direction
double sigma[ncolloid+1];
double rc[ncolloid+1];
double rn[ncolloid+1];
double square_sigma1,square_sigma2;
double rd1=2.05;
double square_rd1=rd1*rd1;
double rbottem1=2.1;
double square_rbottem1=rbottem1*rbottem1;
double rc_att1=2.2;
double square_rc_att1=rc_att1*rc_att1;
double rd2=2.52;
double square_rd2=rd2*rd2;
double rbottem2=2.56;
double square_rbottem2=rbottem2*rbottem2;
double rc_att2=2.64;
double square_rc_att2=rc_att2*rc_att2;
double rw=25.0;
double square_rw=rw*rw;
double rwall[ncolloid+1];
double rcwall[ncolloid+1];
double square_rwall1,square_rwall2;
double square_rc[ncolloid+1];
double square_rca[ncolloid+1];
double square_rn[ncolloid+1];
double rcx0[ncolloid+1],rcy0[ncolloid+1],rcz0[ncolloid+1];
double rcx1[ncolloid+1],rcy1[ncolloid+1],rcz1[ncolloid+1];
double rcx[ncolloid+1],rcy[ncolloid+1],rcz[ncolloid+1];
double fcx0[ncolloid+1],fcy0[ncolloid+1],fcz0[ncolloid+1];
double fcx[ncolloid+1],fcy[ncolloid+1],fcz[ncolloid+1];
double vcx[ncolloid+1],vcy[ncolloid+1],vcz[ncolloid+1];
double omegax[ncolloid+1],omegay[ncolloid+1],omegaz[ncolloid+1];
double Momx0[ncolloid+1],Momy0[ncolloid+1],Momz0[ncolloid+1];
double Momx[ncolloid+1],Momy[ncolloid+1],Momz[ncolloid+1];
double nx[ncolloid+1][3],ny[ncolloid+1][3],nz[ncolloid+1][3];
double nx_tor[ncolloid+1][3],ny_tor[ncolloid+1][3],nz_tor[ncolloid+1][3];
double torquex[ncolloid+1][3],torquey[ncolloid+1][3],torquez[ncolloid+1][3],torque[ncolloid+1][3];
double costheta[ncolloid+1][3],theta_st[ncolloid+1][3];
double square_rea1;
double square_bb[ncolloid+1];
double Mass[ncolloid+1];
double I[ncolloid+1];
double Miu[ncolloid+1];

/*random number generator*/

#define AMOL .23283064365386962890625e-9
#define NBIT 32
static int ir[256];
static int k;
double r250()
{
        int iran_ks;
        k=(k+1)&255;
        ir[k]=ir[(k-250 & 255)] ^ ir[(k-103 & 255)];
        iran_ks=ir[k];
        return iran_ks*AMOL+0.5;
}
void wmup_ks(int nseed)
{
        int ibm, idum, i;
        double rdum;
    	ibm=2*nseed+1;
        for(k=0; k<256; k++)
        {
            idum=0;
            for(i=1; i<=NBIT; i++)
            {
                idum=idum*2;
                ibm=ibm*16807;
                if(ibm < 0) idum++;
            }
            ir[k]=idum;
        }
        k=0;
        for(i=0;i<100000;i++)
	rdum=r250();
}

//initiate velocity and position of MPC particles
void confmpc()
{
    int i,j;
    double phi,theta;
    double m1=0;
    double m2=0;
    double m3=0;
    double m4;
    double m5;
    double m6;
    double m7;

    // 只初始化粒子1和2(活性-惰性雪人对)
    rcx[1]=16.0;
    rcy[1]=16.0;
    rcz[1]=12.5;

    rcx[2]=16.0-dis;
    rcy[2]=16.0;
    rcz[2]=12.5;

    rcx0[1]=16.0;
    rcy0[1]=16.0;
    rcz0[1]=12.5;

    rcx0[2]=16.0-dis;
    rcy0[2]=16.0;
    rcz0[2]=12.5;

    rcx1[1]=16.0;
    rcy1[1]=16.0;
    rcz1[1]=12.5;

    rcx1[2]=16.0-dis;
    rcy1[2]=16.0;
    rcz1[2]=12.5;

    sigma[1]=2.0;
    rc[1]=sigma[1]*pow(2.0,1.0/12.0);
    rn[1]=rc_att1+0.9;
    square_sigma1=sigma[1]*sigma[1];
    square_rc[1]=rc[1]*rc[1];
    square_rca[1]=square_rc_att1;
    square_rn[1]=rn[1]*rn[1];

    sigma[2]=2.4;
    rc[2]=sigma[2]*pow(2.0,1.0/12.0);
    rn[2]=rc_att2+0.9;
    square_sigma2=sigma[2]*sigma[2];
    square_rc[2]=rc[2]*rc[2];
    square_rca[2]=square_rc_att2;
    square_rn[2]=rn[2]*rn[2];

    square_bb[1]=square_sigma1;
    square_bb[2]=square_sigma2;
    square_rea1=square_rd1;

    for(j=1;j<=ncolloid;j++)
    {
        Mass[j]=4*Pi*sigma[j]*sigma[j]*sigma[j]*rho/3;
    }

    for(j=1;j<=ncolloid;j++)
    {
        I[j]=2*Mass[j]*sigma[j]*sigma[j]/5;
    }

    for(j=1;j<=ncolloid;j++)
    {
        Miu[j]=Mass[j]/(Mass[j]+1);
    }

    for(j=1;j<=ncolloid;j++)
    {
        vcx[j]=0.0;
        vcy[j]=0.0;
        vcz[j]=0.0;
        fcx0[j]=0.0;
        fcy0[j]=0.0;
        fcz0[j]=0.0;
        fcx[j]=0.0;
        fcy[j]=0.0;
        fcz[j]=0.0;
        omegax[j]=0.0;
        omegay[j]=0.0;
        omegaz[j]=0.0;
        Momx0[j]=0.0;
        Momy0[j]=0.0;
        Momz0[j]=0.0;
        Momx[j]=0.0;
        Momy[j]=0.0;
        Momz[j]=0.0;
    }

    nx[1][1]=1.0;
    ny[1][1]=0.0;
    nz[1][1]=0.0;

    nx[2][1]=1.0;
    ny[2][1]=0.0;
    nz[2][1]=0.0;

    nx[1][2]=0.0;
    ny[1][2]=0.0;
    nz[1][2]=1.0;

    nx[2][2]=0.0;
    ny[2][2]=0.0;
    nz[2][2]=1.0;

    rwall[1]=2.25;
    rcwall[1]=rwall[1]*pow(2.0,1.0/48.0);
    square_rwall1=rwall[1]*rwall[1];

    rwall[2]=2.7;
    rcwall[2]=rwall[2]*pow(2.0,1.0/48.0);
    square_rwall2=rwall[2]*rwall[2];

    for(i=1;i<=nparmpc;i++)
    {
loop:   rxmpc[i]=xlen*r250();
        rympc[i]=ylen*r250();
        rzmpc[i]=zlen*r250();

        for(j=1;j<=ncolloid;j++)
        {
            m4=rxmpc[i]-rcx0[j];
            m5=rympc[i]-rcy0[j];
            m6=rzmpc[i]-rcz0[j];

            if(m4>xhalf)m4=m4-xlen;
            else if(m4<-xhalf)m4=m4+xlen;
            if(m5>yhalf)m5=m5-ylen;
            else if(m5<-yhalf)m5=m5+ylen;
            m7=m4*m4+m5*m5+m6*m6;

            if(m7<=square_rca[j])
            {
               goto loop;
            }
         }

         phi=2*Pi*r250();
         theta=Pi*r250();
         vxmpc[i]=vavempc*sin(theta)*cos(phi);
         vympc[i]=vavempc*sin(theta)*sin(phi);
         vzmpc[i]=vavempc*cos(theta);

         m1=m1+vxmpc[i];
         m2=m2+vympc[i];
         m3=m3+vzmpc[i];
         fxmpc[i]=0;
         fympc[i]=0;
         fzmpc[i]=0;
    }

    for(i=1;i<=nparmpc;i++)
    {
        vxmpc[i]=vxmpc[i]-m1/nparmpc;
        vympc[i]=vympc[i]-m2/nparmpc;
        vzmpc[i]=vzmpc[i]-m3/nparmpc;
    }
}

double rpotential1(double r2)
{
   double reppotential1,ir2,ir4,ir8,ir12,ir24;
   ir2=square_sigma1/r2;
   ir4=ir2*ir2;
   ir8=ir4*ir4;
   ir12=ir4*ir8;
   ir24=ir12*ir12;
   reppotential1=4*(ir24-ir12)+1;
   return(reppotential1);
}

double apotential1(double r2)
{
   double attpotential1,r1,r3,ir2,ir4,ir8,ir12,ir24;

   if(r2<=square_rd1)
   {
      ir2=square_sigma1/r2;
      ir4=ir2*ir2;
      ir8=ir4*ir4;
      ir12=ir4*ir8;
      ir24=ir12*ir12;
      attpotential1=4*(ir24-ir12)+1;
      return(attpotential1);
   }

   else if(r2>square_rd1 && r2<=square_rbottem1)
   {
       r1=sqrt(r2);
       r3=r1*r2;
       attpotential1=16404.1823135262*r3-102031.228308921*r2+211503.826889518*r1-146120.452031211;
   }

   else if(r2>square_rbottem1 && r2<=square_rc_att1)
   {
      r1=sqrt(r2);
      r3=r1*r2;
      attpotential1=-2000.0*r3+12900.0*r2-27720.0*r1+19844.0;
   }
   return(attpotential1);
}

double rpotential2(double r2)
{
   double reppotential2,ir2,ir4,ir8,ir12,ir24;
   ir2=square_sigma2/r2;
   ir4=ir2*ir2;
   ir8=ir4*ir4;
   ir12=ir4*ir8;
   ir24=ir12*ir12;
   reppotential2=4*(ir24-ir12)+1;
   return(reppotential2);
}

double apotential2(double r2)
{
   double attpotential2,r1,r3,ir2,ir4,ir8,ir12,ir24;

   if(r2<=square_rd2)
   {
      ir2=square_sigma2/r2;
      ir4=ir2*ir2;
      ir8=ir4*ir4;
      ir12=ir4*ir8;
      ir24=ir12*ir12;
      attpotential2=4*(ir24-ir12)+1;
      return(attpotential2);
   }

   else if(r2>square_rd2 && r2<=square_rbottem2)
   {
       r1=sqrt(r2);
       r3=r1*r2;
       attpotential2=30900.25909844916*r3-235444.903281994*r2+597954.09072102*r1-506172.0754475822;
   }

   else if(r2>square_rbottem2 && r2<=square_rc_att2)
   {
      r1=sqrt(r2);
      r3=r1*r2;
      attpotential2=-3906.25*r3+30468.75*r2-79200.0*r1+68607.0;
   }
   return(attpotential2);
}

double rpotentialwall1(double r2)
{
   double reppotentialwall1,ir2,ir4,ir8,ir12,ir24,ir48,ir96;
   ir2=square_rwall1/r2;
   ir4=ir2*ir2;
   ir8=ir4*ir4;
   ir12=ir4*ir8;
   ir24=ir12*ir12;
   ir48=ir24*ir24;
   ir96=ir48*ir48;
   reppotentialwall1=4*(ir96-ir48)+1;
   return(reppotentialwall1);
}

double rpotentialwall2(double r2)
{
   double reppotentialwall2,ir2,ir4,ir8,ir12,ir24,ir48,ir96;
   ir2=square_rwall2/r2;
   ir4=ir2*ir2;
   ir8=ir4*ir4;
   ir12=ir4*ir8;
   ir24=ir12*ir12;
   ir48=ir24*ir24;
   ir96=ir48*ir48;
   reppotentialwall2=4*(ir96-ir48)+1;
   return(reppotentialwall2);
}

double rforce1(double r2)
{
   double repforce1,ir2,ir4,ir8,ir12,ir24;
   ir2=square_sigma1/r2;
   ir4=ir2*ir2;
   ir8=ir4*ir4;
   ir12=ir4*ir8;
   ir24=ir12*ir12;
   repforce1=(96*ir24-48*ir12)/r2;
   return(repforce1);
}

double aforce1(double r2)
{
   double attforce1,r1,r3,ir2,ir4,ir8,ir12,ir24;

   if(r2<=square_rd1)
   {
      ir2=square_sigma1/r2;
      ir4=ir2*ir2;
      ir8=ir4*ir4;
      ir12=ir4*ir8;
      ir24=ir12*ir12;
      attforce1=(96*ir24-48*ir12)/r2;
      return(attforce1);
   }

   else if(r2>square_rd1 && r2<=square_rbottem1)
   {
      r1=sqrt(r2);
      attforce1=-49212.5469405786*r1+204062.456617842-211503.826889518/r1;
   }

   else if(r2>square_rbottem1 && r2<=square_rc_att1)
   {
      r1=sqrt(r2);
      attforce1=6000.0*r1-25800.0+27720.0/r1;
   }
   return(attforce1);
}

double rforce2(double r2)
{
   double repforce2,ir2,ir4,ir8,ir12,ir24;
   ir2=square_sigma2/r2;
   ir4=ir2*ir2;
   ir8=ir4*ir4;
   ir12=ir4*ir8;
   ir24=ir12*ir12;
   repforce2=(96*ir24-48*ir12)/r2;
   return(repforce2);
}

double aforce2(double r2)
{
   double attforce2,r1,r3,ir2,ir4,ir8,ir12,ir24;

   if(r2<=square_rd2)
   {
      ir2=square_sigma2/r2;
      ir4=ir2*ir2;
      ir8=ir4*ir4;
      ir12=ir4*ir8;
      ir24=ir12*ir12;
      attforce2=(96*ir24-48*ir12)/r2;
      return(attforce2);
   }

   else if(r2>square_rd2 && r2<=square_rbottem2)
   {
      r1=sqrt(r2);
      attforce2=-92700.77729534748*r1+470889.806563988-597954.09072102/r1;
   }

   else if(r2>square_rbottem2 && r2<=square_rc_att2)
   {
      r1=sqrt(r2);
      attforce2=11718.75*r1-60937.5+79200.0/r1;
   }
   return(attforce2);
}

double rforcewall1(double r2)
{
   double repforcewall1,ir2,ir4,ir8,ir12,ir24,ir48,ir96;
   ir2=square_rwall1/r2;
   ir4=ir2*ir2;
   ir8=ir4*ir4;
   ir12=ir4*ir8;
   ir24=ir12*ir12;
   ir48=ir24*ir24;
   ir96=ir48*ir48;
   repforcewall1=(384*ir96-192*ir48)/r2;
   return(repforcewall1);
}

double rforcewall2(double r2)
{
   double repforcewall2,ir2,ir4,ir8,ir12,ir24,ir48,ir96;
   ir2=square_rwall2/r2;
   ir4=ir2*ir2;
   ir8=ir4*ir4;
   ir12=ir4*ir8;
   ir24=ir12*ir12;
   ir48=ir24*ir24;
   ir96=ir48*ir48;
   repforcewall2=(384*ir96-192*ir48)/r2;
   return(repforcewall2);
}

/*main function*/
int main()
{
    int nseed;                      // random number seed
    nseed= (unsigned)time(NULL);
    wmup_ks(nseed);
    FILE *fpw1;
    confmpc();
    int i,j,k,l,m;
    int inside[ncolloid+1];
    double record1[ncolloid+1];
    static int record2[nparmpc+1];
    static int record3[nparmpc+1];
    static int kind[nparmpc+1];
    double bx=xlen/nxcell;
    double by=ylen/nycell;
    double bz=zlen/nzcell;
    double p1=1.0;
    double p2=1.0/1000.0;
    double cosphi;
    static int countcellmpc[nxcell*nycell*(nzcell+2)+1];  // number of particle in each collision box
    int me1,me2,me3,me4; //medium variable
    double med1,med2,med3,med4,med5,med6,med7,med8,med9,med10;
    static double vxmc[nxcell*nycell*(nzcell+2)+1];  // center of mass velocity
    static double vymc[nxcell*nycell*(nzcell+2)+1];
    static double vzmc[nxcell*nycell*(nzcell+2)+1];
    static double dirx[nxcell*nycell*(nzcell+2)+1];  // rotational axis
    static double diry[nxcell*nycell*(nzcell+2)+1];
    static double dirz[nxcell*nycell*(nzcell+2)+1];
    double vx,vy,vz; //relative velocity
    double vpx,vpy,vpz; //relative velocity parallel to random vector
    double vvx,vvy,vvz; //relative velocity vertical to random vector
    static int cellmpc[nparmpc+1]; // the collision box that a MPC particle belongs to
    double xshift,yshift,zshift;
    static int npar[nparmpc+1];
    static int mpar[ncolloid+1][a+1];
    double E,E1,E2,E3,E4,E5;
    static double fx[nparmpc+1];
    static double fy[nparmpc+1];
    static double fz[nparmpc+1];
    double nx_st,ny_st,nz_st;
    double vrx,vry,vrz;
    double vrvx,vrvy,vrvz;
    double vrpx,vrpy,vrpz;
    double dt;

    fpw1=fopen("a1_snowman.txt","w");
    if(fpw1==NULL)
    {
       printf("Cannot open this file!! " );
       exit(0);
    }

    for(i=1;i<=nxcell*nycell*(nzcell+2);i++)
    {
	    countcellmpc[i]=0;
    }

    for(j=1;j<=nparmpc;j++)
    {
        kind[j]=1;
    }

	//srd simulation
    for(i=1;i<=eqtime;i++)
    {
	    for(j=1;j<=nparmpc;j++)
        {
            npar[j]=0;
            record2[j]=0;
            record3[j]=0;
        }

	    for(j=1;j<=ncolloid;j++)
        {
            inside[j]=0;

            for(m=1;m<=a;m++)
            {
                mpar[j][m]=0;
            }
        }

	    for(k=1;k<=ncolloid;k++)
        {
            for(j=1;j<=nparmpc;j++)
            {
                med1=rxmpc[j]-rcx[k];
                med2=rympc[j]-rcy[k];
                med3=rzmpc[j]-rcz[k];

                if(med1>xhalf)med1=med1-xlen;
                else if(med1<-xhalf)med1=med1+xlen;
                if(med2>yhalf)med2=med2-ylen;
                else if(med2<-yhalf)med2=med2+ylen;
                med4=med1*med1+med2*med2+med3*med3;

                if(med4<=square_rn[k])
                {
                   npar[j]=1;
                   inside[k]=inside[k]+1;
                   mpar[k][inside[k]]=j;
                }
            }
        }

	    for(j=1;j<=nparmpc;j++)
	    {
            if(npar[j]==0)
            {
               if(rzmpc[j]>=0.9 && rzmpc[j]<=(zlen-0.9))
               {
                  rxmpc[j]=rxmpc[j]+vxmpc[j]*dtmpc;
                  rympc[j]=rympc[j]+vympc[j]*dtmpc;
                  rzmpc[j]=rzmpc[j]+vzmpc[j]*dtmpc;
               }

               else if(rzmpc[j]<0.9)
               {
                  if(vzmpc[j]>=0.0)
                  {
                     rxmpc[j]=rxmpc[j]+vxmpc[j]*dtmpc;
                     rympc[j]=rympc[j]+vympc[j]*dtmpc;
                     rzmpc[j]=rzmpc[j]+vzmpc[j]*dtmpc;
                  }
                  else if(vzmpc[j]<0.0)
                  {
                     dt=-rzmpc[j]/vzmpc[j];
                     if(dt>dtmpc)
                     {
                        rxmpc[j]=rxmpc[j]+vxmpc[j]*dtmpc;
                        rympc[j]=rympc[j]+vympc[j]*dtmpc;
                        rzmpc[j]=rzmpc[j]+vzmpc[j]*dtmpc;
                     }
                     else if(dt<=dtmpc)
                     {
                        rxmpc[j]=rxmpc[j]+vxmpc[j]*dt;
                        rympc[j]=rympc[j]+vympc[j]*dt;
                        rzmpc[j]=rzmpc[j]+vzmpc[j]*dt;
                        vxmpc[j]=-vxmpc[j];
                        vympc[j]=-vympc[j];
                        vzmpc[j]=-vzmpc[j];
                        rxmpc[j]=rxmpc[j]+vxmpc[j]*(dtmpc-dt);
                        rympc[j]=rympc[j]+vympc[j]*(dtmpc-dt);
                        rzmpc[j]=rzmpc[j]+vzmpc[j]*(dtmpc-dt);
                     }
                  }
               }

               else if(rzmpc[j]>(zlen-0.9))
               {
                  if(vzmpc[j]<=0.0)
                  {
                     rxmpc[j]=rxmpc[j]+vxmpc[j]*dtmpc;
                     rympc[j]=rympc[j]+vympc[j]*dtmpc;
                     rzmpc[j]=rzmpc[j]+vzmpc[j]*dtmpc;
                  }
                  else if(vzmpc[j]>0.0)
                  {
                     dt=(zlen-rzmpc[j])/vzmpc[j];
                     if(dt>dtmpc)
                     {
                        rxmpc[j]=rxmpc[j]+vxmpc[j]*dtmpc;
                        rympc[j]=rympc[j]+vympc[j]*dtmpc;
                        rzmpc[j]=rzmpc[j]+vzmpc[j]*dtmpc;
                     }
                     else if(dt<=dtmpc)
                     {
                        rxmpc[j]=rxmpc[j]+vxmpc[j]*dt;
                        rympc[j]=rympc[j]+vympc[j]*dt;
                        rzmpc[j]=rzmpc[j]+vzmpc[j]*dt;
                        vxmpc[j]=-vxmpc[j];
                        vympc[j]=-vympc[j];
                        vzmpc[j]=-vzmpc[j];
                        rxmpc[j]=rxmpc[j]+vxmpc[j]*(dtmpc-dt);
                        rympc[j]=rympc[j]+vympc[j]*(dtmpc-dt);
                        rzmpc[j]=rzmpc[j]+vzmpc[j]*(dtmpc-dt);
                     }
                  }
               }

               if(rxmpc[j]<=0)rxmpc[j]=xlen+rxmpc[j];                   //periodic boundary
               else if(rxmpc[j]>xlen)rxmpc[j]=rxmpc[j]-xlen;
               if(rympc[j]<=0)rympc[j]=ylen+rympc[j];
               else if(rympc[j]>ylen)rympc[j]=rympc[j]-ylen;

               if(kind[j]==2)
               {
                  for(k=1;k<=ncolloid;k++)
                  {
                      med1=rxmpc[j]-rcx[k];
                      med2=rympc[j]-rcy[k];
                      med3=rzmpc[j]-rcz[k];

                      if(med1>xhalf)med1=med1-xlen;
                      else if(med1<-xhalf)med1=med1+xlen;
                      if(med2>yhalf)med2=med2-ylen;
                      else if(med2<-yhalf)med2=med2+ylen;
                      med4=med1*med1+med2*med2+med3*med3;
                      record1[k]=med4;
                  }

                  for(k=1;k<=ncolloid;k++)
                  {
                      if(record1[k]<square_rw)
                      {
                         kind[j]=3;
                      }
                  }
               }
            }
	    }

	    for(j=1;j<=nparmpc;j++)
        {
            if(npar[j]==0)
            {
               if(kind[j]==2)
               {
                  if(r250()<=p2)
                  {
                     kind[j]=1;
                  }
               }

               else if(kind[j]==3)
               {
                  kind[j]=2;
               }
            }
        }

        for(l=1;l<=50;l++)
        {

            for(k=1;k<=ncolloid;k++)
            {
                fcx0[k]=fcx[k];
                fcy0[k]=fcy[k];
                fcz0[k]=fcz[k];
                Momx0[k]=Momx[k];
                Momy0[k]=Momy[k];
                Momz0[k]=Momz[k];
                fcx[k]=0.0;
                fcy[k]=0.0;
                fcz[k]=0.0;
                Momx[k]=0.0;
                Momy[k]=0.0;
                Momz[k]=0.0;
                rcx[k]=rcx[k]+dtmd*vcx[k]+dtmd*dtmd*fcx0[k]/(2*Mass[k]);
                rcy[k]=rcy[k]+dtmd*vcy[k]+dtmd*dtmd*fcy0[k]/(2*Mass[k]);
                rcz[k]=rcz[k]+dtmd*vcz[k]+dtmd*dtmd*fcz0[k]/(2*Mass[k]);
                rcx1[k]=rcx1[k]+dtmd*vcx[k]+dtmd*dtmd*fcx0[k]/(2*Mass[k]);
                rcy1[k]=rcy1[k]+dtmd*vcy[k]+dtmd*dtmd*fcy0[k]/(2*Mass[k]);
                rcz1[k]=rcz1[k]+dtmd*vcz[k]+dtmd*dtmd*fcz0[k]/(2*Mass[k]);

                if(rcx[k]<=0)rcx[k]=xlen+rcx[k];
                else if(rcx[k]>xlen)rcx[k]=rcx[k]-xlen;
                if(rcy[k]<=0)rcy[k]=ylen+rcy[k];
                else if(rcy[k]>ylen)rcy[k]=rcy[k]-ylen;
            }

            // 更新粒子1和2的方向向量
            for(k=1;k<=2;k++)
            {
                for(j=1;j<=2;j++)
                {
                    med1=nx[k][j];
                    med2=ny[k][j];
                    med3=nz[k][j];

                    med4=omegay[k]*med3-omegaz[k]*med2;
                    med5=omegaz[k]*med1-omegax[k]*med3;
                    med6=omegax[k]*med2-omegay[k]*med1;

                    nx[k][j]=med1+med4*dtmd;
                    ny[k][j]=med2+med5*dtmd;
                    nz[k][j]=med3+med6*dtmd;

                    med7=nx[k][j]*nx[k][j]+ny[k][j]*ny[k][j]+nz[k][j]*nz[k][j];
                    med8=sqrt(med7);
                    nx[k][j]=nx[k][j]/med8;
                    ny[k][j]=ny[k][j]/med8;
                    nz[k][j]=nz[k][j]/med8;
                }
            }

            for(k=1;k<=ncolloid;k++)
            {
                for(m=1;m<=a;m++)
                {
                    if(mpar[k][m]==0)break;
                    j=mpar[k][m];

                    if(record2[j]==1)continue;

                    else if(record2[j]==0)
                    {
                       fx[j]=fxmpc[j];
                       fy[j]=fympc[j];
                       fz[j]=fzmpc[j];
                       fxmpc[j]=0.0;
                       fympc[j]=0.0;
                       fzmpc[j]=0.0;

                       if(rzmpc[j]>=0.05 && rzmpc[j]<=(zlen-0.05))
                       {
                          rxmpc[j]=rxmpc[j]+dtmd*vxmpc[j]+dtmd*dtmd*fx[j]/2;
                          rympc[j]=rympc[j]+dtmd*vympc[j]+dtmd*dtmd*fy[j]/2;
                          rzmpc[j]=rzmpc[j]+dtmd*vzmpc[j]+dtmd*dtmd*fz[j]/2;
                       }

                       else if(rzmpc[j]<0.05)
                       {
                          if(vzmpc[j]>=0.0)
                          {
                             rxmpc[j]=rxmpc[j]+dtmd*vxmpc[j]+dtmd*dtmd*fx[j]/2;
                             rympc[j]=rympc[j]+dtmd*vympc[j]+dtmd*dtmd*fy[j]/2;
                             rzmpc[j]=rzmpc[j]+dtmd*vzmpc[j]+dtmd*dtmd*fz[j]/2;
                          }
                          else if(vzmpc[j]<0.0)
                          {
                             dt=-rzmpc[j]/vzmpc[j];
                             if(dt>dtmd)
                             {
                                rxmpc[j]=rxmpc[j]+dtmd*vxmpc[j]+dtmd*dtmd*fx[j]/2;
                                rympc[j]=rympc[j]+dtmd*vympc[j]+dtmd*dtmd*fy[j]/2;
                                rzmpc[j]=rzmpc[j]+dtmd*vzmpc[j]+dtmd*dtmd*fz[j]/2;
                             }
                             else if(dt<=dtmd)
                             {
                                rxmpc[j]=rxmpc[j]+dt*vxmpc[j]+dt*dt*fx[j]/2;
                                rympc[j]=rympc[j]+dt*vympc[j]+dt*dt*fy[j]/2;
                                rzmpc[j]=rzmpc[j]+dt*vzmpc[j]+dt*dt*fz[j]/2;
                                vxmpc[j]=-vxmpc[j];
                                vympc[j]=-vympc[j];
                                vzmpc[j]=-vzmpc[j];
                                rxmpc[j]=rxmpc[j]+(dtmd-dt)*vxmpc[j]+(dtmd-dt)*(dtmd-dt)*fx[j]/2;
                                rympc[j]=rympc[j]+(dtmd-dt)*vympc[j]+(dtmd-dt)*(dtmd-dt)*fy[j]/2;
                                rzmpc[j]=rzmpc[j]+(dtmd-dt)*vzmpc[j]+(dtmd-dt)*(dtmd-dt)*fz[j]/2;
                             }
                          }
                       }

                       else if(rzmpc[j]>(zlen-0.05))
                       {
                          if(vzmpc[j]<=0.0)
                          {
                             rxmpc[j]=rxmpc[j]+dtmd*vxmpc[j]+dtmd*dtmd*fx[j]/2;
                             rympc[j]=rympc[j]+dtmd*vympc[j]+dtmd*dtmd*fy[j]/2;
                             rzmpc[j]=rzmpc[j]+dtmd*vzmpc[j]+dtmd*dtmd*fz[j]/2;
                          }
                          else if(vzmpc[j]>0.0)
                          {
                             dt=(zlen-rzmpc[j])/vzmpc[j];
                             if(dt>dtmd)
                             {
                                rxmpc[j]=rxmpc[j]+dtmd*vxmpc[j]+dtmd*dtmd*fx[j]/2;
                                rympc[j]=rympc[j]+dtmd*vympc[j]+dtmd*dtmd*fy[j]/2;
                                rzmpc[j]=rzmpc[j]+dtmd*vzmpc[j]+dtmd*dtmd*fz[j]/2;
                             }
                             else if(dt<=dtmd)
                             {
                                rxmpc[j]=rxmpc[j]+dt*vxmpc[j]+dt*dt*fx[j]/2;
                                rympc[j]=rympc[j]+dt*vympc[j]+dt*dt*fy[j]/2;
                                rzmpc[j]=rzmpc[j]+dt*vzmpc[j]+dt*dt*fz[j]/2;
                                vxmpc[j]=-vxmpc[j];
                                vympc[j]=-vympc[j];
                                vzmpc[j]=-vzmpc[j];
                                rxmpc[j]=rxmpc[j]+(dtmd-dt)*vxmpc[j]+(dtmd-dt)*(dtmd-dt)*fx[j]/2;
                                rympc[j]=rympc[j]+(dtmd-dt)*vympc[j]+(dtmd-dt)*(dtmd-dt)*fy[j]/2;
                                rzmpc[j]=rzmpc[j]+(dtmd-dt)*vzmpc[j]+(dtmd-dt)*(dtmd-dt)*fz[j]/2;
                             }
                          }
                       }

                       if(rxmpc[j]<=0)rxmpc[j]=xlen+rxmpc[j];              //periodic boundary
                       else if(rxmpc[j]>xlen)rxmpc[j]=rxmpc[j]-xlen;
                       if(rympc[j]<=0)rympc[j]=ylen+rympc[j];
                       else if(rympc[j]>ylen)rympc[j]=rympc[j]-ylen;
                       record2[j]=1;
                    }
                }
            }

            for(k=1;k<=ncolloid;k++)
            {
                for(m=1;m<=a;m++)
                {
                    if(mpar[k][m]==0)break;
                    j=mpar[k][m];
                    record2[j]=0;
                }
            }

            //chemical reaction - 只针对粒子1
            for(m=1;m<=a;m++)
            {
                if(mpar[1][m]==0)break;
                j=mpar[1][m];

                med1=rxmpc[j]-rcx[1];
                med2=rympc[j]-rcy[1];
                med3=rzmpc[j]-rcz[1];

                if(med1>xhalf)med1=med1-xlen;
                else if(med1<-xhalf)med1=med1+xlen;
                if(med2>yhalf)med2=med2-ylen;
                else if(med2<-yhalf)med2=med2+ylen;
                med4=med1*med1+med2*med2+med3*med3;

                if(med4<=square_rea1)
                {
                   if(kind[j]==1)
                   {
                      if(r250()<=p1)
                      {
                         record3[j]=1;
                      }
                   }
                }

                med5=rxmpc[j]-rcx[2];
                med6=rympc[j]-rcy[2];
                med7=rzmpc[j]-rcz[2];
                if(med5>xhalf)med5=med5-xlen;
                else if(med5<-xhalf)med5=med5+xlen;
                if(med6>yhalf)med6=med6-ylen;
                else if(med6<-yhalf)med6=med6+ylen;
                med8=med5*med5+med6*med6+med7*med7;

                if(med8>square_rc_att2)
                {
                   if(record3[j]==1)
                   {
                      kind[j]=2;
                      record3[j]=0;
                   }
                }
            }

            // calculate new force - 粒子1与MPC粒子的相互作用
            for(m=1;m<=a;m++)
            {
                if(mpar[1][m]==0)break;
                j=mpar[1][m];

                med1=rxmpc[j]-rcx[1];
                med2=rympc[j]-rcy[1];
                med3=rzmpc[j]-rcz[1];

                if(med1>xhalf)med1=med1-xlen;
                else if(med1<-xhalf)med1=med1+xlen;
                if(med2>yhalf)med2=med2-ylen;
                else if(med2<-yhalf)med2=med2+ylen;
                med4=med1*med1+med2*med2+med3*med3;

                if(kind[j]==1)
                {
                   if(med4<=square_rc[1])
                   {
                      med5=rforce1(med4);
                      med6=med1*med5;
                      med7=med2*med5;
                      med8=med3*med5;
                      fxmpc[j]=fxmpc[j]+med6;
                      fympc[j]=fympc[j]+med7;
                      fzmpc[j]=fzmpc[j]+med8;
                      fcx[1]=fcx[1]-med6;
                      fcy[1]=fcy[1]-med7;
                      fcz[1]=fcz[1]-med8;
                   }
                }

                else if(kind[j]==2)
                {
                   if(med4<=square_rca[1])
                   {
                      med5=aforce1(med4);
                      med6=med1*med5;
                      med7=med2*med5;
                      med8=med3*med5;
                      fxmpc[j]=fxmpc[j]+med6;
                      fympc[j]=fympc[j]+med7;
                      fzmpc[j]=fzmpc[j]+med8;
                      fcx[1]=fcx[1]-med6;
                      fcy[1]=fcy[1]-med7;
                      fcz[1]=fcz[1]-med8;
                   }
                }
            }

            // 粒子2与MPC粒子的相互作用
            for(m=1;m<=a;m++)
            {
                if(mpar[2][m]==0)break;
                j=mpar[2][m];

                med1=rxmpc[j]-rcx[2];
                med2=rympc[j]-rcy[2];
                med3=rzmpc[j]-rcz[2];

                if(med1>xhalf)med1=med1-xlen;
                else if(med1<-xhalf)med1=med1+xlen;
                if(med2>yhalf)med2=med2-ylen;
                else if(med2<-yhalf)med2=med2+ylen;
                med4=med1*med1+med2*med2+med3*med3;

                if(kind[j]==1)
                {
                   if(med4<=square_rc[2])
                   {
                      med5=rforce2(med4);
                      med6=med1*med5;
                      med7=med2*med5;
                      med8=med3*med5;
                      fxmpc[j]=fxmpc[j]+med6;
                      fympc[j]=fympc[j]+med7;
                      fzmpc[j]=fzmpc[j]+med8;
                      fcx[2]=fcx[2]-med6;
                      fcy[2]=fcy[2]-med7;
                      fcz[2]=fcz[2]-med8;
                   }
                }

                else if(kind[j]==2)
                {
                   if(med4<=square_rca[2])
                   {
                      med5=aforce2(med4);
                      med6=med1*med5;
                      med7=med2*med5;
                      med8=med3*med5;
                      fxmpc[j]=fxmpc[j]+med6;
                      fympc[j]=fympc[j]+med7;
                      fzmpc[j]=fzmpc[j]+med8;
                      fcx[2]=fcx[2]-med6;
                      fcy[2]=fcy[2]-med7;
                      fcz[2]=fcz[2]-med8;
                   }
                }
            }

            for(k=1;k<=ncolloid;k++)
            {
                for(m=1;m<=a;m++)
                {
                    if(mpar[k][m]==0)break;
                    j=mpar[k][m];

                    if(record2[j]==1)continue;

                    else if(record2[j]==0)
                    {
                       vxmpc[j]=vxmpc[j]+dtmd*(fxmpc[j]+fx[j])/2;
                       vympc[j]=vympc[j]+dtmd*(fympc[j]+fy[j])/2;
                       vzmpc[j]=vzmpc[j]+dtmd*(fzmpc[j]+fz[j])/2;
                       record2[j]=1;
                    }
                }
            }

            for(k=1;k<=ncolloid;k++)
            {
                for(m=1;m<=a;m++)
                {
                    if(mpar[k][m]==0)break;
                    j=mpar[k][m];
                    record2[j]=0;
                }
            }

            // 粒子1-2之间的约束力
            med1=rcx[1]-rcx[2];
            med2=rcy[1]-rcy[2];
            med3=rcz[1]-rcz[2];

            if(med1>xhalf)med1=med1-xlen;
            else if(med1<-xhalf)med1=med1+xlen;
            if(med2>yhalf)med2=med2-ylen;
            else if(med2<-yhalf)med2=med2+ylen;
            med4=med1*med1+med2*med2+med3*med3;
            med5=sqrt(med4);

            nx_st=med1/med5;
            ny_st=med2/med5;
            nz_st=med3/med5;

            fcx[1]=fcx[1]-med1*k1*(med5-dis)/med5;
            fcy[1]=fcy[1]-med2*k1*(med5-dis)/med5;
            fcz[1]=fcz[1]-med3*k1*(med5-dis)/med5;
            fcx[2]=fcx[2]+med1*k1*(med5-dis)/med5;
            fcy[2]=fcy[2]+med2*k1*(med5-dis)/med5;
            fcz[2]=fcz[2]+med3*k1*(med5-dis)/med5;

            for(k=1;k<=2;k++)
            {
                costheta[k][1]=nx[k][1]*nx_st+ny[k][1]*ny_st+nz[k][1]*nz_st;
                theta_st[k][1]=acos(costheta[k][1]);
                torque[k][1]=k1*theta_st[k][1];

                nx_tor[k][1]=ny[k][1]*nz_st-nz[k][1]*ny_st;
                ny_tor[k][1]=nz[k][1]*nx_st-nx[k][1]*nz_st;
                nz_tor[k][1]=nx[k][1]*ny_st-ny[k][1]*nx_st;

                torquex[k][1]=torque[k][1]*nx_tor[k][1];
                torquey[k][1]=torque[k][1]*ny_tor[k][1];
                torquez[k][1]=torque[k][1]*nz_tor[k][1];

                omegax[k]=omegax[k]+dtmd*torquex[k][1]/(2*I[k]);
                omegay[k]=omegay[k]+dtmd*torquey[k][1]/(2*I[k]);
                omegaz[k]=omegaz[k]+dtmd*torquez[k][1]/(2*I[k]);
            }

            costheta[1][2]=nx[1][2]*nx[2][2]+ny[1][2]*ny[2][2]+nz[1][2]*nz[2][2];
            theta_st[1][2]=acos(costheta[1][2]);
            torque[1][2]=k1*theta_st[1][2];

            nx_tor[1][2]=ny[1][2]*nz[2][2]-nz[1][2]*ny[2][2];
            ny_tor[1][2]=nz[1][2]*nx[2][2]-nx[1][2]*nz[2][2];
            nz_tor[1][2]=nx[1][2]*ny[2][2]-ny[1][2]*nx[2][2];

            torquex[1][2]=torque[1][2]*nx_tor[1][2];
            torquey[1][2]=torque[1][2]*ny_tor[1][2];
            torquez[1][2]=torque[1][2]*nz_tor[1][2];
            torquex[2][2]=-torque[1][2]*nx_tor[1][2];
            torquey[2][2]=-torque[1][2]*ny_tor[1][2];
            torquez[2][2]=-torque[1][2]*nz_tor[1][2];

            for(k=1;k<=2;k++)
            {
                omegax[k]=omegax[k]+dtmd*torquex[k][2]/(2*I[k]);
                omegay[k]=omegay[k]+dtmd*torquey[k][2]/(2*I[k]);
                omegaz[k]=omegaz[k]+dtmd*torquez[k][2]/(2*I[k]);
            }

            // 墙壁作用力
            if(rcz[1]<rcwall[1])
            {
               fcz[1]=fcz[1]+rcz[1]*rforcewall1(rcz[1]*rcz[1]);
            }

            else if(rcz[1]>=(zlen-rcwall[1]))
            {
               fcz[1]=fcz[1]-(zlen-rcz[1])*rforcewall1((zlen-rcz[1])*(zlen-rcz[1]));
            }

            if(rcz[2]<rcwall[2])
            {
               fcz[2]=fcz[2]+rcz[2]*rforcewall2(rcz[2]*rcz[2]);
            }

            else if(rcz[2]>=(zlen-rcwall[2]))
            {
               fcz[2]=fcz[2]-(zlen-rcz[2])*rforcewall2((zlen-rcz[2])*(zlen-rcz[2]));
            }

            // 重力
            fcz[1]=fcz[1]-k2*Mass[1];
            fcz[2]=fcz[2]-k3*Mass[2];

            for(k=1;k<=ncolloid;k++)
            {
                vcx[k]=vcx[k]+dtmd*(fcx0[k]+fcx[k])/(2*Mass[k]);
                vcy[k]=vcy[k]+dtmd*(fcy0[k]+fcy[k])/(2*Mass[k]);
                vcz[k]=vcz[k]+dtmd*(fcz0[k]+fcz[k])/(2*Mass[k]);
            }

        }

        xshift=(r250()-0.5);
        yshift=(r250()-0.5);
        zshift=(r250()-0.5);

        for(i=1;i<=nxcell*nycell*(nzcell+2);i++)
        {
            countcellmpc[i]=0;
            vxmc[i]=0.0;
            vymc[i]=0.0;
            vzmc[i]=0.0;
        }

        for(j=1;j<=nparmpc;j++)
        {
            if(npar[j]==0)
            {
               me1=(int)((rxmpc[j]+xshift)/bx);
               me2=(int)((rympc[j]+yshift)/by);
               me3=(int)((rzmpc[j]+zshift)/bz);

               if(me1<0)me1=nxcell+me1;
               else if(me1>=nxcell)me1=me1-nxcell;
               if(me2<0)me2=nycell+me2;
               else if(me2>=nycell)me2=me2-nycell;

               me4=me1*nycell*(nzcell+2)+me2*(nzcell+2)+me3+1;
               cellmpc[j]=me4;
               countcellmpc[me4]=countcellmpc[me4]+1;
            }
        }

        for(j=1;j<=nparmpc;j++)
        {
            if(npar[j]==0)
            {
               vxmc[cellmpc[j]]=vxmc[cellmpc[j]]+vxmpc[j];
               vymc[cellmpc[j]]=vymc[cellmpc[j]]+vympc[j];
               vzmc[cellmpc[j]]=vzmc[cellmpc[j]]+vzmpc[j];
            }
        }

        for(i=1;i<=nxcell*nycell*(nzcell+2);i++)
        {
            if(countcellmpc[i]!=0)
            {
               vxmc[i]=vxmc[i]/countcellmpc[i];
               vymc[i]=vymc[i]/countcellmpc[i];
               vzmc[i]=vzmc[i]/countcellmpc[i];
            }
        }

        for(i=1;i<=nxcell*nycell*(nzcell+2);i++)
        {
            if(countcellmpc[i]!=0)
            {
               cosphi=1.0-2*r250();
               dirx[i]=cosphi;
               diry[i]=sqrt(1.0-cosphi*cosphi)*cos(2*Pi*r250());
               dirz[i]=sqrt(1.0-cosphi*cosphi)*sin(2*Pi*r250());
            }
        }

        for(j=1;j<=nparmpc;j++)
        {
            if(npar[j]==0)
            {
            vx=vxmpc[j]-vxmc[cellmpc[j]];
            vy=vympc[j]-vymc[cellmpc[j]];
            vz=vzmpc[j]-vzmc[cellmpc[j]];
            vpx=vx*dirx[cellmpc[j]]+vy*diry[cellmpc[j]]+vz*dirz[cellmpc[j]];
            vpy=vpx*dirx[cellmpc[j]];
            vpz=vpx*dirz[cellmpc[j]];
            vpy=vpx*diry[cellmpc[j]];
            vvx=vx-vpy;
            vvy=vy-vpy;
            vvz=vz-vpz;
            vxmpc[j]=vxmc[cellmpc[j]]+vpx+vvx*cosang+(vvy*dirz[cellmpc[j]]-vvz*diry[cellmpc[j]])*sinang;
            vympc[j]=vymc[cellmpc[j]]+vpy+vvy*cosang+(vvz*dirx[cellmpc[j]]-vvx*dirz[cellmpc[j]])*sinang;
            vzmpc[j]=vzmc[cellmpc[j]]+vpz+vvz*cosang+(vvx*diry[cellmpc[j]]-vvy*dirx[cellmpc[j]])*sinang;
            }
        }


        if(i%10000==0)
        {
           cout<<i<<endl;
        }

        if(i%20==0)
        {
           for(k=1;k<=ncolloid;k++)
           {
               fprintf(fpw1,"%d    %lf    %lf    %lf\n",i,rcx1[k],rcy1[k],rcz1[k]);
           }
        }
    }

    fclose(fpw1);
    return 0;
}
