#include <iostream>
#include <fstream>
#include <math.h> 
#include <utility>
#include <boost/math/tools/roots.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
using namespace std;
#define thread_num 1
#define hb (0.658*pow(10,-21))
#define v (3.0*pow(10,23))
#define qq (1.0/3.0*0.029*pow(10,-15))
#define m (1.66*pow(10,-45))
#define BC (1.17912*pow(10,19))
#define la (1.312)
#define ronum 170 //number of points in rho
#define bnum 100 //number of points in magnetic field
#define gnum 5  //number of points in polarization
long double jmax(double ef,double B)
{
    return BC*(ef*ef-1)/(2.0*B);
}
double hkf(double ef,double B)
{
    long double jj=floor(jmax(ef,B));
    double wyn=0.0;
    double gg,Z;
    for (int a=0;a<=jj;a++)
    {
        if (a==0)
        {
            gg=1;
        }
        else
        {
            gg=3;
        }
        Z=gg*sqrt(ef*ef-1-2.0*a*B/BC);
        //cout<<B/BC<<endl;
        if (isnan(Z))
        {
            wyn+=0;
        }
        else
        {
            wyn+=Z;
        }
    }
    return wyn;
}
double roo(double ef, double B)
{
    double wyn=(2.0*B)/(BC*4.0*M_PI*M_PI*la*la*la)*hkf(ef,B);
    return wyn;
}
void gen_tab(double *tab,double start, double end, int n)
{
    double h=(end-start)/n;
    for (int i=0;i<n;i++)
    {
        *tab=start+h*i;
        tab++;
    }
}

void gen_tab_cheb(double *tab,double start, double end, int n)
{
    for (int i=0;i<n;i++)
    {
        *tab=0.5*(start+end)+0.5*(end-start)*cos(M_PI*(1-(2*i+1)/(2 * (double) n)));
        tab++;
    }
}

class funcion
{
    private:
        double ro,B;
    public:
        funcion(double a,double c)
        {
            ro=a;
            B=c;
        }
        double operator()(double x)
        {
            return ro-roo(x,B);
        }
};

double find_ef(double ro,double B)
{
    double a=1;
    double b=100;
    double fa;
    double fb;
    fa=(ro-roo(a,B));
    fb=(ro-roo(b,B));
    while (fa*fb>0)
    {
        a=b;
        fa=fb;
        b*=2;
        fb=(ro-roo(b,B));
    }
    double x=(a+b)/2;
    while (true)
    {
        double c=roo(x,B)-ro;
        if (c<0)
        {
            a=x;
        }
        else
        {
            b=x;
        }
        if (b-a<pow(10,-1) && abs(c)<0.1)
        {
            return x;
        }
        x=(a+b)/2;
    }    
}
class war
{
    public:
        double B;
        double ro;
        war(double a,double b)
        {
            B=a;
            ro=b;
        }
        double operator()(double a, double b)
        {
            
            if (abs(a-b)<0.00001)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
};

double find_ef_pala(double ro, double B)
{
    double wyn=1;
    while (abs(roo(wyn,B)-ro)>0.001)
    {
        wyn+=0.00001;
    }
    return wyn;
}
double find_ef_prim(double ro, double B)
{
    funcion test(ro, B);
    //Scout<<B<<endl;
    war temp(B,ro);
    boost::uintmax_t k=100000000000;
    double a=1;
    double b=2;
    double fa;
    double fb;
    fa=test(a);
    fb=test(b);
    while (fa*fb>0)
    {
        a=b;
        fa=fb;
        b*=2;
        fb=(ro-roo(b,B));
    }
    pair<double,double> wyn=boost::math::tools::toms748_solve<funcion,double,war>(test,a,b,test(a),test(b),temp,k);
    return  (wyn.first+wyn.second)/2;
}
double roo(double ef)
{
    double pf=sqrt((pow(ef,2)-pow(m*v*v,2))/pow(v,2));
    double xf=pf/(m*v);
    return pow(xf,3)/(pow(M_PI,2)*pow(la,3));
}

class fun0
{
    public:
        double rr;
        fun0(double a)
        {
            rr=a;
        }
        double operator()(double ef)
        {
            return rr-roo(ef);
        }
};

double find_ef_prim(double ro)
{
    fun0 test(ro);
    war temp(0,ro);
    boost::uintmax_t k=1000000000;
    double a=m*v*v;
    double b=100*a;
    double fa;
    double fb;
    fa=test(a);
    fb=test(b);
    while (fa*fb>0)
    {
        a=b;
        fa=fb;
        b*=2;
        fb=(ro-roo(b));
    }
    pair<double,double> wyn=boost::math::tools::toms748_solve<fun0,double,war>(test,a,b,test(a),test(b),temp,k);
    return  (wyn.first+wyn.second)/2;
}

int main()
{
    int deg=3;
    int ro_n=ronum; //length of ro table
    double ro_start=0.84;    //ro starting value
    double ro_end=2;        //ro ending value
    double *ro=new double[ro_n];   
    gen_tab(ro,ro_start,ro_end,ro_n);  //fill ro table
    double B_start=0.1;
    double B_end=10;
    int B_n=bnum;
    double *BB=new double[B_n];
    gen_tab(BB,B_start,B_end,B_n);
    double gi_start=0.005;
    double gi_end=0.05;
    int gi_n=gnum;
    double *gi=new double[gi_n];
    gen_tab(gi,gi_start,gi_end,gi_n);
    double Jmaxp[ronum][gnum][bnum];
    double Efp[ronum][gnum][bnum];
    for (int b=0;b<B_n;b++)
    {
        double B=BB[b]*pow(10,17);
        for (int i=0;i<ro_n;i++)
        {
            for (int G=0;G<gi_n;G++)
            {   
                double rrp=(ro[i]/2)*(1+gi[G]);
                double ef=find_ef_prim(rrp,B);
                Efp[i][G][b]=ef;
                Jmaxp[i][G][b]=floor(jmax(ef,B));
            }
        }
    }
    double ekp[ronum][gnum][bnum];

    for (int b=0;b<B_n;b++)
    {
        double B=BB[b]*pow(10,17);
        for (int i=0;i<ro_n;i++)
        {
            for (int G=0;G<gi_n;G++)
            {
                double eki=0;
                double pq=0;
                int gg;
                double jj=Jmaxp[i][G][b];
                double efi=Efp[i][G][b];
                double xf,fu,w;
                for (int a=0;a<=jj;a++)
                {
                    if (a==0)
                    {
                        gg=1;
                    }
                    else
                    {
                        gg=3;
                    }
                    xf=sqrt(efi*efi-1-2*B/BC*a)/sqrt(1+2*a*B/BC);
                    fu=(xf*sqrt(1+xf*xf)+log(xf+sqrt(1+xf*xf)))/2;
                    w=((2*(B/BC)*m*v*v)/(4*M_PI*M_PI*la*la*la))*gg*(1+2*a*B/BC)*fu;
                    if (isnan(w))
                    {
                        eki+=0;
                    }
                    else
                    {
                        eki+=w;
                    }
                }
            ekp[i][G][b]=eki;

            }
        }
    }
    double Jmaxn[ronum][gnum][bnum];
    double Efn[ronum][gnum][bnum];
    for (int b=0;b<B_n;b++)
    {
        double B=BB[b]*pow(10,17);
        for (int i=0;i<ro_n;i++)
        {
            for (int G=0;G<gi_n;G++)
            {
                double rrn=(ro[i]/2)*(1-gi[G]);
                double ef=find_ef_prim(rrn,B);
                Efn[i][G][b]=ef;
                Jmaxn[i][G][b]=floor(jmax(ef,B));
            }
        }
    }
    double ekn[ronum][gnum][bnum];

    for (int b=0;b<B_n;b++)
    {
        double B=BB[b]*pow(10,17);
        for (int i=0;i<ro_n;i++)
        {
            for (int G=0;G<gi_n;G++)
            {
                double eki=0;
                double pq=0;
                int gg;
                double jj=Jmaxn[i][G][b];
                double efi=Efn[i][G][b];
                double xf,fu,w;
                for (int a=1;a<=jj;a++)
                {
                    if (a==0)
                    {
                        gg=1;
                    }
                    else
                    {
                        gg=3;
                    }
                    xf=sqrt(efi*efi-1-2*B/BC*a)/sqrt(1+2*a*B/BC);
                    fu=(xf*sqrt(1+xf*xf)+log(xf+sqrt(1+xf*xf)))/2;
                    w=((2*(B/BC)*m*v*v)/(4*M_PI*M_PI*la*la*la))*gg*(1+2*a*B/BC)*fu;
                    if (isnan(w))
                    {
                        eki+=0;
                    }
                    else
                    {
                        eki+=w;
                    }
                }
            ekn[i][G][b]=eki;

            }
        }
    }

    double EM[ronum][gnum][bnum];
    double M[ronum][gnum][bnum];
    for (int b=0;b<B_n;b++)
    {
        double B=BB[b]*pow(10,17);
        for (int i=0;i<ro_n;i++)
        {
            for (int G=0;G<gi_n;G++)
            {
                EM[i][G][b]=(0.299*gi[G]*ro[i]*B*pow(10,-4)*5.05*pow(10,-27))/(1.6*pow(10,-13));
                M[i][G][b]=EM[i][G][b]/B;
            }
        }
    }
    double Bint=8.99;
    double B0=400;
    double r0=0.17;
    double y=0.17;
    double Bag[ronum][gnum][bnum];

    for (int b=0;b<B_n;b++)
    {
        double B=BB[b]*pow(10,17);
        for (int i=0;i<ro_n;i++)
        {
            for (int G=0;G<gi_n;G++)
            {
                Bag[i][G][b]=Bint+(B0-Bint)*exp(-y*pow(ro[i]/r0,2));
            }
        }
    }
    double ekk[ronum][gnum][bnum];
    for (int b=0;b<B_n;b++)
    {
        for (int i=0;i<ro_n;i++)
        {
            for (int G=0;G<gi_n;G++)
            {
                ekk[i][G][b]=ekn[i][G][b]+ekp[i][G][b]+Bag[i][G][b];
            }
        }
    } 
    double Ek[ronum][bnum];
    double Gm[ronum][bnum];
    for (int b=0;b<B_n;b++)
    {
        for (int i=0;i<ro_n;i++)
        {
            double min=pow(10,175);
            double index=0;
            for (int G=0;G<gi_n;G++)
            {
                if (ekk[i][G][b]+EM[i][G][b]<=min)
                {
                    index=G;
                    min=ekk[i][G][b]+EM[i][G][b];
                }
            }
            Ek[i][b]=min;
            Gm[i][b]=index;

        }
    }
    double pp[ronum][bnum];
    double mu[ronum][bnum];
    for (int b=0;b<B_n;b++)
    {
        for (int i=0;i<ro_n-1;i++)
        {
            if (i<2 || i>ro_n-4)
            {
                pp[i][b]=(ro[i]*(Ek[i+1][b]-Ek[i][b])/(ro[i+1]-ro[i]))-Ek[i][b];
            }
            else
            {
                pp[i][b]=ro[i]*(8*Ek[i+1][b]-8*Ek[i-1][b]-Ek[i+2][b]+Ek[i-2][b])/(12*(ro[i+1]-ro[i]))-Ek[i][b];
            }
            mu[i][b]=(pp[i][b]+Ek[i][b])/ro[i];
        }
    }
    //**************************************************************
    //calculating for zero magnetic field
    double Ef[ronum];
    double ek0[ronum];
    double ekk0[ronum];
    double pp0[ronum];
    double mu0[ronum];
    double Bag0[ronum];

    for (int i=0;i<ronum;i++)
    {
        double rr=ro[i];
        Ef[i]=find_ef_prim(rr);
        //cout<<Ef[i]/(m*v*v)<<endl;
    }
    for (int i=0;i<ronum;i++)
    {
        double rr=ro[i];
        Bag0[i]=Bint+(B0-Bint)*exp(-y*pow((rr/r0),2));
        //Bag0[i]=B0;
    }

    for (int i=0;i<ronum;i++)
    {
        double rr=ro[i];
        double ef=Ef[i];
        double pf=sqrt((pow(ef,2)-pow(m*v*v,2)))/v;
        double xf=pf/(m*v);
        double X=(1/(8*M_PI*M_PI))*(xf*sqrt(1+pow(xf,2))*(1+2*pow(xf,2))-log(xf+sqrt(1+pow(xf,2))));
        ek0[i]=(3*(m*v*v)/(la*la*la))*X+Bag0[i];
    }

    for (int i=0;i<ronum;i++)
    {
        pp0[i]=((ro[i]*(ek0[i+1]-ek0[i])/(ro[i+1]-ro[i]))-ek0[i]);
        mu0[i]=(ek0[i]+pp0[i])/ro[i];
    }
    
    


    //******************************************************************
    ofstream f;
    f.precision(11);
    f.open("tab.txt");
    for (int i=0;i<5;i++)
    {f<<"#"<<endl;}
    f<<ro_n-1<<" "<<B_n+1<<endl;
    for (int i=0;i<2;i++)
    {f<<"#"<<endl;}
    f<<endl;
    for (int j=0;j<ronum-1;j++)
    {
        f<<j+1<<"\t"<<ro[j]<<"\t"<<ek0[j]*1.78*pow(10,12)<<"\t"<<pp0[j]*pow(10,33)*1.6<<"\t"<<mu0[j]<<"\t"<<0<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"<<endl;
    }
    for (int i=0;i<B_n;i++)
    {
        for (int j=0;j<ro_n-1;j++)
        {
            f<<j+1<<"\t"<<ro[j]<<"\t"<<Ek[j][i]*1.78*pow(10,12)<<"\t"<<pp[j][i]*pow(10,33)*1.6<<"\t"<<mu[j][i]<<"\t"<<BB[i]*pow(10,2)<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"<<endl;
        }
    }
    f.close();
    
    /*
    for (int i=0;i<5;i++)
    {cout<<"#"<<endl;}
    cout<<ro_n-1<<" "<<B_n<<endl;
    for (int i=0;i<2;i++)
    {cout<<"#"<<endl;}
    cout<<endl;
    for (int i=0;i<B_n;i++)
    {
        for (int j=0;j<ro_n-1;j++)
        {
            cout<<j+1<<"\t"<<ro[j]<<"\t"<<Ek[j][i]*1.78*pow(10,12)<<"\t"<<pp[j][i]*pow(10,33)*1.6<<"\t"<<mu[j][i]<<"\t"<<BB[i]*pow(10,2)<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"<<endl;
        }
    }
    */
    /*
    delete[] BB;
    delete[] gi;
    delete[] ro;
    dealloc_3d(&ekk,ro_n,B_n,gi_n);
    dealloc_3d(&Jmaxn,ro_n,B_n,gi_n);
    dealloc_3d(&Jmaxp,ro_n,B_n,gi_n);
    dealloc_3d(&ekn,ro_n,B_n,gi_n);
    */
    return 0;
}