#include <stdio.h>
#include <math.h>

#define TEST 0

#define RECURS_LEVEL_MAX  10
#define N_MAX             10

double max(double a, double b) {
    return (a>b)? a : b;
}
double min(double a, double b) {
    return (a<b)? a : b;
}
double dabs(double a) {
    return (a<0)? -1*a : a;
}

/////////////////////////////////////////////////////////////////////
// 7.1 Calki jednokrotne

typedef double (*Func1vFp)(double); 

double f_poly(double x) {  //    polynomial  a[0]+a[1]x+ ... + a[n]x^n
    double a[6] = {-6.25, 1.35, 3.5, 0, -4, 2};
    double sum = 0, curr = 1;
    for (int i = 0; i < 6; ++i) {
        sum += curr*a[i];
        curr *= x;
    }
    return sum;
}

double f_rat(double x) {  // f(x) = 1/((x-0.5)*(x-0.5)+0.01)

    double val = x-0.5;
    val *= val;
    val += 0.01;
    val = 1/val;

    return val;
}

double f_exp(double x) {  
    double val = 2*x*exp(-1.5*x)-1;

    return val;
}

double f_trig(double x) {  
    double val = x*tan(x);

    return val;
}

double quad_rect_left(Func1vFp f1, double a, double b, int n) {
    const double interval = (b-a)/n;
    double result = 0, c = a;

    for (int i = 0; i < n; ++i, c += interval) {
        result += f1(c);
    }
    result *= interval;

    return result;
}

double quad_rect_right(Func1vFp f1, double a, double b, int n) {
    const double interval = (b-a)/n;
    double result = 0, c = interval+a;

    for (int i = 0; i < n; ++i, c += interval) {
        result += f1(c);
    }
    result *= interval;

    return result;
}

double quad_rect_mid(Func1vFp f1, double a, double b, int n) {
    const double interval = (b-a)/n;
    const double step = interval/2;
    double result = 0, c = a+step;

    for (int i = 0; i < n; ++i, c += interval) {
        result += f1(c);
    }
    result *= interval;

    return result;
}

double quad_trap(Func1vFp func, double a, double b, int n) {
    const double interval = (b-a)/n;
    const double step = interval/2;
    double result = 0, c1 = a, c2 = c1+interval;

    double func_val[n+1];
    func_val[0] = func(c1);

    for (int i = 0; i < n; ++i, c1 = c2, c2 += interval) {
        func_val[i+1] = func(c2);
        result += func_val[i]+func_val[i+1];
    }
    result *= step;

    return result;
}

double quad_simpson(Func1vFp f, double a, double b, int n) {
    const double interval = (b-a)/n;
    const double step = interval/2;
    const double h = interval/6;
    double result = 0, c1 = a, c2 = c1+interval, c_mid = c1+step;

    double func_val[n+1];
    func_val[0] = f(c1);

    for (int i = 0; i < n; ++i, c1 = c2, c2 += interval, c_mid += interval) {
        func_val[i+1] = f(c2);
        result += func_val[i]+4*f(c_mid)+func_val[i+1];
    }
    result *= h;

    return result;
}


typedef double (*QuadratureFp)(Func1vFp, double, double, int);  


Func1vFp func_tab[4]={f_poly, f_rat, f_trig, f_exp};

QuadratureFp quad_tab[5]={quad_rect_left, quad_rect_right, quad_rect_mid, quad_trap, quad_simpson};


double quad_select(int fun_no, int quad_no, double a, double b, int n) {
    QuadratureFp quadf = quad_tab[quad_no];
    Func1vFp fun = func_tab[fun_no];

    return quadf(fun, a, b, n);
}


double recurs(Func1vFp f, double a, double b, double S, double delta, QuadratureFp quad, int level) {
    double S1 = quad(f, a, (a+b)/2, 1), S2 = quad(f, (a+b)/2, b, 1);

    if (dabs(S-S1-S2) <= delta) return S1+S2;
    if (level == RECURS_LEVEL_MAX) return NAN;

    double rek1 = recurs(f, a, (a+b)/2, S1, delta/2, quad, level+1);
    double rek2 = recurs(f, (a+b)/2, b, S2, delta/2, quad, level+1);

    if (rek1 == NAN || rek2 == NAN) return NAN;
    return rek1+rek2;
}


double init_recurs(Func1vFp f, double a, double b, double delta, QuadratureFp quad) {
    double S = quad(f, a, b, 1);
    return recurs(f, a, b, S, delta, quad, 0);
}

///////////////////////////////////////////////////////////////
// 7.2 Calki dwukrotne:

typedef double (*Func2vFp)(double, double);


double func2v_1(double x, double y) {
    return sin(x/y);
}
double func2v_2(double x, double y) {
    return 2-x*x-y*y*y;
}


double lower_bound1(double x) {
    return sin(x);
}
double upper_bound1(double x) {
    return x*x+1;
}
double lower_bound2(double x) {
    return 0.7*exp(-2*x*x);
}
double upper_bound2(double x) {
    return sin(10*x);
}


double dbl_integr(Func2vFp f, double x1, double x2, int nx, double y1, double y2, int ny)  {
    const double intervaly = (y2-y1)/ny, intervalx = (x2-x1)/nx;
    double result = 0;

    for (int i = 0; i < ny; ++i) {
        double y_const = intervaly*i+y1;

        for (int j = 0; j < nx; ++j) {
            double c = intervalx*j+x1;
            result += f(c, y_const);
        }

    }

    result *= intervaly;
    result *= intervalx;
    return result;
}


double dbl_integr_normal_1(Func2vFp f, double x1, double x2, int nx, double hy, Func1vFp fg, Func1vFp fh)  {
    const double intervalx = (x2-x1)/nx;
    const double stepx = intervalx/2;
    double result = 0;

    for (int i = 0; i < nx; ++i) {
        double a_x1 = intervalx*i+x1;
        a_x1 += stepx;
        double y1 = fg(a_x1), y2 = fh(a_x1);
        int ny = ceil((y2-y1)/hy);
        double intervaly = (y2-y1)/ny;
        double stepy = intervaly/2;
        double resulty = 0;
        for (int j = 0; j < ny; ++j) {
            double c = intervaly*j+y1+stepy;
            resulty += f(a_x1, c);
        }
        resulty *= intervaly;
        result += resulty;
    }

    result *= intervalx;

    return result;
}


double dbl_integr_normal_n(Func2vFp f, double x1, double x2, int nx, double y1, double y2, int ny, Func1vFp fg, Func1vFp fh)  {
    const double dx = (x2-x1)/nx;

    double result = 0;
    double hy = (y2-y1)/ny;

    double a_x1 = x1;

    for (int i = 0; i < nx; ++i, a_x1 += dx) {
        double y_lower = max(y1, fg(a_x1)), y_upper = min(y2, fh(a_x1));
        if (y_lower > y_upper) continue;
        int dny = ceil((y_upper-y_lower)/hy);
        double dy = (y_upper-y_lower)/dny;

        double c = y_lower;
        double resulty = 0;
        for (int j = 0; j < dny; ++j, c += dy) {
            resulty += f(a_x1, c);
        }
        resulty *= dy;
        result += resulty;
    }

    result *= dx;

    return result;
}

///////////////////////////////////////////////////////////////
// 7.3 Calki wielokrotne:

typedef double (*FuncNvFp)(const double*, int);
typedef int (*BoundNvFp)(const double*, int);


double func3v(const double v[], int n) {
    return v[0]-v[1]+2*v[2];
}


int bound3v(const double v[], int n) {
    if(v[0]>0 && v[0]<0.5 && v[1]*v[1]+(v[2]-1)*(v[2]-1)<1) return 1; 
    return 0;
}


double funcNv(const double v[], int n) {
    double fv=1.;
    for(int i=1; i<n; ++i) fv += sin(i*v[i]);   
    return fv;
}

int boundNv(const double v[], int n) {
    double r=0.0;
    for(int i=0; i<n; ++i) r += (v[i]-1)*(v[i]-1); 
    if(r > 1.) return 0;
    return 1;
}

// Obliczanie calek wielokrotnych

double trpl_quad_rect(FuncNvFp f, const double variable_lim[][2], const int tn[], BoundNvFp boundary) {
    int n = 3;
    double dvar[n];
    double coord[n];

    for (int i = 0; i < n; ++i) {
        dvar[i] = (variable_lim[i][1]-variable_lim[i][0])/tn[i];
    }

    double result = 0;

    coord[0] = variable_lim[0][0]+dvar[0];
    for (int i = 0; i < tn[0]; ++i, coord[0] += dvar[0]) {
        coord[1] = variable_lim[1][0]+dvar[1];
        for (int j = 0; j < tn[1]; ++j, coord[1] += dvar[1]) {
            coord[2] = variable_lim[2][0]+dvar[2];
            for (int k = 0; k < tn[2]; ++k, coord[2] += dvar[2]) {
                if (boundary == NULL || boundary(coord, n)) {
                    result += f(coord, n);
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) result *= dvar[i];

    return result;
}

void recur_quad_rect_mid(double *psum, FuncNvFp f, int variable_no, double tvariable[], const double variable_lim[][2], const int tn[], int level, BoundNvFp boundary) {
    double dvar = (variable_lim[level][1]-variable_lim[level][0])/tn[level];

    double dvar_product = dvar;
    for (int i = 0; i < level; ++i) dvar_product *= (variable_lim[i][1]-variable_lim[i][0])/tn[i];

    tvariable[level] = variable_lim[level][0]+dvar/2;
    for (int i = 0; i < tn[level]; ++i, tvariable[level] += dvar) {
        if (level == variable_no-1) {
            if (boundary == NULL || boundary(tvariable, variable_no)) *psum += dvar_product*f(tvariable, variable_no);
        } else {
            recur_quad_rect_mid(psum, f, variable_no, tvariable, variable_lim, tn, level+1, boundary);
        }
    }
}

int main(void)
{
    int to_do, n, nx, ny, integrand_fun_no, method_no, flag;
    double a,b,x1,x2,y1,y2,hy,sum,delta;
    double tvariable[N_MAX], variable_lim[N_MAX][2];
    int tn[N_MAX];

    if (TEST) printf("Wpisz numer testu [1, 7]: ");
    scanf("%d", &to_do);
    switch (to_do) {
    case 1:
        if (TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow: ");
        scanf("%lf %lf %d", &a, &b, &n);
        for(int q=0; q<5; ++q) {
            for(int f=0; f<4; ++f) {
                printf("%.5f ",quad_select(f, q, a, b, n));
            }
            printf("\n");
        }
        break;
    case 2:
        if(TEST) printf("Nr funkcji (0-3) i metody (0-4): ");
        scanf("%d %d",&integrand_fun_no,&method_no);
        if (TEST) printf("Wpisz przedzial calkowania i tolerancje bledu: ");
        scanf("%lf %lf %lf", &a, &b, &delta);
        printf("%.5f\n",init_recurs(func_tab[integrand_fun_no],a,b,delta,quad_tab[method_no]));
        break;
    case 3:
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow wzdluz x: ");
        scanf("%lf %lf %d",&x1,&x2,&nx);
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow wzdluz y: ");
        scanf("%lf %lf %d",&y1,&y2,&ny);
        printf("%.5f\n",dbl_integr(func2v_2, x1, x2, nx, y1, y2, ny));
        break;
    case 4:
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow zmiennej x: ");
        scanf("%lf %lf %d",&x1,&x2,&nx);
        if(TEST) printf("Wpisz dlugosc podprzedzialu wzdluz y: ");
        scanf("%lf",&hy);
        printf("%.5f\n",dbl_integr_normal_1(func2v_2, x1, x2, nx, hy, lower_bound2, upper_bound2));
        break;
    case 5:
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow wzdluz x: ");
        scanf("%lf %lf %d",&x1,&x2,&nx);
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow wzdluz y: ");
        scanf("%lf %lf %d",&y1,&y2,&ny);
        printf("%.5f\n",dbl_integr_normal_n(func2v_2, x1, x2, nx, y1, y2, ny, lower_bound2, upper_bound2));
        break;
    case 6:
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow 1. zmiennej: ");
        scanf("%lf %lf %d",&variable_lim[0][0],&variable_lim[0][1],tn);
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow 2. zmiennej: ");
        scanf("%lf %lf %d",&variable_lim[1][0],&variable_lim[1][1],tn+1);
        if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow 3. zmiennej: ");
        scanf("%lf %lf %d",&variable_lim[2][0],&variable_lim[2][1],tn+2);
        if(TEST) printf("Wpisz 1 gdy ograniczenie obszaru całkowania ma byc aktywne; 0 - w przeciwnym przypadku: ");
        scanf("%d",&flag);
        printf("%.5f\n",trpl_quad_rect(func3v, variable_lim, tn, flag?bound3v:NULL));
        break;
    case 7:
        if(TEST) printf("Wpisz liczbe zmiennych <= %d: ", N_MAX);
        scanf("%d",&n);
        if(n>N_MAX) break;
        for(int i=0; i<n; ++i) {
            if(TEST) printf("Wpisz przedzial calkowania i liczbe podprzedzialow %d. zmiennej: ",i+1);
            scanf("%lf %lf %d",&variable_lim[i][0],&variable_lim[i][1],tn+i);
        }
        if(TEST) printf("Wpisz 1 gdy ograniczenie obszaru całkowania ma byc aktywne; 0 - w przeciwnym przypadku: ");
        scanf("%d",&flag);
        sum=0.;
        recur_quad_rect_mid(&sum, funcNv, n, tvariable, variable_lim, tn, 0, flag?boundNv:NULL);
        printf("%.5f", sum);
        break;
    default:
        printf("Numer testu spoza zakresu [1, 7] %d", to_do);
    }
    return 0;
}

