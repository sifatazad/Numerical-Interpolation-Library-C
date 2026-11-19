#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#define Abs(a) ((a)>=0?(a):-(a))

struct data {
    double x;
    double y;
};

static int compare(const void *a, const void *b) {
    const struct data *p= (const struct data*) a;
    const struct data *q= (const struct data*) b;

    if (p->x > q->x)
        return 1;
    if (p->x < q->x)
        return -1;
    return 0;
}

static int total_point(FILE *fp) {
    int n=0;
    char line[500];
    while (fgets(line,500,fp)!=NULL) {
        n++;
    }
    rewind(fp);
    return n;
}

static double sum_x_pow(const double *x,int n,int p) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=pow(x[i],p*1.0);
    }
    return s;
}

static double sum_y_x_pow(const double *x,const double *y,const int n,const int p) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=pow(x[i],p*1.0)*y[i];
    }
    return s;
}

static void lu_decomposition_without_file(double **matrix, double *soln, const int n) {
    for (int step=0; step<n-1; step++) {
        double col_max=Abs(matrix[step][step]);

        for (int i=step+1; i<n; i++) {
            if (Abs(matrix[i][step])>Abs(col_max)) {
                col_max=Abs(matrix[i][step]);
            }
        }

        if (col_max<1e-10) {
            printf("Singular Matrix. No unique Solution!\n");
            return;
        }

        if (Abs(matrix[step][step])!=col_max) {
            for (int i=step+1; i<n; i++) {
                if (col_max==Abs(matrix[i][step])) {
                    double *temp= malloc(sizeof(double)*(n+1));
                    if (!temp) {
                        return;
                    }

                    for (int it=0; it<n+1; it++) {
                        temp[it]=matrix[i][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[i][it]=matrix[step][it];
                    }
                    for (int it=0; it<n+1; it++) {
                        matrix[step][it]=temp[it];
                    }
                    // pivoting done. now gaussian kora need
                    free(temp);
                    break;
                }
            }
        }

        // gaussian elimination here
        for (int i=step+1; i<n; i++) {
            double factor= matrix[i][step]/matrix[step][step];
            for(int j=step; j<n; j++) {
                if (j==step) {
                    matrix[i][j]= factor;
                }
                else {
                    matrix[i][j]+= -factor*matrix[step][j];
                }
            }
        }
    }

    if (Abs(matrix[n-1][n-1])<1e-10) {
        printf("Singular Matrix. No unique Solution!\n");
        return;
    }

    double *y= malloc(sizeof(double)*n);
    if (!y) {
        return;
    }
    for (int i=0; i<n; i++) {
        y[i]=0;
    }

    for (int i=0; i<n; i++) {
        double s=0;
        for (int j=0; j<i; j++) {
            s+=matrix[i][j]*y[j];
        }
        y[i]=matrix[i][n]-s;
    }

    for (int i=0; i<n; i++) {
        soln[i]=0;
    }

    for (int i=n-1; i>=0; i--) {
        double s=0;
        for (int j=n-1; j>=i+1; j--) {
            s+=matrix[i][j]*soln[j];
        }
        soln[i]=(y[i]-s)/matrix[i][i];
    }

    free(y);
}

static void polyfit(FILE *fp,double *soln, const int fit_rank) {
    const int n= total_point(fp);
    double **matrix= malloc(sizeof(double*)*(fit_rank+1));
    if (!matrix) {
        return;
    }
    for (int i=0; i<=fit_rank; i++) {
        matrix[i]=malloc(sizeof(double)*(fit_rank+2));
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    double *x,*y;
    x=malloc(sizeof(double)*n);
    if (!x) {
        for (int i=0; i<=fit_rank; i++) {
            free(matrix[i]);
        }
        free(matrix);
        return;
    }
    y=malloc(sizeof(double)*n);
    if (!y) {
        for (int i=0; i<=fit_rank; i++) {
            free(matrix[i]);
        }
        free(matrix);
        free(x);
        return;
    }

    for (int i=0; i<n; i++) {
        fscanf(fp,"%lf\t%lf",&x[i],&y[i]);
    }
    rewind(fp);

    for (int i=0; i<=fit_rank; i++) {
        for (int j=0; j<=fit_rank; j++) {
            matrix[i][j]=sum_x_pow(x,n,i+j);
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        matrix[i][fit_rank+1]=sum_y_x_pow(x,y,n,i);
    }

    lu_decomposition_without_file(matrix,soln,fit_rank+1);

    // freeing here
    for (int i=0; i<=fit_rank; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(x);
    free(y);
}

static void polyfit_without_file(const double *x, const double *y, double *soln,
    int n, int fit_rank) {
    double **matrix= malloc(sizeof(double*)*(fit_rank+1));
    if (!matrix) {
        return;
    }
    for (int i=0; i<=fit_rank; i++) {
        matrix[i]=malloc(sizeof(double)*(fit_rank+2));
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        for (int j=0; j<=fit_rank; j++) {
            matrix[i][j]=sum_x_pow(x,n,i+j);
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        matrix[i][fit_rank+1]=sum_y_x_pow(x,y,n,i);
    }

    lu_decomposition_without_file(matrix,soln,fit_rank+1);

    // freeing here
    for (int i=0; i<=fit_rank; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

static double interpolate(FILE *fp, const double x_r, char *s) {
    int n= total_point(fp);

    double *x=malloc(sizeof(double)*n);
    if (!x) {
        return NAN;
    }

    double *y=malloc(sizeof(double)*n);
    if (!y) {
        free(x);
        return NAN;
    }

    for (int i=0; i<n; i++) {
        fscanf(fp,"%lf\t%lf",&x[i],&y[i]);
    }
    rewind(fp);

    struct data *point= malloc(sizeof(struct data)*n);
    if (!point) {
        free(x);
        free(y);
        return NAN;
    }

    for (int i=0; i<n; i++) {
        (point+i)->x=x[i];
        (point+i)->y=y[i];
    }

    qsort(point,n,sizeof(struct data),compare);

    for (int i=0; i<n; i++) {
        x[i]= point[i].x;
        y[i]= point[i].y;
    }

    int len= strlen(s);
    for (int i=0; i<len; i++) {
        if (s[i]>=65 && s[i]<=90)
            s[i]+=32;
    }

    if (x_r<x[0] || x_r>x[n-1]) {
        printf("Error! Extrapolation can't be done!\n");
        free(x);
        free(y);
        free(point);
        return NAN;
    }

    if (!strcmp(s,"linear")) {
        if (n<=1) {
            printf("Linear Interpolation not possible!\n");
            printf("At least 2 points required.\n");
            free(x);
            free(y);
            free(point);
            return NAN;
        }
        else {
            for (int i=0; i<n; i++) {
                if (x_r<x[i]) {
                    double val= y[i-1]+(y[i]-y[i-1])/(x[i]-x[i-1])*(x_r-x[i-1]);
                    free(x);
                    free(y);
                    free(point);
                    return val;
                }
                else if (x_r==x[i]) {
                    double val= y[i];
                    free(x);
                    free(y);
                    free(point);
                    return val;
                }
            }
        }
    }

    else if (!strcmp(s,"spline")) {
        if (n<3) {
            printf("Spline Interpolation not possible!\n");
            printf("At least 3 points required.\n");
            free(x);
            free(y);
            free(point);
            return NAN;
        }
        else {
            double **matrix= malloc(sizeof(double*)*n);
            if (!matrix) {
                free(x);
                free(y);
                free(point);
                return NAN;
            }

            for (int i=0; i<n; i++) {
                matrix[i]=malloc(sizeof(double)*(n+1));
                if (!matrix[i]) {
                    free(x);
                    free(y);
                    free(point);
                    free(matrix);
                    return NAN;
                }
            }

            double *h= malloc(sizeof(double)*(n-1));
            if (!h) {
                free(x);
                free(y);
                free(point);
                for (int i=0; i<n; i++) {
                    free(matrix[i]);
                }
                free(matrix);
                return NAN;
            }

            for (int i=0; i<n-1; i++) {
                h[i]=x[i+1]-x[i];
            }

            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    if (i==0 || i==n-1) {
                        if (i==j) {
                            matrix[i][j]=1;
                        }
                        else {
                            matrix[i][j]=0;
                        }
                    }
                    else {
                        if (i==j) {
                            matrix[i][j]=2*(h[i]+h[i-1]);
                            matrix[i][j-1]=h[i-1];
                            matrix[i][j+1]=h[i];
                            j++;
                        }
                        else {
                            matrix[i][j]=0;
                        }
                    }
                }
            }

            // for (int i=0; i<n; i++) {
            //     for (int j=0; j<n; j++) {
            //         printf("%.3lf\t",matrix[i][j]);
            //     }
            //     printf("\n");
            // }

            for (int i=0; i<n; i++) {
                if (i==0 || i==n-1) {
                    matrix[i][n]=0;
                }
                else {
                    matrix[i][n]= (3.0/h[i])*(y[i+1]-y[i])- (3.0/h[i-1])*(y[i]-y[i-1]);
                }
            }

            double *c= malloc(sizeof(double)*n);
            if (!c) {
                free(x);
                free(y);
                free(point);
                for (int i=0; i<n; i++) {
                    free(matrix[i]);
                }
                free(matrix);
                free(h);
                return NAN;
            }

            lu_decomposition_without_file(matrix,c,n);

            // for (int i=0; i<n; i++) {
            //     printf("%.6lf\n",c[i]);
            // }

            double *b= malloc(sizeof(double)*n);
            if (!b) {
                free(x);
                free(y);
                free(point);
                for (int i=0; i<n; i++) {
                    free(matrix[i]);
                }
                free(matrix);
                free(h);
                free(c);
                return NAN;
            }

            double *d= malloc(sizeof(double)*n);
            if (!d) {
                free(x);
                free(y);
                free(point);
                for (int i=0; i<n; i++) {
                    free(matrix[i]);
                }
                free(matrix);
                free(h);
                free(b);
                free(c);
                return NAN;
            }

            for (int i=0; i<n-1; i++) {
                b[i]= (y[i+1]-y[i])/(h[i])-(h[i]/3.0)*(2*c[i]+c[i+1]);
                d[i]= (c[i+1]-c[i])/(3*h[i]);
            }
            b[n-1]=0;
            d[n-1]=0;

            for (int i=0; i<n; i++) {
                if (x[i]>x_r) {
                    i--;
                    double val= y[i]+b[i]*(x_r-x[i])+c[i]*(x_r-x[i])*(x_r-x[i])
                    +d[i]*(x_r-x[i])*(x_r-x[i])*(x_r-x[i]);
                    free(x);
                    free(y);
                    free(point);
                    for (int ii=0; ii<n; ii++) {
                        free(matrix[ii]);
                    }
                    free(matrix);
                    free(h);
                    free(b);
                    free(c);
                    free(d);
                    return val;
                }
                else if (x[i]==x_r) {
                    double val= y[i];
                    free(x);
                    free(y);
                    free(point);
                    for (int ii=0; ii<n; ii++) {
                        free(matrix[ii]);
                    }
                    free(matrix);
                    free(h);
                    free(b);
                    free(c);
                    free(d);
                    return val;
                }
            }
        }
    }
    else if (!strcmp(s,"pchip")) {
        if (n<3) {
            printf("Pchip Interpolation not possible!\n");
            printf("At least 3 points required.\n");
            free(x);
            free(y);
            free(point);
            return NAN;
        }

        double *h= malloc(sizeof(double)*(n-1));
        if (!h) {
            free(x);
            free(y);
            free(point);
            return NAN;
        }

        for (int i=0; i<n-1; i++) {
            h[i]=x[i+1]-x[i];
        }

        double *delta= malloc(sizeof(double)*(n-1));
        if (!delta) {
            free(x);
            free(y);
            free(point);
            free(h);
            return NAN;
        }

        double *m= malloc(sizeof(double)*n);
        if (!m) {
            free(x);
            free(y);
            free(point);
            free(h);
            free(delta);
            return NAN;
        }

        int interval=0;
        for (int i=0; i<n; i++) {
            if (x[i]>x_r) {
                interval=i-1;
                break;
            }
            else if (x[i]==x_r) {
                double val=y[i];
                free(x);
                free(y);
                free(point);
                free(h);
                free(delta);
                free(m);
                return val;
            }
        }

        for (int i=0; i<n-1; i++) {
            delta[i]= (y[i+1]-y[i])/(x[i+1]-x[i]);
        }

        double c,d;
        for (int i=0; i<n; i++) {
            if (i==0) {
                m[0]=((2*h[0]+h[1])*delta[0]-h[0]*delta[1])/(h[0]+h[1]);
                if (m[0]*delta[0]<=0) {
                    m[0]=0;
                }
                else {
                    if (Abs(m[0])>3.0*Abs(delta[0])) {
                        m[0]=3.0*delta[0];
                    }
                }
            }
            else if (i==n-1) {
                m[n-1]=((2*h[n-2]+h[n-3])*delta[n-2]-h[n-2]*delta[n-3])/(h[n-2]+h[n-3]);
                if (m[n-1]*delta[n-2]<=0) {
                    m[n-1]=0;
                }
                else {
                    if (Abs(m[n-1])>3.0*Abs(delta[n-2])) {
                        m[n-1]=3.0*delta[n-2];
                    }
                }
            }
            else {
                if (delta[i-1]*delta[i]<=0) {
                    m[i]=0;
                }
                else {
                    m[i]=3*(h[i]+h[i-1])/((2*h[i]+h[i-1])/
                        delta[i-1]+(h[i]+2*h[i-1])/delta[i]);
                }
            }
        }

        c= (3*delta[interval]-2*m[interval]-m[interval+1])/h[interval];
        d= (m[interval]+m[interval+1]-2*delta[interval])/(h[interval]*h[interval]);

        double val= y[interval]+ m[interval]*(x_r-x[interval])+
            c*(x_r-x[interval])*(x_r-x[interval])+
                d*(x_r-x[interval])*(x_r-x[interval])*(x_r-x[interval]);

        //
        free(x);
        free(y);
        free(point);
        free(delta);
        free(h);
        free(m);
        return val;
    }
    else if (!strcmp(s,"nearest")) {
        if (n<2) {
            printf("Error! At least two points required!\n");
            free(x);
            free(y);
            free(point);
            return NAN;
        }

        for(int i=0; i<n; i++) {
            if (x_r<x[i]) {
                i--;
                if (Abs(x_r-x[i])<=Abs(x[i+1]-x_r)) {
                    double val= y[i];
                    free(x);
                    free(y);
                    free(point);
                    return val;
                }
                else {
                    double val= y[i+1];
                    free(x);
                    free(y);
                    free(point);
                    return val;
                }
            }
            else if (x_r==x[i]) {
                double val= y[i];
                free(x);
                free(y);
                free(point);
                return val;
            }
        }
    }
    else if (!strcmp(s,"makima")) {
        if (n<3) {
            printf("Makima Interpolation not possible!\n");
            printf("At least 3 points required.\n");
            free(x);
            free(y);
            free(point);
            return NAN;
        }

        double *h= malloc(sizeof(double)*(n-1));
        if (!h) {
            free(x);
            free(y);
            free(point);
            return NAN;
        }

        for (int i=0; i<n-1; i++) {
            h[i]=x[i+1]-x[i];
        }

        double *delta= malloc(sizeof(double)*(n+3));
        if (!delta) {
            free(x);
            free(y);
            free(point);
            free(h);
            return NAN;
        }

        double *m= malloc(sizeof(double)*n);
        if (!m) {
            free(x);
            free(y);
            free(point);
            free(h);
            free(delta);
            return NAN;
        }

        int interval=0;
        for (int i=0; i<n; i++) {
            if (x[i]>x_r) {
                interval=i-1;
                break;
            }
            else if (x[i]==x_r) {
                double val=y[i];
                free(x);
                free(y);
                free(point);
                free(h);
                free(delta);
                free(m);
                return val;
            }
        }

        for (int i=2; i<n+1; i++) {
            delta[i]= (y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
        }

        delta[0]= 3*delta[2]-2*delta[3];
        delta[1]= 2*delta[2]-delta[3];
        delta[n+1]= 2*delta[n]-delta[n-1];
        delta[n+2]= 3*delta[n]-2*delta[n-1];

        double c,d;
        for (int i=2; i<n+2; i++) {
            double w1= Abs(delta[i+1]-delta[i])+ Abs(delta[i+1]+delta[i])/2;
            double w2= Abs(delta[i-1]-delta[i-2])+ Abs(delta[i-1]+delta[i-2])/2;

            if (w1+w2!=0) {
                m[i-2]= (w1*delta[i-1]+w2*delta[i])/(w1+w2);
            }
            else {
                m[i-2]= (delta[i-1]+delta[i])/2;
            }
        }

        c= (3*delta[interval+2]-2*m[interval]-m[interval+1])/h[interval];
        d= (m[interval]+m[interval+1]-2*delta[interval+2])/(h[interval]*h[interval]);

        double val= y[interval]+ m[interval]*(x_r-x[interval])+
            c*(x_r-x[interval])*(x_r-x[interval])+
                d*(x_r-x[interval])*(x_r-x[interval])*(x_r-x[interval]);

        //
        free(x);
        free(y);
        free(point);
        free(delta);
        free(h);
        free(m);
        return val;
    }
    else if (!strcmp(s,"polyfit")) {
        double *co_eff= malloc(sizeof(double)*4);
        if (!co_eff) {
            free(x);
            free(y);
            free(point);
            return NAN;
        }

        polyfit_without_file(x,y,co_eff,n,3);

        double val= 0;
        for (int i=0; i<4; i++) {
            val+= co_eff[i]*pow(x_r,i);
        }

        free(x);
        free(y);
        free(point);
        free(co_eff);
        return val;
    }
    else {
        printf("Error! Type a valid method.\n");
        // free
        free(x);
        free(y);
        free(point);
        return NAN;
    }
    free(x);
    free(y);
    free(point);
    return NAN;
}

static double interpolate_without_file(double *x, double *y,
    const int n, const double x_r, char *s) {
    struct data *point= malloc(sizeof(struct data)*n);
    if (!point) {
        return NAN;
    }

    for (int i=0; i<n; i++) {
        (point+i)->x=x[i];
        (point+i)->y=y[i];
    }

    qsort(point,n,sizeof(struct data),compare);

    for (int i=0; i<n; i++) {
        x[i]= point[i].x;
        y[i]= point[i].y;
    }

    int len= strlen(s);
    for (int i=0; i<len; i++) {
        if (s[i]>=65 && s[i]<=90)
            s[i]+=32;
    }

    if (x_r<x[0] || x_r>x[n-1]) {
        printf("Error! Extrapolation can't be done!\n");
        free(point);
        return NAN;
    }

    if (!strcmp(s,"linear")) {
        if (n<=1) {
            printf("Linear Interpolation not possible!\n");
            printf("At least 2 points required.\n");
            free(point);
            return NAN;
        }
        else {
            for (int i=0; i<n; i++) {
                if (x_r<x[i]) {
                    double val= y[i-1]+(y[i]-y[i-1])/(x[i]-x[i-1])*(x_r-x[i-1]);
                    free(point);
                    return val;
                }
                else if (x_r==x[i]) {
                    double val= y[i];
                    free(point);
                    return val;
                }
            }
        }
    }

    else if (!strcmp(s,"spline")) {
        if (n<3) {
            printf("Spline Interpolation not possible!\n");
            printf("At least 3 points required.\n");
            free(point);
            return NAN;
        }
        else {
            double **matrix= malloc(sizeof(double*)*n);
            if (!matrix) {
                free(point);
                return NAN;
            }

            for (int i=0; i<n; i++) {
                matrix[i]=malloc(sizeof(double)*(n+1));
                if (!matrix[i]) {
                    free(point);
                    free(matrix);
                    return NAN;
                }
            }

            double *h= malloc(sizeof(double)*(n-1));
            if (!h) {
                free(point);
                for (int i=0; i<n; i++) {
                    free(matrix[i]);
                }
                free(matrix);
                return NAN;
            }

            for (int i=0; i<n-1; i++) {
                h[i]=x[i+1]-x[i];
            }

            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    if (i==0 || i==n-1) {
                        if (i==j) {
                            matrix[i][j]=1;
                        }
                        else {
                            matrix[i][j]=0;
                        }
                    }
                    else {
                        if (i==j) {
                            matrix[i][j]=2*(h[i]+h[i-1]);
                            matrix[i][j-1]=h[i-1];
                            matrix[i][j+1]=h[i];
                            j++;
                        }
                        else {
                            matrix[i][j]=0;
                        }
                    }
                }
            }

            for (int i=0; i<n; i++) {
                if (i==0 || i==n-1) {
                    matrix[i][n]=0;
                }
                else {
                    matrix[i][n]= (3.0/h[i])*(y[i+1]-y[i])- (3.0/h[i-1])*(y[i]-y[i-1]);
                }
            }

            double *c= malloc(sizeof(double)*n);
            if (!c) {
                free(point);
                for (int i=0; i<n; i++) {
                    free(matrix[i]);
                }
                free(matrix);
                free(h);
                return NAN;
            }

            lu_decomposition_without_file(matrix,c,n);

            // for (int i=0; i<n; i++) {
            //     printf("%.6lf\n",c[i]);
            // }

            double *b= malloc(sizeof(double)*n);
            if (!b) {
                free(point);
                for (int i=0; i<n; i++) {
                    free(matrix[i]);
                }
                free(matrix);
                free(h);
                free(c);
                return NAN;
            }

            double *d= malloc(sizeof(double)*n);
            if (!d) {
                free(point);
                for (int i=0; i<n; i++) {
                    free(matrix[i]);
                }
                free(matrix);
                free(h);
                free(b);
                free(c);
                return NAN;
            }

            for (int i=0; i<n-1; i++) {
                b[i]= (y[i+1]-y[i])/(h[i])-(h[i]/3.0)*(2*c[i]+c[i+1]);
                d[i]= (c[i+1]-c[i])/(3*h[i]);
            }
            b[n-1]=0;
            d[n-1]=0;

            for (int i=0; i<n; i++) {
                if (x[i]>x_r) {
                    i--;
                    double val= y[i]+b[i]*(x_r-x[i])+c[i]*(x_r-x[i])*(x_r-x[i])
                    +d[i]*(x_r-x[i])*(x_r-x[i])*(x_r-x[i]);
                    free(point);
                    for (int ii=0; ii<n; ii++) {
                        free(matrix[ii]);
                    }
                    free(matrix);
                    free(h);
                    free(b);
                    free(c);
                    free(d);
                    return val;
                }
                else if (x[i]==x_r) {
                    double val= y[i];
                    free(point);
                    for (int ii=0; ii<n; ii++) {
                        free(matrix[ii]);
                    }
                    free(matrix);
                    free(h);
                    free(b);
                    free(c);
                    free(d);
                    return val;
                }
            }
        }
    }
    else if (!strcmp(s,"pchip")) {
        if (n<3) {
            printf("Pchip Interpolation not possible!\n");
            printf("At least 3 points required.\n");
            free(point);
            return NAN;
        }

        double *h= malloc(sizeof(double)*(n-1));
        if (!h) {
            free(point);
            return NAN;
        }

        for (int i=0; i<n-1; i++) {
            h[i]=x[i+1]-x[i];
        }

        double *delta= malloc(sizeof(double)*(n-1));
        if (!delta) {
            free(point);
            free(h);
            return NAN;
        }

        double *m= malloc(sizeof(double)*n);
        if (!m) {
            free(point);
            free(h);
            free(delta);
            return NAN;
        }

        int interval=0;
        for (int i=0; i<n; i++) {
            if (x[i]>x_r) {
                interval=i-1;
                break;
            }
            else if (x[i]==x_r) {
                double val=y[i];
                free(point);
                free(h);
                free(delta);
                free(m);
                return val;
            }
        }

        for (int i=0; i<n-1; i++) {
            delta[i]= (y[i+1]-y[i])/(x[i+1]-x[i]);
        }

        double c,d;
        for (int i=0; i<n; i++) {
            if (i==0) {
                m[0]=((2*h[0]+h[1])*delta[0]-h[0]*delta[1])/(h[0]+h[1]);
                if (m[0]*delta[0]<=0) {
                    m[0]=0;
                }
                else {
                    if (Abs(m[0])>3.0*Abs(delta[0])) {
                        m[0]=3.0*delta[0];
                    }
                }
            }
            else if (i==n-1) {
                m[n-1]=((2*h[n-2]+h[n-3])*delta[n-2]-h[n-2]*delta[n-3])/(h[n-2]+h[n-3]);
                if (m[n-1]*delta[n-2]<=0) {
                    m[n-1]=0;
                }
                else {
                    if (Abs(m[n-1])>3.0*Abs(delta[n-2])) {
                        m[n-1]=3.0*delta[n-2];
                    }
                }
            }
            else {
                if (delta[i-1]*delta[i]<=0) {
                    m[i]=0;
                }
                else {
                    m[i]=3*(h[i]+h[i-1])/((2*h[i]+h[i-1])/
                        delta[i-1]+(h[i]+2*h[i-1])/delta[i]);
                }
            }
        }

        c= (3*delta[interval]-2*m[interval]-m[interval+1])/h[interval];
        d= (m[interval]+m[interval+1]-2*delta[interval])/(h[interval]*h[interval]);

        double val= y[interval]+ m[interval]*(x_r-x[interval])+
            c*(x_r-x[interval])*(x_r-x[interval])+
                d*(x_r-x[interval])*(x_r-x[interval])*(x_r-x[interval]);

        //
        free(point);
        free(delta);
        free(h);
        free(m);
        return val;
    }
    else if (!strcmp(s,"nearest")) {
        if (n<2) {
            printf("Error! At least two points required!\n");
            free(point);
            return NAN;
        }

        for(int i=0; i<n; i++) {
            if (x_r<x[i]) {
                i--;
                if (Abs(x_r-x[i])<=Abs(x[i+1]-x_r)) {
                    double val= y[i];
                    free(point);
                    return val;
                }
                else {
                    double val= y[i+1];
                    free(point);
                    return val;
                }
            }
            else if (x_r==x[i]) {
                double val= y[i];
                free(point);
                return val;
            }
        }
    }
    else if (!strcmp(s,"makima")) {
        if (n<3) {
            printf("Makima Interpolation not possible!\n");
            printf("At least 3 points required.\n");
            free(point);
            return NAN;
        }

        double *h= malloc(sizeof(double)*(n-1));
        if (!h) {
            free(point);
            return NAN;
        }

        for (int i=0; i<n-1; i++) {
            h[i]=x[i+1]-x[i];
        }

        double *delta= malloc(sizeof(double)*(n+3));
        if (!delta) {
            free(point);
            free(h);
            return NAN;
        }

        double *m= malloc(sizeof(double)*n);
        if (!m) {
            free(point);
            free(h);
            free(delta);
            return NAN;
        }

        int interval=0;
        for (int i=0; i<n; i++) {
            if (x[i]>x_r) {
                interval=i-1;
                break;
            }
            else if (x[i]==x_r) {
                double val=y[i];
                free(point);
                free(h);
                free(delta);
                free(m);
                return val;
            }
        }

        for (int i=2; i<n+1; i++) {
            delta[i]= (y[i-1]-y[i-2])/(x[i-1]-x[i-2]);
        }

        delta[0]= 3*delta[2]-2*delta[3];
        delta[1]= 2*delta[2]-delta[3];
        delta[n+1]= 2*delta[n]-delta[n-1];
        delta[n+2]= 3*delta[n]-2*delta[n-1];

        double c,d;
        for (int i=2; i<n+2; i++) {
            double w1= Abs(delta[i+1]-delta[i])+ Abs(delta[i+1]+delta[i])/2;
            double w2= Abs(delta[i-1]-delta[i-2])+ Abs(delta[i-1]+delta[i-2])/2;

            if (w1+w2!=0) {
                m[i-2]= (w1*delta[i-1]+w2*delta[i])/(w1+w2);
            }
            else {
                m[i-2]= (delta[i-1]+delta[i])/2;
            }
        }

        c= (3*delta[interval+2]-2*m[interval]-m[interval+1])/h[interval];
        d= (m[interval]+m[interval+1]-2*delta[interval+2])/(h[interval]*h[interval]);

        double val= y[interval]+ m[interval]*(x_r-x[interval])+
            c*(x_r-x[interval])*(x_r-x[interval])+
                d*(x_r-x[interval])*(x_r-x[interval])*(x_r-x[interval]);

        //
        free(point);
        free(delta);
        free(h);
        free(m);
        return val;
    }
    else if (!strcmp(s,"polyfit")) {
        if (n<4) {
            printf("Error! At least 4 pints required.\n");
            free(point);
            return NAN;
        }
        double *co_eff= malloc(sizeof(double)*4);
        if (!co_eff) {
            free(point);
            return NAN;
        }

        polyfit_without_file(x,y,co_eff,n,3);

        double val= 0;
        for (int i=0; i<4; i++) {
            val+= co_eff[i]*pow(x_r,i);
        }

        free(point);
        free(co_eff);
        return val;
    }
    else {
        printf("Error! Type a valid method.\n");
        // free
        free(point);
        return NAN;
    }
    free(point);
    return NAN;
}

int main() {
    // FILE *fp= fopen("C:/Users/Rifat/Downloads/50_data_points.txt","r");
    // if (!fp) {
    //     return 1;
    // }
    double x[5]= {0.0,1.0,2.0,3.0,4.0};
    double y[5]= {0.0,1.0,2.0,3.0,4.0};

    char s[15]="polyfit";
    printf("%.6lf\n",interpolate_without_file(x,y,5,6,s));

    // int n=total_point(fp);
    // fclose(fp);
}