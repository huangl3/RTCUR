#include "mex.h"
#include "stdio.h"
#include "stdlib.h"
#include "matrix.h"
// r:rlogn; n:n; s:|J|
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6)
        mexErrMsgTxt("A requires Five input arguments");
    else if  (nlhs != 1)
        mexErrMsgTxt("A requires One output arguments");
    double *core = mxGetPr(prhs[0]);
    double *X1 = mxGetPr(prhs[1]);
    double *X2 = mxGetPr(prhs[2]);
    double *X3 = mxGetPr(prhs[3]);
    double *J = mxGetPr(prhs[4]);
    int dim = static_cast<int>(mxGetScalar(prhs[5]));
    size_t n1, n2, n3, d;
    n1 = mxGetM(prhs[1]);
    n2 = mxGetM(prhs[2]);
    n3 = mxGetM(prhs[3]);
    if (dim == 1)
        d = mxGetM(prhs[2]);
    else
        d = mxGetM(prhs[1]);
    size_t s = mxGetN(prhs[4]);
    int** p = new int*[s];
    for(int i = 0; i < s; i++)
        p[i] = new int[2];
    for (int i=0;i<s;i++){
        p[i][0] = (static_cast<int>(J[i])-1) % static_cast<int>(d);
        p[i][1] = (static_cast<int>(J[i])-1) / static_cast<int>(d);
    }
    size_t r = mxGetDimensions(prhs[0])[dim-1];
    plhs[0] = mxCreateDoubleMatrix(n1*s, 1, mxREAL);

    double* C = new double[r*s];
    double *C_out= mxGetPr(plhs[0]);
    mxSetM(plhs[0],n1);
    mxSetN(plhs[0],s);
    // double* C = mxGetPr(plhs[0]);
    // double *C_out= new double[n1*s];
    // mxSetM(plhs[0],r);
    // mxSetN(plhs[0],s);

    double* CX = new double[r*r];
    double temp;
    if (dim == 1){
        for (int i=0;i<s;i++){
//          C(:,i) = ttm(core,{X2(p(i,1),:),X3(p(i,2),:)},[2,3]);
            for (int j=0;j<r;j++){
                for (int l=0;l<r;l++){
                    temp = 0;
                    for (int k=0;k<r;k++){
                        temp += core[j+k*r+l*r*r]*X2[p[i][0]+k*n2];}
                    CX[j+l*r] = temp;
                    // mexPrintf("CX:%f,", temp);  
                }
            }
            for (int j=0;j<r;j++){
                temp = 0;
                for (int k=0;k<r;k++){
                    temp += CX[j+k*r]*X3[p[i][1]+k*n3];
                }
                C[j+i*r] = temp;    
            }
        }

     // C_out = C*X1
        for (int i=0;i<n1;i++){
            for (int j=0;j<s;j++){
                temp = 0;
                for (int k=0;k<r;k++)
                    temp += X1[i+k*n1] * C[k+j*r];
                C_out[i+j*n1] = temp;
            }
        }
    }

    if (dim == 2){
        for (int i=0;i<s;i++){
            for (int k=0;k<r;k++){
                for (int l=0;l<r;l++){
                    temp = 0;
                    for (int j=0;j<r;j++){
                        temp += core[j+k*r+l*r*r]*X1[p[i][0]+j*n1];}
                    CX[k+l*r] = temp;
                    // mexPrintf("CX:%f,", temp);
                }
            }
            for (int j=0;j<r;j++){
                temp = 0;
                for (int k=0;k<r;k++){
                    temp += CX[j+k*r]*X3[p[i][1]+k*n3];
                }
                C[j+i*r] = temp;    
            }
        }
        for (int i=0;i<n2;i++){
            for (int j=0;j<s;j++){
                temp = 0;
                for (int k=0;k<r;k++)
                    temp += X2[i+k*n2] * C[k+j*r];
                C_out[i+j*n2] = temp;
            }
        }
    }

    if (dim == 3){
        for (int i=0;i<s;i++){
            for (int k=0;k<r;k++){
                for (int l=0;l<r;l++){
                    temp = 0;
                    for (int j=0;j<r;j++){
                        temp += core[j+k*r+l*r*r]*X1[p[i][0]+j*n1];}
                    CX[k+l*r] = temp;
                    // mexPrintf("CX:%f,", temp);
                }
            }
            for (int j=0;j<r;j++){
                temp = 0;
                for (int k=0;k<r;k++){
                    temp += CX[j*r+k]*X2[p[i][1]+k*n2];
                    // mexPrintf("CX:%i,", j+k*r); 
                    // mexPrintf("X2:%f,", X2[p[i][1]+k*n2]); 
                }
                C[j+i*r] = temp; 
                // mexPrintf("CX:%f,", temp);   
            }
        }
        for (int i=0;i<n3;i++){
            for (int j=0;j<s;j++){
                temp = 0;
                for (int k=0;k<r;k++)
                    temp += X3[i+k*n3] * C[k+j*r];
                C_out[i+j*n3] = temp;
            }
        }
    }

    delete [] C, CX;
    return;
}