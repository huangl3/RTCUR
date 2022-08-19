/*
 * cpp version of ttm_reduce (3D only)
 * input _core_ should be converted to 'double' instead of 'tensor'
 */
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::mex;
using namespace matlab::data;

class MexFunction : public Function {
private:
std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;

public:
MexFunction()
{
matlabPtr = getEngine();
}
ArrayFactory factory;
std::ostringstream stream;
void displayOnMATLAB(std::ostringstream stream)
{
ArrayFactory factory;
matlabPtr->feval(u"fprintf", 0, std::vector<Array>
        ({ factory.createScalar(stream.str())}));
}
  
void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
    //checkArguments(outputs, inputs);
    ArrayDimensions dims = inputs[0].getDimensions();
    TypedArray<double> core = std::move(inputs[0]);
    TypedArray<double> X1 = std::move(inputs[1]);
    TypedArray<double> X2 = std::move(inputs[2]);
    TypedArray<double> X3 = std::move(inputs[3]);
    TypedArray<double> J  = std::move(inputs[4]);
    int dim = inputs[5][0];
    size_t d;
    if (dim == 1)
        d = X2.getDimensions()[0];
    else
        d = X1.getDimensions()[0];
    
    size_t s = J.getDimensions()[1];
    int sz = static_cast<int>(s);
    int** p = new int*[s];
    for(int i = 0; i < s; i++)
        p[i] = new int[2];
    index2D(J,d,p);
    size_t C_dim1 = core.getDimensions()[dim-1];
    TypedArray<double> C = factory.createArray<double>({C_dim1,s});
    TypedArray<double> C_out = factory.createArray<double>({d,s});
    int flop_count = 0;
    if (dim == 1){
        for (int i=0;i<sz;i++){
//             C(:,i) = ttm(core,{X2(p(i,1),:),X3(p(i,2),:)},[2,3]);
            for (int j=0;j<C_dim1;j++){
                double temp = 0;
                for (int k=0;k<C_dim1;k++){
                    for (int l=0;l<C_dim1;l++){
                        temp += core[j][k][l]*X2[p[i][0]][k]*X3[p[i][1]][l];
                        flop_count += 2;}
                }
                C[j][i] = temp;
            }
        }
        for (int i=0;i<d;i++){
            for (int j=0;j<s;j++){
                C_out[i][j] = 0;
                for (int k=0;k<C_dim1;k++)
                    C_out[i][j] += C[k][j] * X1[i][k];
                    flop_count += 1;
            }
        }
    }
    
    if (dim == 2){
        for (int i=0;i<sz;i++){
//             C(:,i) = ttm(core,{X1(p(i,1),:),X3(p(i,2),:)},[1,3]);
            for (int j=0;j<C_dim1;j++){
                C[j][i] = 0;
                for (int k=0;k<C_dim1;k++){
                    for (int l=0;l<C_dim1;l++)
                        C[j][i] += core[k][j][l]*X1[p[i][0]][k]*X3[p[i][1]][l];
                }
            }
        }
        for (int i=0;i<d;i++){
            for (int j=0;j<s;j++){
                C_out[i][j] = 0;
                for (int k=0;k<C_dim1;k++)
                    C_out[i][j] += C[k][j] * X2[i][k];
            }
        }
    }
    
    if (dim == 3){
        for (int i=0;i<sz;i++){
//             C(:,i) = ttm(core,{X1(p(i,1),:),X2(p(i,2),:)},[1,2]);
            for (int j=0;j<C_dim1;j++){
                C[j][i] = 0;
                for (int k=0;k<C_dim1;k++){
                    for (int l=0;l<C_dim1;l++)
                        C[j][i] += core[k][l][j]*X1[p[i][0]][k]*X2[p[i][1]][l];
                }
            }
        }
        for (int i=0;i<d;i++){
            for (int j=0;j<s;j++){
                C_out[i][j] = 0;
                for (int k=0;k<C_dim1;k++)
                    C_out[i][j] += C[k][j] * X3[i][k];
            }
        }
    }
    
    
    
    
//     ArrayFactory factory;
//     TypedArray<double> out = factory.createArray<double>({ d, s });
//     for (int i=0;i<d;i++){
//         for (int j=0;j<s;j++)
//             out[i][j] = C_out[i][j];
//     }
    stream<<"Total FLOPS: "<< flop_count <<std::endl;
    displayOnMATLAB(std::move(stream));
    outputs[0] = std::move(C_out);
}

void index2D(const matlab::data::TypedArray<double>& J, const size_t d1, int **p) {
    int sz = static_cast<int>(J.getDimensions()[1]);
    for (int i=0;i<sz;i++){
        p[i][0] = (static_cast<int>(J[0][i])-1) % static_cast<int>(d1);
        p[i][1] = (static_cast<int>(J[0][i])-1) / static_cast<int>(d1);
//         stream<<"p[i][1]: "<<static_cast<int>(J[0][i])-1 / static_cast<int>(d1)<<std::endl;
//     displayOnMATLAB(std::move(stream));
    }
}



};
