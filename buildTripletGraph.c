// Copyright (C) Yoni Kasten, Weizmann Institute, 2019
#include "mex.h"



void getMatrixIndex(size_t * matrixInd,size_t row,size_t column,size_t numRows,size_t numCols)
{
    *matrixInd =(column-1)*(numRows)+(row-1);
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    int *inputs; double *output;//double *outputB;double *q_2s;double *q_3s;double *q_4s;
    size_t numtriplets;
    int tuplei_1,tuplei_2,tuplei_3;
    int tuplej_1,tuplej_2,tuplej_3;
    int tuplei[3];
    int tuplej[3];
    int numequals;
    int indd1,indd2;
    
    int i,j,k;
    
    /*  create a pointer to the input matrix y */
    inputs =(int *) mxGetPr(prhs[0]);
    
    /*  get the dimensions of the matrix input y */
    numtriplets = mxGetN(prhs[0]);
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix( (mwSize)numtriplets, (mwSize)numtriplets, mxREAL);
    
    /*  create a C pointer to a copy of the output matrix */
    output = mxGetPr(plhs[0]);
    
    for(i=0;i<numtriplets-1;i++)
    {
        tuplei[0]=inputs[i*3];
        tuplei[1]=inputs[i*3+1];
        tuplei[2]=inputs[i*3+2];
        
        
        for(j=i+1;j<numtriplets;j++)
        {
            tuplej[0]=inputs[j*3];
            tuplej[1]=inputs[j*3+1];
            tuplej[2]=inputs[j*3+2];
            numequals=0;
            indd1=0;
            indd2=0;
            
            for (k=0;k<4;k++){
                if(tuplei[indd1]==tuplej[indd2]){
                    numequals=numequals+1;
                    indd1=indd1+1;
                    indd2=indd2+1;
                }
                else if  (tuplei[indd1]>tuplej[indd2]){
                    
                    indd2=indd2+1;
                }else {
                    indd1=indd1+1;
                }
                
                if (indd1>2 || indd2>2){
                    break;
                }
            }
            
            
            if (numequals>=2){
                output[i*numtriplets+j]=1;
                output[j*numtriplets+i]=1;
            }
        }
    }
    
    
}
