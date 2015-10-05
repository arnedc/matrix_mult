#include <stdio.h>
#include <mkl.h>
#include <mkl_cblas.h>
#include <stdlib.h>
#include <stddef.h>
#include "matrix.h"



void multiply_matrix(int n, dense_matrix matA, dense_matrix matB, int lld, dense_matrix matprod)
{
    int i,j,k;
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j) {
            for (k=0; k<n; ++k) {
                *(matprod+j+i*lld) += *(matA + i*lld + k) * *(matB + j + k*lld);
            }
        }
    }
    return;
}

void multiply_matrix_blocked(int n, dense_matrix matA, dense_matrix matB, int blocksize, int lld, dense_matrix matprod)
{
    int i,j,k;
    int blocks;
    if (n%blocksize!=0){
      printf("blocksize not correct");
      return;
    }
    blocks= n/blocksize;
    for (i=0; i<blocks; ++i) {
        for (j=0; j<blocks; ++j) {
            for (k=0; k<blocks; ++k) {
                multiply_matrix(blocksize,matA+i*blocksize*lld+k*blocksize,matB+j*blocksize+k*blocksize*lld,n,matprod+j*blocksize+i*blocksize*lld);
            }
        }
    }
    return;
}

void multiply_blas_matrix(int n, dense_matrix matA, dense_matrix matB, dense_matrix matprod)
{
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1,matA,n,matB,n,0,matprod,n);
    return;
}

struct sparse_CRS *multiply_CRS_dense(int n, struct sparse_CRS *matA, struct sparse_CRS *matB)
{
    int i,j,Acol,pcol,k,count;
    double Aval,Bval,*prow;
    prow=calloc(n,sizeof(double));
    if (prow==NULL) {
        printf("Unable to create dense row (function: multiply_CRS_dense)");
        return NULL;
    }
    for (i=0; i<n; ++i)
        *(prow+i)=0;
    struct sparse_CRS *matprod;
    matprod=(struct sparse_CRS *) malloc(sizeof(struct sparse_CRS));
    if (matprod==NULL) {
        printf("Unable to create sparse matrix (function: multiply_CRS_dense)");
        return NULL;
    }
    matprod->rowptr=(int *) calloc(n+1,sizeof(int));
    if (matprod->rowptr==NULL) {
        printf("Unable to create sparse matrix (function: multiply_CRS_dense)");
        return NULL;
    }
    *(matprod->rowptr)=0;
    matprod->colnumber=malloc(0);
    matprod->value=malloc(0);
    for(i=0; i<n; ++i) {
        *(matprod->rowptr+i+1)=*(matprod->rowptr+i);
        for(j=*(matA->rowptr+i); j<*(matA->rowptr+i+1); ++j) {
            Acol=*(matA->colnumber+j);
            Aval=*(matA->value+j);
            for(k=*(matB->rowptr+Acol); k<*(matB->rowptr+Acol+1); ++k) {
                Bval=*(matB->value+k);
                pcol=*(matB->colnumber+k);
                *(prow + pcol) += Aval * Bval;
            }
        }
        for (j=0; j<n; ++j) {
            if (*(prow+j)!=0) {
                ++(*(matprod->rowptr+i+1));
                count=*(matprod->rowptr+i+1);
                matprod->colnumber=realloc(matprod->colnumber,count * sizeof(int));
                if (matprod->colnumber==NULL) {
                    printf("Unable to create sparse matrix (function: multiply_CRS_dense)");
                    return NULL;
                }
                matprod->value=realloc(matprod->value,count * sizeof(double));
                if (matprod->value==NULL) {
                    printf("Unable to create sparse matrix (function: multiply_CRS_dense)");
                    return NULL;
                }
                *(matprod->colnumber+count)=j;
                *(matprod->value+count)=*(prow+j);
                *(prow+j)=0;
            }
        }
    }
    free(prow);
    return matprod;
}

sparse_matrix multiply_sparse_matrix(int n, sparse_matrix matA, sparse_matrix matB)
{
    int i,j,k,pcol,count, *A_col, *B_col;
    double *prow, *A_val, *B_val;
    struct row* matp;
    matp = (struct row *) calloc(n+1, sizeof(struct row));
    if (matp==NULL) {
        printf("Unable to create sparse matrix (function: multiply_sparse_matrix)");
        return NULL;
    }
    prow=(double *) calloc(n,sizeof(double));
    if (prow==NULL) {
        printf("Unable to create dense row (function: multiply_sparse_matrix)");
        return NULL;
    }
    for (k=0; k<n; ++k)
        prow[k]=0;
    for (i=0; i<n; ++i) {
        A_col=(matA+i)->colnumber;
        A_val=(matA+i)->value;
        count=0;
        for (j=*A_col; j>=0; j=*(++A_col),++A_val) {
            B_col=(matB+j)->colnumber;
            B_val=(matB+j)->value;
            while((pcol=*B_col) != -1) {
                if (prow[pcol]==0)
                    ++count;
                prow[pcol]+= *A_val * *B_val;
                ++B_col;
                ++B_val;
            }
        }
        (matp+i)->colnumber=(int *) calloc(count+1,sizeof(int));
        if ((matp+i)->colnumber==NULL) {
            printf("Unable to create sparse matrix (function: multiply_sparse_matrix)");
            return NULL;
        }
        (matp+i)->value=(double *) calloc(count+1,sizeof(double));
        if ((matp+i)->value==NULL) {
            printf("Unable to create sparse matrix (function: multiply_sparse_matrix)");
            return NULL;
        }
        count=0;
        for (k=0; k<n; ++k) {
            if (prow[k]!=0) {
                *((matp+i)->colnumber+count)=k;
                *((matp+i)->value+count)=prow[k];
                prow[k]=0;
                ++count;
            }
        }
        *((matp+i)->colnumber+count)=-1;
    }
    free(prow);
    return matp;
}

sparse2_matrix multiply_sparse2_matrix(int n, sparse2_matrix matA, sparse2_matrix matB)
{
    int i,pcol;
    double prod;
    struct cell *huidigA, *huidigB, *huidigp, *prev;
    struct cellrow *matp;
    matp = (struct cellrow *) calloc(n+1,sizeof(struct cellrow));
    if (matp==NULL) {
        printf("Unable to create sparse matrix (function: multiply_sparse2_matrix)");
        return NULL;
    }
    for (i=0; i<n; ++i) {
        huidigA= (matA+i)->head->next;
        (matp+i)->head= malloc(sizeof(struct cell));
        if ((matp+i)->head==NULL) {
            printf("Unable to create sparse matrix (function: multiply_sparse2_matrix)");
            return NULL;
        }
        huidigp=(matp+i)->head;
        huidigp->next=NULL;
        huidigp->colnumber=-1;
        while (huidigA != NULL) {
            huidigB=(matB + huidigA->colnumber)->head->next;
            while(huidigB !=NULL) {
                pcol=huidigB->colnumber;
                prod=huidigA->value * huidigB->value;
                for(prev=(matp+i)->head; huidigp != NULL && huidigp->colnumber < pcol; huidigp=huidigp->next)
                    prev=huidigp;
                if (huidigp != NULL && huidigp->colnumber == pcol)
                    huidigp->value +=prod;
                else if(huidigp == NULL) {
                    huidigp= malloc(sizeof(struct cell));
                    if (huidigp==NULL) {
                        printf("Unable to create sparse matrix (function: multiply_sparse2_matrix)");
                        return NULL;
                    }
                    huidigp->colnumber=pcol;
                    huidigp->value=prod;
                    huidigp->next=NULL;
                    prev->next=huidigp;
                }
                else {
                    prev->next=malloc(sizeof (struct cell));
                    if (prev->next==NULL) {
                        printf("Unable to create sparse matrix (function: multiply_sparse2_matrix)");
                        return NULL;
                    }
                    prev->next->colnumber=pcol;
                    prev->next->value=prod;
                    prev->next->next=huidigp;
                }
                huidigB=huidigB->next;
                huidigp=(matp+i)->head;
            }
            huidigA=huidigA->next;
        }
    }
    return matp;
}

sparse2_matrix multiply_sparse2B_matrix(int n, sparse2_matrix matA, sparse2_matrix matB)
{
    int i,j,k,pcol;
    double prod, *prow;
    struct cell *huidigA, *huidigB, *huidigp;
    struct cellrow *matp;
    matp = (struct cellrow *) calloc(n+1,sizeof(struct cellrow));
    if (matp==NULL) {
        printf("Unable to create sparse matrix (function: multiply_sparse2B_matrix)");
        return NULL;
    }
    prow=(double *) calloc(n,sizeof(double));
    if (prow==NULL) {
        printf("Unable to create dense row (function: multiply_sparse2B_matrix)");
        return NULL;
    }
    for (k=0; k<n; ++k)
        prow[k]=0;
    for (i=0; i<n; ++i) {
        huidigA= (matA+i)->head->next;
        while (huidigA != NULL) {
            huidigB=(matB + huidigA->colnumber)->head->next;
            while(huidigB !=NULL) {
                pcol=huidigB->colnumber;
                prod=huidigA->value * huidigB->value;
                prow[pcol]+=prod;
                huidigB=huidigB->next;
            }
            huidigA=huidigA->next;
        }
        (matp+i)->head=malloc(sizeof(struct cell));
        if ((matp+i)->head==NULL) {
            printf("Unable to create sparse matrix (function: multiply_sparse2B_matrix)");
            return NULL;
        }
        (matp+i)->head->colnumber=-1;
        (matp+i)->head->next=NULL;
        huidigp=(matp+i)->head;
        for (j=0; j<n; ++j) {
            if (prow[j]!=0) {
                huidigp->next=malloc(sizeof(struct cell));
                if (huidigp->next==NULL) {
                    printf("Unable to create sparse matrix (function: multiply_sparse2B_matrix)");
                    return NULL;
                }
                huidigp=huidigp->next;
                huidigp->next=NULL;
                huidigp->colnumber=j;
                huidigp->value=prow[j];
                prow[j]=0;
            }
        }
    }
    free(prow);
    return matp;
}
