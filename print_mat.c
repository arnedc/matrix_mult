#include <stdio.h>
#include "matrix.h"

void printdense(int n, dense_matrix mat, char *filename)
{
  FILE *fd;
  fd = fopen(filename,"w");
  if (fd==NULL)
    printf("error creating file");
  int i,j;
  for (i=0;i<n;++i) {
    fprintf(fd,"[\t");
    for (j=0;j<n;++j) {
      fprintf(fd,"%7.3f\t",*mat++);
    }
    fprintf(fd,"]\n");
  }
  fclose(fd);
}

void printsparse_CRS(int n, struct sparse_CRS *mat, char *filename)
{
  FILE *fd;
  double *v;
  int *c,min,max;
  fd = fopen(filename,"w");
  if (fd==NULL)
    printf("error creating file");
  int i,j;
  for (i=0;i<n;++i) {
    fprintf(fd,"[\t");
    for (j=0;j<n;++j) {
      min=*(mat->rowptr+i);
      max=*(mat->rowptr+i+1);
      c=mat->colnumber+min;
      v=mat->value+min;
      if (j== *c && c<mat->colnumber+max){
	fprintf(fd,"%7.3f\t",*v);
	v++;
	c++;
      }
      else
	fprintf(fd,"%7.3f\t",0.0);
    }
    fprintf(fd,"]\n");
  }
  fclose(fd);
}

void printsparse(int n, sparse_matrix mat, char *filename)
{
  FILE *fd;
  struct row huidig;
  fd = fopen(filename,"w");
  if (fd==NULL)
    printf("error creating file");
  int i,j;
  for (i=0;i<n;++i) {
    fprintf(fd,"[\t");
    huidig=*(mat+i);
    for (j=0;j<n;++j) {
      if (j== *huidig.colnumber){
	fprintf(fd,"%7.3f\t",*huidig.value);
	++(huidig.value);
	++(huidig.colnumber);
      }
      else
	fprintf(fd,"%7.3f\t",0.0);
    }
    fprintf(fd,"]\n");
  }
  fclose(fd);
}

void printsparse2(int n, sparse2_matrix mat, char *filename)
{
  FILE *fd;
  struct cell *huidig;
  fd = fopen(filename,"w");
  if (fd==NULL)
    printf("error creating file");
  int i,j;
  for (i=0;i<n;++i) {
    fprintf(fd,"[\t");
    huidig=(mat+i)->head->next;
    for (j=0;j<n;++j) {
      if (huidig != NULL && j== huidig->colnumber) {
	fprintf(fd,"%7.3f\t",huidig->value);
	huidig=huidig->next;
      }
      else
	fprintf(fd,"%7.3f\t",0.0);
    }
    fprintf(fd,"]\n");
  }
  fclose(fd);
}