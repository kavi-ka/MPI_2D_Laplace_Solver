#include <stdio.h>
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#define n 32
#define arrSize (n/2)+1
// matrix is n x n; only works for 4 processes

void printMatrix(double arr[arrSize][arrSize],int ghostIndexR,int ghostIndexC,bool ghostUpdated ){
  for(int i=0;i<=n/2;i++){
    for(int j=0;j<=n/2;j++){
      if(i!=ghostIndexR && j!=ghostIndexC){
        std::cout<<arr[i][j]<<"    ";
      }
      if(ghostUpdated && (i==ghostIndexR || j== ghostIndexC)){
        std::cout<<"["<<arr[i][j]<<"] "; //prints ghost elements in square brackets
      }

    }
    std::cout<<"\n";
  }

}
int rowPartner(int rank){
  if(rank<2){
    return rank+2;
  }
  return rank-2;
}
int columnPartner(int rank){
  if(rank%2==0){
    return rank+1;
  }
  return rank-1;
}
int rowToSend(int rank){
  if(rank<2){
    return (n/2)-1;
  }
  return 1;
}
int rowToReceive(int rank){
  if(rank<2){
    return (n/2);
  }
  return 0;

}
int columnToSend(int rank){
  if(rank%2==0){
    return (n/2)-1;
  }
  return 1;
}
int columnToRecieve(int rank){
  if(rank%2==0){
    return n/2;
  }
  return 0;
}

int main(int argc, char** argv){
  MPI_Request requests[4];
  double up,down,left,right; //neighbors to each element; 5 point stencil
  int ghostRow, ghostColumn; //stores index of column & row which need to be received by neighboring processes
  double localDiff = 0.0,diff=0.0,totDiff=0.0; // diff is global,
  double temp[n/2],tempRec[n/2]; //temp arrays used for sends and receives respectively
  int rankOfProcess,processNum;
  long iter = 0;
  double localArr[arrSize][arrSize]; //each process has an n/2 x n/2 matrix, with additional ghost row and column
  double newArr[arrSize][arrSize]; //used to calculate difference between iterations

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &processNum);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfProcess);

  //process 0 and 1 need their last row filled in, process 2 and 3 need their first row filled in
  ghostRow = (rankOfProcess < 2)? n/2 : 0;
  //process 0 and 2 need their last column filled in, process 1 and 3 need their first column filled in
  ghostColumn = (rankOfProcess %2 == 0)? n/2 : 0;

  //fills in each process's array with double from 0-49 inclusive
  for(int i=0;i<=n/2;i++){
    for(int j=0;j<=n/2;j++){
      if(i!=ghostRow && j!=ghostColumn){
        localArr[i][j]= rand() % 50;
      }
    }
  }

  //prints out original matrix for process 0
  if(rankOfProcess == 0){
    std::cout<<"Original Matrix: "<<"\n";
    printMatrix(localArr,ghostRow,ghostColumn,false);
  }

  do{ //runs until difference is less than 0.1
    localDiff = 0.0;
    totDiff =0.0;
    iter=iter+1;
    //sends and receives between processes
    MPI_Irecv(localArr[rowToReceive(rankOfProcess)],n/2,MPI_DOUBLE,rowPartner(rankOfProcess),0,MPI_COMM_WORLD,&requests[0]);
    MPI_Isend(localArr[rowToSend(rankOfProcess)],n/2,MPI_DOUBLE,rowPartner(rankOfProcess),0,MPI_COMM_WORLD,&requests[1]);
    MPI_Irecv(&tempRec,n/2,MPI_DOUBLE,columnPartner(rankOfProcess),0,MPI_COMM_WORLD,&requests[2]);
    for(int i=0;i<n/2;i++){
      localArr[i][columnToRecieve(rankOfProcess)]=tempRec[i];
    }
    for(int i=0;i<n/2;i++){
      temp[i]=localArr[i][columnToSend(rankOfProcess)];
    }
    MPI_Isend(&temp,n/2,MPI_DOUBLE,columnPartner(rankOfProcess),0,MPI_COMM_WORLD,&requests[3]);
    MPI_Waitall(4,requests,MPI_STATUSES_IGNORE);

    for(int i=0;i<=n/2;i++){
      for(int j=0;j<=n/2;j++){
        if(i!=ghostRow && j!=ghostColumn){
          //checks if the neighbor in the specified direction exists
          //if neighbor d.n.e, it considers it to be 25.0
          up = (i-1<0) ? 25.0:(localArr[i-1][j]);
          down = (i+1>(n/2))? 25.0:(localArr[i+1][j]);
          left = (j-1<0)? 25.0:(localArr[i][j-1]);
          right = (j+1>(n/2))? 25.0:(localArr[i][j+1]);
          newArr[i][j]=(left+right+up+down)/4.0; //calculates average of neighbors
          totDiff = totDiff + fabs((newArr[i][j]-localArr[i][j])); //calculates difference for each element
          localArr[i][j]=newArr[i][j]; //assigns calculated value to local array
        }
      }
    }

    localDiff = totDiff / ((double)((n/2-1)*(n/2-1))); //calculates average difference across all elements
    MPI_Allreduce(&localDiff,&diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); //combines differences across processes
    if(rankOfProcess==0){ //prints iteration #,difference, and matrix for process 0
      std::cout<<"Iter: "<<iter<<", Diff: "<<diff<<"\n";
      printMatrix(localArr,ghostRow,ghostColumn,true);
    }
  }while(diff>0.01);
  //barrier  to sync processes
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
