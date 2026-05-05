// #include "em_ctx.h"
#include <iostream>
#include <fstream>
#include "PetscAMGInterface.h"

using namespace std;

int main(int argc, char **argv){
    

    PetscCall(PetscInitialize(NULL, NULL, NULL, NULL));   //初始化

    solve_eg1();

    PetscCall(PetscFinalize());

    return 0;
}