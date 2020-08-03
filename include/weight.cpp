#include <iostream>
#include <algorithm>
#include <sstream>
#include <string.h>
#include <mpi.h>

#include "../include/weight.hpp"

float calc_max_weight(framework_t<param_t> **accepted, uint A_per_proc, uint np ) {

    float max_weight=0,max_weights=0;

    for (uint a=0;a<A_per_proc;a++) {
        max_weights=MAX(*(accepted[np*A_per_proc+a]->w),max_weights);
    }
    MPI_Allreduce(&max_weights,&max_weight,1,MPI_FLOAT,MPI_MAX, MPI_COMM_WORLD);

    return max_weight;
}

float calc_sum_weight(framework_t<param_t> **current, uint A_per_proc, uint np) {

    float sum_weight=0, sum_weights=0;
    for (uint a=0;a<A_per_proc;a++) {
        sum_weights+=*(current[np*A_per_proc+a]->w);
    }
    MPI_Allreduce(&sum_weights,&sum_weight,1,MPI_FLOAT,MPI_SUM, MPI_COMM_WORLD);

    return sum_weight;

}

bool sortMethod(framework_t<param_t> *p1,framework_t<param_t> *p2) { //how to sort our frameworks

    return *(p1->d)<*(p2->d);

}
