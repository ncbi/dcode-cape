/* 
 * File:   SVMPredict.cpp
 * Author: veraalva
 * 
 * Created on February 26, 2016, 4:31 PM
 */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <iostream>
#include <stdbool.h>
#include "berror.h"
#include "bmemory.h"
#include "Global.h"
#include "TimeUtils.h"
#include "svm.h"
#include "SVMPredict.h"

using namespace std;
using namespace svm;

SVMPredict::SVMPredict() {
    this->predict_probability = 0;
    this->correct = 0;
    this->total = 0;
    this->error = 0;
    this->sump = 0;
    this->sumt = 0;
    this->sumpp = 0;
    this->sumtt = 0;
    this->sumpt = 0;
    this->model = NULL;
}

SVMPredict::SVMPredict(const SVMPredict& orig) {
}

SVMPredict::~SVMPredict() {
    if (this->model) svm_free_and_destroy_model(&this->model);
}

void SVMPredict::SVMLoadModel(char* fileName) {
    if ((this->model = svm_load_model(fileName)) == 0) {
        fprintf(stderr, "can't open model file %s\n", fileName);
        exit(1);
    }
    if (this->predict_probability) {
        if (svm_check_probability_model(this->model) == 0) {
            fprintf(stderr, "Model does not support probability estimates\n");
            exit(1);
        }
    } else {
        if (svm_check_probability_model(this->model) != 0)
            cout << "Model supports probability estimates, but disabled in prediction.\n";
    }

    this->svm_type = svm_get_svm_type(this->model);
    this->nr_class = svm_get_nr_class(this->model);
    this->prob_estimates = NULL;

    if (this->predict_probability) {
        if (this->svm_type == NU_SVR || this->svm_type == EPSILON_SVR)
            cout << "Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=" << svm_get_svr_probability(model) << endl;
        else {
            this->prob_estimates = (double *) allocate(nr_class * sizeof (double), __FILE__, __LINE__);
        }
    }
}

void SVMPredict::SVMPredictCalulation(struct svm_node *x, double target_label) {
    double predict_label;

    if (this->predict_probability && (this->svm_type == C_SVC || this->svm_type == NU_SVC)) {
        predict_label = svm_predict_probability(this->model, x, this->prob_estimates);
    } else {
        predict_label = svm_predict(this->model, x);
    }

    if (predict_label == target_label)
        ++this->correct;
    this->error += (predict_label - target_label)*(predict_label - target_label);
    this->sump += predict_label;
    this->sumt += target_label;
    this->sumpp += predict_label*predict_label;
    this->sumtt += target_label*target_label;
    this->sumpt += predict_label*target_label;
    ++this->total;
}

