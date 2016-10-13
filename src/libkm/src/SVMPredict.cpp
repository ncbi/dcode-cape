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
#include <memory>
#include <stdbool.h>
#include <set>

#include "berror.h"
#include "bmemory.h"
#include "Global.h"
#include "TimeUtils.h"
#include "svm.h"
#include "SVMPredict.h"

using namespace std;
using namespace svm;

SVMPredict::SVMPredict() {
    this->correct = 0;
    this->total = 0;
    this->error = 0;
    this->sump = 0;
    this->sumt = 0;
    this->sumpp = 0;
    this->sumtt = 0;
    this->sumpt = 0;
    this->model = NULL;
    this->nr_class = 0;
    this->prob_estimates = NULL;
    this->prob = (svm_problem *) allocate(sizeof (svm_problem), __FILE__, __LINE__);

    this->param.svm_type = C_SVC;
    this->param.kernel_type = LINEAR;
    this->param.degree = 3;
    this->param.gamma = 0; // 1/num_features
    this->param.coef0 = 0;
    this->param.nu = 0.5;
    this->param.cache_size = 100;
    this->param.C = 1;
    this->param.eps = 1e-3;
    this->param.p = 0.1;
    this->param.shrinking = 1;
    this->param.probability = 0;
    this->param.nr_weight = 2;
    this->param.weight_label = (int *) allocate(sizeof (int) * param.nr_weight, __FILE__, __LINE__);
    this->param.weight_label[0] = 1;
    this->param.weight_label[1] = -1;
    this->param.weight = (double *) allocate(sizeof (double) * param.nr_weight, __FILE__, __LINE__);
    this->param.weight[0] = 4.0;
    this->param.weight[1] = 1.0;
}

SVMPredict::~SVMPredict() {
    if (this->model) svm_free_and_destroy_model(&this->model);
    if (this->prob) free(this->prob);
    if (this->param.weight_label) free(this->param.weight_label);
    if (this->param.weight) free(this->param.weight);
}

void SVMPredict::svmLoadModel(std::string fileName) {
    if ((this->model = svm_load_model(fileName.c_str())) == 0) {
        cerr << "Can't open model file " << fileName << endl;
        exit(1);
    }
    if (this->getPredictProbability()) {
        if (svm_check_probability_model(this->model) == 0) {
            cerr << "Model does not support probability estimates" << endl;
            exit(1);
        }
    } else {
        if (svm_check_probability_model(this->model) != 0)
            cout << "Model supports probability estimates, but disabled in prediction.\n";
    }

    this->param.svm_type = svm_get_svm_type(this->model);
    this->nr_class = svm_get_nr_class(this->model);
    this->prob_estimates = NULL;

    if (this->getPredictProbability()) {
        if (this->param.svm_type == NU_SVR || this->param.svm_type == EPSILON_SVR)
            cout << "Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma=" << svm_get_svr_probability(model) << endl;
        else {
            this->prob_estimates = (double *) allocate(nr_class * sizeof (double), __FILE__, __LINE__);
        }
    }
}

void SVMPredict::svmPredictCalulation(struct svm_node *x, double target_label) {
    double predict_label;

    if (this->getPredictProbability() && (this->param.svm_type == C_SVC || this->param.svm_type == NU_SVC)) {
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

void SVMPredict::svmTrainModel(int l, double *y, struct svm_node **x, const char *model_file_name) {
    prob->l = l;
    prob->x = x;
    prob->y = y;
    model = svm_train(prob, &param);
    if (svm_save_model(model_file_name, model)) {
        fprintf(stderr, "can't save model to file %s\n", model_file_name);
        exit(1);
    }
    svm_free_and_destroy_model(&model);
}
