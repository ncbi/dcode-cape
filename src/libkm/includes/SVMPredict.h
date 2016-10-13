/* 
 * File:   SVMPredict.h
 * Author: veraalva
 *
 * Created on February 26, 2016, 4:31 PM
 */

#ifndef SVMPREDICT_H
#define SVMPREDICT_H

#include "svm.h"


namespace svm {

    class SVMPredict {
    public:
        SVMPredict();
        virtual ~SVMPredict();

        svm_model* getModel() {
            return model;
        }

        void setModel(svm_model* model) {
            this->model = model;
        }

        svm_problem *getProb() {
            return prob;
        }

        void setSVMProbY(double *y) {
            this->prob->y = y;
        }

        void setSVMProbX(struct svm_node **x) {
            this->prob->x = x;
        }

        int getPredictProbability() const {
            return this->param.probability;
        }

        void setPredictProbability(int predict_probability) {
            this->param.probability = predict_probability;

        }

        int getNrClass() const {
            return nr_class;
        }

        double* getProbEstimates() {
            return prob_estimates;
        }

        int getSVMType() const {
            return this->param.svm_type;
        }

        void setSVMType(int svm_type) {
            this->param.svm_type = svm_type;
        }

        void setSVMKernelType(int kernel_type) {
            this->param.kernel_type = kernel_type;
        }

        void setSVMw1(double w1) {
            this->param.weight[0] = w1;
        }

        void setSVMwMinus1(double wminus1) {
            this->param.weight[1] = wminus1;
        }

        void svmLoadModel(std::string fileName);
        void svmPredictCalulation(struct svm_node *x, double target_label);
        void svmTrainModel(int l, double *y, struct svm_node **x, const char *model_file_name);
    private:
        struct svm_model *model;
        struct svm_parameter param;
        struct svm_problem *prob;
        int nr_class;
        double *prob_estimates;

        int correct;
        int total;
        double error;
        double sump;
        double sumt;
        double sumpp;
        double sumtt;
        double sumpt;

    };
}

#endif /* SVMPREDICT_H */

