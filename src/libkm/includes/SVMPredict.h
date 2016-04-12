/* 
 * File:   SVMPredict.h
 * Author: veraalva
 *
 * Created on February 26, 2016, 4:31 PM
 */

#ifndef SVMPREDICT_H
#define SVMPREDICT_H

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

        int getPredictProbability() const {
            return predict_probability;
        }

        void setPredictProbability(int predict_probability) {
            this->predict_probability = predict_probability;
        }

        int getNrClass() const {
            return nr_class;
        }

        double* getProbEstimates() {
            return prob_estimates;
        }

        int getSVMType() const {
            return svm_type;
        }

        void svmLoadModel(std::string fileName);
        void svmPredictCalulation(struct svm_node *x, double target_label);
    private:
        struct svm_model* model;
        int predict_probability;
        int svm_type;
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

