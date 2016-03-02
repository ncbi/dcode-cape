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
        SVMPredict(const SVMPredict& orig);
        virtual ~SVMPredict();

        svm_model* GetModel() {
            return model;
        }

        void SetModel(svm_model* model) {
            this->model = model;
        }

        int GetPredict_probability() const {
            return predict_probability;
        }

        void SetPredict_probability(int predict_probability) {
            this->predict_probability = predict_probability;
        }

        int GetNr_class() const {
            return nr_class;
        }

        double* GetProb_estimates() {
            return prob_estimates;
        }

        int GetSvm_type() const {
            return svm_type;
        }

        void SVMLoadModel(char *fileName);
        void SVMPredictCalulation(struct svm_node *x, double target_label);
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

