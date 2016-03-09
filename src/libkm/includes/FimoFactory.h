/* 
 * File:   FimoFactory.h
 * Author: veraalva
 *
 * Created on March 4, 2016, 3:48 PM
 */

#ifndef FIMOFACTORY_H
#define FIMOFACTORY_H

namespace fimo {

    class Fimo {
    public:

        unsigned long int GetEnd() const {
            return end;
        }

        void SetEnd(unsigned long int end) {
            this->end = end;
        }

        std::string GetId() const {
            return id;
        }

        void SetId(std::string id) {
            this->id = id;
        }

        std::string GetMotif() const {
            return motif;
        }

        void SetMotif(std::string motif) {
            this->motif = motif;
        }

        double GetValue() const {
            return pValue;
        }

        void SetValue(double Value) {
            pValue = Value;
        }

        double GetScore() const {
            return score;
        }

        void SetScore(double score) {
            this->score = score;
        }

        std::string GetSeq() const {
            return seq;
        }

        void SetSeq(std::string seq) {
            this->seq = seq;
        }

        unsigned long int GetStart() const {
            return start;
        }

        void SetStart(unsigned long int start) {
            this->start = start;
        }

        char GetStrand() const {
            return strand;
        }

        void SetStrand(char strand) {
            this->strand = strand;
        }

        double GetExpression() const {
            return expression;
        }

        void SetExpression(double expression) {
            this->expression = expression;
        }
        
        std::string GetExpEnsembl() const {
            return expEnsembl;
        }

        void SetExpEnsembl(std::string expEnsembl) {
            this->expEnsembl = expEnsembl;
        }

        bool operator==(const Fimo& right) const {
            return this->GetMotif().compare(right.GetMotif()) == 0 &&
                    this->start == right.GetStart() &&
                    this->end == right.GetEnd();
        }

        bool operator!=(const Fimo& right) const {
            bool result = !(*this == right);
            return result;
        }

        bool operator>(const Fimo& right) const {
            if (this->start == right.GetStart()) {
                if (this->end == right.GetEnd()) {
                    return this->expression > right.GetExpression();
                }
                return this->end > right.GetEnd();
            }
            return this->start > right.GetStart();
        }

        bool operator<(const Fimo& right) const {
            return right > * this;
        }

        friend std::ostream& operator<<(std::ostream& os, const Fimo& obj) {
            os << obj.GetId() << "\t" << obj.GetMotif() << "\t" << obj.GetStart() << "\t" << obj.GetEnd()
                    << "\t" << obj.GetValue() << "\t" << obj.GetExpression() << "\t" << obj.GetExpEnsembl();
            return os;
        }

    private:
        std::string motif;
        std::string id;
        unsigned long int start;
        unsigned long int end;
        char strand;
        double score;
        double pValue;
        double expression;
        std::string expEnsembl;
        std::string seq;
    };

    class FimoFactory {
    public:
        FimoFactory();
        FimoFactory(const FimoFactory& orig);
        virtual ~FimoFactory();

        void CreateTissueIndexFromFiles(char *pwm_EnsembleID, char *tissue_file);
        void CreateCutoffIndexFromFile(char *cutoffFileName, size_t column);
        void ParseFimoOutput(char *fimoOuputName, char *tissueCode, unsigned long int snpPos);

        std::pair<std::string, double> GetTissueValue(std::string motifName, std::string tissueName) {
            std::unordered_map<std::string, std::unordered_map<std::string, std::pair < std::string, double>>>::iterator it = tissueIndex.find(motifName);
            if (it != tissueIndex.end()) {
                std::unordered_map<std::string, std::pair < std::string, double>>::iterator it1 = it->second.find(tissueName);
                if (it1 != it->second.end()) return it1->second;
            }
            return std::pair<std::string, double>("", 0.0000);
        }

        double GetCutoffValue(std::string motifName) {
            std::unordered_map<std::string, double>::iterator it = cutoffIndex.find(motifName);
            if (it != cutoffIndex.end()) {
                return it->second;
            }
            return -1.0;
        }
        
        std::unordered_map<std::string, std::vector<double> >& GetSnpIDMap() {
            return snpIDMap;
        }
    private:
        std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, double>>> tissueIndex;
        std::unordered_map<std::string, double> cutoffIndex;

        /**
         * This map is used to store the SNP id with a vector of two elements.
         * The first element is the sum of the co-occupied and the second the 
         * sum of the neighbors 
         */
        std::unordered_map<std::string, std::vector<double>> snpIDMap;

        struct PointerCompare {

            bool operator()(const Fimo* l, const Fimo* r) {
                return *l != *r;
            }
        };
    };
}

#endif /* FIMOFACTORY_H */

