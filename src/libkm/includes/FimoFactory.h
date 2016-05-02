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

        Fimo() {
            this->end = 0;
            this->expression = 0.0;
            this->pValue = 0.0;
            this->score = 0.0;
            this->start = 0;
            this->strand = 0;
        }

        virtual ~Fimo() {
        }

        unsigned long int getEnd() const {
            return end;
        }

        void setEnd(unsigned long int end) {
            this->end = end;
        }

        std::string getId() const {
            return id;
        }

        void setId(std::string id) {
            this->id = id;
        }

        std::string getMotif() const {
            return motif;
        }

        void setMotif(std::string motif) {
            this->motif = motif;
        }

        double getValue() const {
            return pValue;
        }

        void setValue(double Value) {
            pValue = Value;
        }

        double getScore() const {
            return score;
        }

        void setScore(double score) {
            this->score = score;
        }

        std::string getSeq() const {
            return seq;
        }

        void setSeq(std::string seq) {
            this->seq = seq;
        }

        unsigned long int getStart() const {
            return start;
        }

        void setStart(unsigned long int start) {
            this->start = start;
        }

        char getStrand() const {
            return strand;
        }

        void setStrand(char strand) {
            this->strand = strand;
        }

        double getExpression() const {
            return expression;
        }

        void setExpression(double expression) {
            this->expression = expression;
        }

        std::string getExpEnsembl() const {
            return expEnsembl;
        }

        void setExpEnsembl(std::string expEnsembl) {
            this->expEnsembl = expEnsembl;
        }

        bool operator==(const Fimo& right) const {
            return this->getMotif().compare(right.getMotif()) == 0 &&
                    this->start == right.getStart() &&
                    this->end == right.getEnd();
        }

        bool operator!=(const Fimo& right) const {
            bool result = !(*this == right);
            return result;
        }

        bool operator>(const Fimo& right) const {
            if (this->start == right.getStart()) {
                if (this->end == right.getEnd()) {
                    return this->expression > right.getExpression();
                }
                return this->end > right.getEnd();
            }
            return this->start > right.getStart();
        }

        bool operator<(const Fimo& right) const {
            return right > * this;
        }

        friend std::ostream& operator<<(std::ostream& os, const Fimo& obj) {
            os << obj.getId() << "\t" << obj.getMotif() << "\t" << obj.getStart() << "\t" << obj.getEnd()
                    << "\t" << obj.getValue() << "\t" << obj.getExpression() << "\t" << obj.getExpEnsembl();
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

        void createTissueIndexFromFiles(std::string pwm_EnsembleID, std::string tissue_file);
        void createCutoffIndexFromFile(std::string cutoffFileName, size_t column);
        void parseFimoOutput(std::string fimoOuputName, std::string tissueCode, unsigned long int snpPos);

        std::pair<std::string, double> getTissueValue(std::string motifName, std::string tissueName) {
            std::unordered_map<std::string, std::unordered_map<std::string, std::pair < std::string, double>>>::iterator it = tissueIndex.find(motifName);
            if (it != tissueIndex.end()) {
                std::unordered_map<std::string, std::pair < std::string, double>>::iterator it1 = it->second.find(tissueName);
                if (it1 != it->second.end()) return it1->second;
            }
            return std::make_pair("", 0.0000);
        }

        double getCutoffValue(std::string motifName) {
            std::unordered_map<std::string, double>::iterator it = cutoffIndex.find(motifName);
            if (it != cutoffIndex.end()) {
                return it->second;
            }
            return -1.0;
        }

        std::unordered_map<std::string, std::vector<double> >& getSnpIDContainer() {
            return snpIDContainer;
        }
    private:
        std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, double>>> tissueIndex;
        std::unordered_map<std::string, double> cutoffIndex;

        /**
         * This map is used to store the SNP id with a vector of two elements.
         * The first element is the sum of the co-occupied and the second the 
         * sum of the neighbors 
         */
        std::unordered_map<std::string, std::vector<double>> snpIDContainer;

        struct PointerCompare {

            bool operator()(const std::shared_ptr<Fimo> l, const std::shared_ptr<Fimo> r) {
                return *l != *r;
            }
        };
    };
}

#endif /* FIMOFACTORY_H */

