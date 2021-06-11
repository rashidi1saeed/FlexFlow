#ifndef _COMM_OPTIMIZER_H_
#define _COMM_OPTIMIZER_H_
#include <iostream>
#include <vector>
#include <map>
#include <stdint.h>
#include "metis.h"

namespace VirToPhyMapper {
   class CommOptimizer{
        public:
            idx_t numNPUs;
            std::vector<uint32_t> randomMap;
            unsigned seed;
            std::vector<std::map<idx_t,idx_t>> commGraph;
            uint64_t commGraphEdges;

            CommOptimizer(idx_t numNPUs);
            std::vector<uint32_t>& mapRandom();
            ~CommOptimizer();

            idx_t *options;
            void resetCommGraph();
            void addComm(idx_t src,idx_t dst,idx_t commSize);
            std::vector<uint32_t> mapByClustering(std::vector<uint32_t> &clusterHierarchy,int64_t &objective);
            int64_t getInitialObjective(std::vector<uint32_t> &clusterHierarchy);

        private:
            void convertCommGraphToCSR(idx_t **xadj,idx_t **adjncy,idx_t **adjwgt);


   };
}
#endif