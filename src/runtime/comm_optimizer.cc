#include "comm_optimizer.h"
#include <numeric>
#include <random>
#include <iostream>
#include <assert.h>
#include <algorithm>
//#include <bits/stdc++.h>

namespace VirToPhyMapper {
    CommOptimizer::CommOptimizer(idx_t numNPUs) {
        this->numNPUs=numNPUs;
        this->randomMap.resize(numNPUs);
        std::iota (std::begin(randomMap), std::end(randomMap), 0);
        this->seed=0;
        this->commGraphEdges=0;
        options=new idx_t[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_OBJTYPE]=METIS_OBJTYPE_CUT;
        //options[METIS_OPTION_CTYPE]=METIS_CTYPE_SHEM;
        //options[METIS_OPTION_IPTYPE]=METIS_IPTYPE_EDGE;
        //options[METIS_OPTION_RTYPE]=METIS_RTYPE_SEP1SIDED;
        //options[METIS_OPTION_NCUTS]=32;
        //options[METIS_OPTION_NSEPS]=32;
        //options[METIS_OPTION_NITER]=20000;
        //options[METIS_OPTION_MINCONN]=1;
        options[METIS_OPTION_NO2HOP]=0;
        options[METIS_OPTION_NCUTS]=32;
        options[METIS_OPTION_NSEPS]=32;
        options[METIS_OPTION_NITER]=10000;
        
        //options[METIS_OPTION_DBGLVL]=METIS_DBG_TIME;
    }
    int64_t CommOptimizer::getInitialObjective(std::vector <uint32_t> &clusterHierarchy) {
        idx_t nparts=clusterHierarchy[0];
        idx_t NPUsPerNode=numNPUs/nparts;
        int64_t objective=0;
        for(int i=0;i<commGraph.size();i++){
            idx_t nodeID=i/NPUsPerNode;
            idx_t startRange=nodeID*NPUsPerNode;
            idx_t endRange=(nodeID*NPUsPerNode)+NPUsPerNode-1;
            for(auto &a:commGraph[i]){
                if(a.first<startRange || a.first>endRange){
                    objective+=a.second;
                }
            }
        }
        return objective/2;
    }
    CommOptimizer::~CommOptimizer() {
        delete [] options;
    }
    std::vector<uint32_t> & CommOptimizer::mapRandom() {
        std::shuffle(randomMap.begin(), randomMap.end(), std::default_random_engine(seed));
        return randomMap;
    }
    void CommOptimizer::resetCommGraph() {
        commGraph.clear();
        commGraph.resize(numNPUs);
        this->commGraphEdges=0;

    }
    void CommOptimizer::addComm(idx_t src, idx_t dst, idx_t commSize) {
        assert(src>=0 && src<numNPUs && dst>=0 && dst<numNPUs && src!=dst);
        if(commGraph[src].find(dst)==commGraph[src].end()){
            commGraph[src][dst]=commSize;
            commGraph[dst][src]=commSize;
            commGraphEdges++;
        }
        else{
            commGraph[src][dst]+=commSize;
            commGraph[dst][src]+=commSize;
        }
        return;
    }
    void CommOptimizer::convertCommGraphToCSR(idx_t **xadj, idx_t **adjncy, idx_t **adjwgt) {
        *xadj=new idx_t[numNPUs+1];
        *adjncy=new idx_t[commGraphEdges*2];
        *adjwgt=new idx_t[commGraphEdges*2];
        idx_t xadj_pointer=0;
        for(int i=0;i<numNPUs;i++){
            (*xadj)[i]=xadj_pointer;
            for(auto &node:commGraph[i]){
                (*adjncy)[xadj_pointer]=node.first;
                (*adjwgt)[xadj_pointer++]=node.second;
                //std::cout<<node.first<<std::endl;
                //std::cout<<node.second<<std::endl;
            }
            (*xadj)[i+1]=xadj_pointer;
        }
        return;
    }
    std::vector<uint32_t> CommOptimizer::mapByClustering(std::vector<uint32_t> &clusterHierarchy,int64_t &objective) {
        idx_t *xadj,*adjncy,*adjwgt;
        idx_t ncon=1;
        idx_t nparts=clusterHierarchy[0];
        idx_t NPUsPerNode=numNPUs/nparts;
        idx_t edgeCut;
        idx_t *partition=new idx_t[numNPUs];

        convertCommGraphToCSR(&xadj,&adjncy,&adjwgt);
        /*std::cout<<"xadj: ";
        for(int i=0;i<numNPUs+1;i++){
            std::cout<<xadj[i]<<", ";
        }
        std::cout<<std::endl;
        std::cout<<"adjncy: ";
        for(int i=0;i<2*commGraphEdges;i++){
            std::cout<<adjncy[i]<<", ";
        }
        std::cout<<std::endl;
        std::cout<<"adjwgt: ";
        for(int i=0;i<2*commGraphEdges;i++){
            std::cout<<adjwgt[i]<<", ";
        }
        std::cout<<std::endl;*/
        auto result=METIS_PartGraphKway(&numNPUs,&ncon,xadj,adjncy,NULL,NULL,adjwgt,&nparts,NULL,NULL,options,&edgeCut,partition);
        assert(result==METIS_OK);
        /*if(result==METIS_OK){
            std::cout<<"partition done successfully!"<<std::endl;
        }*/
        //std::cout<<"edge cut is:"<<edgeCut<<std::endl;
        /*std::cout<<"partition: ";
        for(int i=0;i<numNPUs;i++){
            std::cout<<partition[i]<<", ";
        }
        std::cout<<std::endl;*/
        objective=edgeCut;
        std::vector<uint32_t> partitionPointer(nparts,0);
        std::vector<uint32_t> mapper;
        mapper.resize(numNPUs);

        for(int i=0;i<numNPUs;i++){
            int part=partition[i];
            mapper[i]=(part*NPUsPerNode)+(partitionPointer[part]);
            partitionPointer[part]++;
            if(partitionPointer[part]>NPUsPerNode){
                delete [] xadj;
                delete [] adjncy;
                delete [] adjwgt;
                delete [] partition;
                mapper.resize(numNPUs);
                std::iota (std::begin(mapper), std::end(mapper), 0);
                objective=-1;
                return mapper;
            }
        }
        /*std::cout<<"mapper: ";
        for(int i=0;i<numNPUs;i++){
            std::cout<<mapper[i]<<", ";
        }
        std::cout<<std::endl;*/
        delete [] xadj;
        delete [] adjncy;
        delete [] adjwgt;
        delete [] partition;
        return mapper;
    }

}
