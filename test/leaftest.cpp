#include "leaf.h"
#include "node.h"
#include <iostream>

void Leaf::mpoleToLocalTest() {
    buildLocalCoeffs();
    auto iList = getInteractionList();

    std::ofstream outFile, outAnlFile;
    outFile.open("out/local.txt");
    outAnlFile.open("out/localAnl.txt");

    for (const auto& obs : getPsn()) {
        cmplx phi;
        for (size_t k = 0; k < P_; ++k)
            phi -= localCoeffs[k] * pow(obs-zk, k);
        outFile << phi << " ";
    }
    outFile << "\n";

    for (const auto& obs : getPsn()) {
        cmplx phi;
        for (const auto& iNode : iList) {
            cmplxVec srcs = iNode->getPsn();
            for (size_t src = 0; src < srcs.size(); ++src)
                phi -= (iNode->getQs())[src] * std::log(obs - srcs[src]);
        }
        outAnlFile << phi << " ";
    }
    outAnlFile << "\n";

    setNodeStat(3);

    for (const auto& iNode : iList)
        iNode->setNodeStat(2);

    std::ofstream psnFile, nodeFile;
    psnFile.open("out/srcs.txt");
    nodeFile.open("out/nodes.txt");

    printPsn(psnFile);
    printNode(nodeFile);
}

/*
void Leaf::iListTest() {
    setNodeStat(3);
    buildInteractionList();
    auto nbors = getInteractionList();

    for (const auto& nbor : nbors)
        nbor->setNodeStat(2);

    std::ofstream psnFile, nodeFile;
    psnFile.open("out/srcs.txt");
    nodeFile.open("out/nodes.txt");

    printPsn(psnFile);
    printNode(nodeFile);
}*/