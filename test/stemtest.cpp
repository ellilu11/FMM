#include <cassert>
#include <iostream>
#include "leaf.h"
#include "math.h"
#include "stem.h"

void Stem::mpoleToLocalTest() {

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<> branchIdx(0, 3);
    std::shared_ptr<Node> node = std::make_shared<Stem>(*this);

    // while (node->isNodeType<Stem>())
    while (node->getLvl() > 2)
        node = (node->getBranches())[branchIdx(gen)];

    std::ofstream outFile, outAnlFile;
    outFile.open("out/local.txt");
    outAnlFile.open("out/localAnl.txt");

    const int P = Node::getP();
    for (int p = 1; p <= P; ++p) {
        Node::setP(p);
        Node::buildBinomTable();

        std::cout << " Computing upward pass...   (" << " p = " << p << " )" << std::endl;
        buildMpoleCoeffs();

        std::cout << " Computing downward pass... (" << " p = " << p << " )" << std::endl;
        buildLocalCoeffs();

        auto localCoeffs = node->getLocalCoeffs();

        // add contribution from node's base's interaction list
        localCoeffs += node->shiftBaseLocalCoeffs(); 

        for (const auto& obs : node->getPsn()) {
            cmplx phi;
            for (size_t k = 0; k <= p; ++k)
                phi -= localCoeffs[k] * pow(obs-(node->getCenter()), k);
            outFile << phi.real() << " ";
        }
        outFile << '\n';
        if (p < P) {
            resetNode();
            binomTable.clear();
        }
    }

    auto iList = node->getInteractionList();
    auto iListBase = node->getBase()->getInteractionList();
    iList.insert(iList.end(), iListBase.begin(), iListBase.end());
    for (const auto& obs : node->getPsn()) {
        cmplx phi;
        for (const auto& iNode : iList) {
            cmplxVec srcs = iNode->getPsn();
            for (size_t src = 0; src < srcs.size(); ++src)
                phi -= (iNode->getQs())[src] * std::log(obs - srcs[src]);
        }
        outAnlFile << phi.real() << " ";
    }
    outAnlFile << '\n';
    
    // Label the chosen node and interaction list
    node->setNodeStat(3);

    for (const auto& iNode : iList)
        iNode->setNodeStat(2);

    std::ofstream psnFile, nodeFile;
    psnFile.open("out/srcs.txt");
    nodeFile.open("out/nodes.txt");

    printPsn(psnFile);
    printNode(nodeFile);
}

//void Stem::iListTest() {
//    std::random_device rd;
//    std::mt19937 gen(rd());
//
//    std::uniform_int_distribution<> branchIdx(0, 3);
//    std::shared_ptr<Node> node = std::make_shared<Stem>(*this);
//
//    while (node->isNodeType<Stem>())
//        // while (node->getLvl() > 1)
//        node = (node->getBranches())[branchIdx(gen)];
//    node->setNodeStat(3);
//
//    node->buildNearNeighbors();
//    node->buildInteractionList();
//    auto iList = node->getInteractionList();
//
//    for (const auto& iNode : iList)
//        iNode->setNodeStat(2);
//
//    std::ofstream psnFile, nodeFile;
//    psnFile.open("out/srcs.txt");
//    nodeFile.open("out/nodes.txt");
//
//    printPsn(psnFile);
//    printNode(nodeFile);
//}