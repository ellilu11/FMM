#include <cassert>
#include <iostream>
#include "../leaf.h"
#include "../math.h"
#include "../stem.h"

using namespace std;

void Stem::mpoleToLocalTest() {

    random_device rd;
    mt19937 gen(rd());

    uniform_int_distribution<> branchIdx(0, 3);
    shared_ptr<Node> node = make_shared<Stem>(*this);

    while (node->isNodeType<Stem>())
    // while (node->getLvl() > 1)
        node = (node->getBranches())[branchIdx(gen)];

    auto leafNode = dynamic_pointer_cast<Leaf>(node);

    ofstream outFile, outAnlFile;
    outFile.open("out/local.txt");
    outAnlFile.open("out/localAnl.txt");

    outFile << setprecision(15) << scientific;
    outAnlFile << setprecision(15) << scientific;

    const int P = Node::getP();
    for (int p = 1; p <= P; ++p) {
        Node::setP(p);
        Node::buildBinomTable();

        cout << " Computing upward pass...   (" << " P = " << p << " )\n";
        buildMpoleCoeffs();

        cout << " Computing downward pass... (" << " P = " << p << " )\n";
        buildLocalCoeffs();

        leafNode->evaluatePhi();

        for (const auto& phi : leafNode->getPhi())
            outFile << phi.real() << " ";
        outFile << '\n';

        if (p < P) {
            resetNode();
            binomTable.clear();
        }
    }

    auto iList = node->getInteractionList();
    //cout << iList.size() << '\n';
    //
    //auto iListBase = node->getBase()->getInteractionList();
    //cout << iListBase.size() << '\n';
    //iList.insert(iList.end(), iListBase.begin(), iListBase.end());

    //auto iListBase2 = node->getBase()->getBase()->getInteractionList();
    //cout << iListBase2.size() << '\n';
    //iList.insert(iList.end(), iListBase2.begin(), iListBase2.end());

    //auto iListBase3 = node->getBase()->getBase()->getBase()->getInteractionList();
    //cout << iListBase3.size() << '\n';
    //iList.insert(iList.end(), iListBase3.begin(), iListBase3.end());

    for (const auto& obs : node->getPsn()) {
        cmplx phi;

        for (size_t src = 0; src < psn.size(); ++src) {
            auto dist = obs - psn[src];
            if (abs(dist) > 1.0E-9)
                phi -= qs[src] * std::log(dist);
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