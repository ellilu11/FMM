#include "../leaf.h"

void Leaf::mpoleToLocalTest() {
    buildLocalCoeffs();
    auto iList = getInteractionList();

    std::ofstream outFile, outAnlFile;
    outFile.open("out/local.txt");
    outAnlFile.open("out/localAnl.txt");

    for (const auto& obs : particles) {
        cmplx phi;
        for (size_t k = 0; k < order; ++k)
            phi -= localCoeffs[k] * pow(obs->getPos()-center, k);
        outFile << phi << " ";
    }
    outFile << "\n";

    for (const auto& obs : particles) {
        cmplx phi;
        for (const auto& iNode : iList) {
            auto srcs = iNode->getParticles();
            // for (size_t src = 0; src < srcs.size(); ++src)
            for (const auto& src : srcs)
                phi -= src->getCharge() * std::log( obs->getPos() - src->getPos() );
        }
        outAnlFile << phi << " ";
    }
    outAnlFile << "\n";

    setNodeStat(3);

    for (const auto& iNode : iList)
        iNode->setNodeStat(2);

    std::ofstream nodeFile;
    nodeFile.open("out/nodes.txt");

    printNode(nodeFile);
}

/*
void Leaf::iListTest() {
    setNodeStat(3);
    buildInteractionList();
    auto nbors = getInteractionList();

    for (const auto& nbor : nbors)
        nbor->setNodeStat(2);

    std::ofstream posFile, nodeFile;
    posFile.open("out/srcs.txt");
    nodeFile.open("out/nodes.txt");

    printPos(posFile);
    printNode(nodeFile);
}*/