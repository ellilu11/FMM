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

    const int P = Node::getExpansionOrder();
    for (int p = 1; p <= P; ++p) {
        Node::setExpansionOrder(p);
        Node::buildBinomTable();

        cout << " Computing upward pass...   (" << " P = " << p << " )\n";
        buildMpoleCoeffs();

        cout << " Computing downward pass... (" << " P = " << p << " )\n";
        buildLocalCoeffs();

        leafNode->evaluatePhi();
        leafNode->printPhi(outFile);
        outFile << '\n';

        if (p < P) {
            resetNode();
            binomTable.clear();
        }
    }

    //auto iList = node->getInteractionList();
    //cout << iList.size() << '\n';
    //
    //auto iListBase = node->getBase()->getInteractionList();
    //cout << iListBase.size() << '\n';
    //iList.insert(iList.end(), iListBase.begin(), iListBase.end());

    for (const auto& obs : node->getParticles()) {
        cmplx phi;

        // for (size_t src = 0; src < pos.size(); ++src) {
        for (const auto& src : particles) {
            // auto dist = obs - pos[src];
            // if (abs(dist) > 1.0E-9)
            if (src != obs)
                phi -= src->getCharge() * log(obs->getPos() - src->getPos());
        }
        outAnlFile << phi.real() << " ";
    }
    outAnlFile << '\n';

    // Label the chosen node and interaction list
    node->setNodeStat(3);

    for (const auto& iNode : iList)
        iNode->setNodeStat(2);

    std::ofstream nodeFile;
    nodeFile.open("out/nodes.txt");

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
//    std::ofstream posFile, nodeFile;
//    posFile.open("out/srcs.txt");
//    nodeFile.open("out/nodes.txt");
//
//    printPos(posFile);
//    printNode(nodeFile);
//}