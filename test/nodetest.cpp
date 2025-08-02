#include "../leaf.h"
#include "../stem.h"

using namespace std;

shared_ptr<Node> Stem::getRandNode(int minLvl) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution branchIdx(0, 7);

    shared_ptr<Node> node = make_shared<Stem>(*this);

    // while (minLvl ? node->getLvl() > 1 : node->isNodeType<Stem>() )
    while (node->isNodeType<Stem>())
        node = (node->getBranches())[branchIdx(gen)];

    return node;
}

shared_ptr<Node> Leaf::getRandNode(int minLvl) {
    return make_shared<Leaf>(*this);
}

void Node::setRandNodeStats() {
    auto node = getRandNode(0);
    node->setNodeStat(3);

    node->buildNearNeighbors();
    for (const auto& nbor: node->getNearNeighbors())
        nbor->setNodeStat(1);

    node->buildInteractionList();
    for (const auto& iNode : node->getInteractionList())
        iNode->setNodeStat(2);
}
