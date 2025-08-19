#include "../leaf.h"
#include "../stem.h"

using namespace std;

shared_ptr<Node> Stem::getRandNode(int maxLvl) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution branchIdx(0, 3);

    shared_ptr<Node> node = make_shared<Stem>(*this);

    // while (node->getLvl() < maxLvl)
    while (node->isNodeType<Stem>())
        node = (node->getBranches())[branchIdx(gen)];

    return node;
}

shared_ptr<Node> Leaf::getRandNode(int minLvl) {
    return make_shared<Leaf>(*this);
}


void Node::labelNodes() {
    auto node = getRandNode(1);

    node->labelNode(1); // self

    auto leaf = dynamic_pointer_cast<Leaf>(node);

    for (const auto& node : leaf->getNearNbors())     // list 1
        node->labelNode(2);

    for (const auto& node : leaf->getIlist())         // list 2
        node->labelNode(3);

    for (const auto& node : leaf->getNearNonNbors())  // list 3
        if (node->isNodeType<Leaf>())
            node->labelNode(4);
        else
            node->labelNode(5);

    for (const auto& node : leaf->getLeafIlist())     // list 4
        node->labelNode(6);
}