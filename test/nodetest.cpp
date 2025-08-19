#include "../node.h"

using namespace std;

void Node::labelNodes() {
    auto node = getRandNode(2);

    node->labelNode(1);

    //for (const auto& node : node->nbors)
    //    node->labelNode(2);

    //for (const auto& iNode : node->getInteractionList())
    //    iNode->setNodeStat(3);

    //for (int dir = 0; dir < 6; ++dir) {
    //    auto iList = (node->dirList)[dir];
    //    for (const auto& iNode : iList) 
    //        iNode->setNodeStat(3+dir);
    //}

    auto leaf = dynamic_pointer_cast<Leaf>(node);

    for (const auto& node : leaf->getFarNbors())
        node->labelNode(3);

    for (int i = 0; i < numDir; ++i)
        for (const auto& node : leaf->getNearNbors(i)) {
            // assert(node->label != 3); // check that node in list 1 is not also in list 3
            node->labelNode(4);
        }
}