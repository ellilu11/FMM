#pragma once

#include "node.h"
#include "stem.h"

class Leaf;

using LeafVec = std::vector<std::shared_ptr<Leaf>>;
using LeafPair = std::pair<std::shared_ptr<Leaf>, std::shared_ptr<Leaf>>;

class Leaf final : public Node, public std::enable_shared_from_this<Leaf> {
public:
    Leaf(
        const ParticleVec&,
        const int,
        Stem* const);

    void evalFarSols();

    void evalNearNonNborSols();

    static std::vector<LeafPair> findNearNborPairs();

    static void evaluateSols();
    
    // template <typename T, typename Func> std::vector<T> evalNearNborSols(Func);
    // void evalNearNborSols();

    static const int getNumLeaves() { return leaves.size(); }

    void pushToNearNonNbors(const std::shared_ptr<Node>& node) {
        nearNonNbors.push_back(node);
    }

    const NodeVec getNearNonNbors() const { return nearNonNbors; }

    const NodeVec getNearNbors() const { return nearNbors; }

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << label << '\n';
    }

    std::shared_ptr<Node> getSelf() override {
        return shared_from_this();
    }

    void buildNeighbors() override;

    void buildLists() override;

    void buildMpoleCoeffs() override;

    void propagateExpCoeffs() override;

    void buildLocalCoeffs() override;

    // test methods
    std::shared_ptr<Node> getRandNode(int);

    const cmplx getPhiFromBranchMpole(const vec3d&, const int);

    void printMpoleCoeffs(std::ofstream&);

    void resetNode();

private:
    static LeafVec leaves;

    NodeVec nearNbors; // list 1
    NodeVec nearNonNbors; // list 3
};
