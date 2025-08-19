#pragma once

#include <iostream>
#include <type_traits>
#include "node.h"
#include "stem.h"

class Leaf final : public Node {
public:
    Leaf(
        const ParticleVec&,
        const int,
        Stem* const);

    void getFarNeighbors();
    void buildLists();
    void buildMpoleCoeffs();
    void propagateExpCoeffs();
    void buildLocalCoeffs();

    realVec getFarPhis();
    std::vector<vec3d> getFarFlds();
    // template <typename Func> realVec getNearPhis(Func);
    template <typename T, typename Func> std::vector<T> getNearSols(Func);
    void evaluateSolAtParticles();


    void printPhis(std::ofstream& f) {
        //for (const auto& p : particles)
        //    f << p->getPhi().real() << ' ';
    }

    void printNode(std::ofstream& f) {
        f << center << " " << nodeLeng << " " << nodeStat << '\n';
    }

    // test methods
    std::shared_ptr<Node> getRandNode(int);
    const cmplx getPhiFromBranchMpole(const vec3d&, const int);
    void printMpoleCoeffs(std::ofstream&);
    void resetNode();

private:
    NodeVec farNbors; // list 3
};
