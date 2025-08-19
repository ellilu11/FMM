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

const cmplx Node::getPhiFromMpole(const vec3d& X) {
    auto dR = toSph(X - center);

    double r = dR[0], th = dR[1], ph = dR[2];
    cmplx phi(0, 0);

    for (int l = 0; l <= order; ++l) {
        realVec legendreCosCoeffs;
        for (int m = 0; m <= l; ++m)
            legendreCosCoeffs.push_back(legendreCos(th,l,m));

        for (int m = -l; m <= l; ++m) {
            int m_ = m + l;
            phi += coeffs[l][m_] / pow(r, l+1) *
                legendreCosCoeffs[std::abs(m)] * expI(static_cast<double>(m)*ph)
                // * (m < 0 ? pm(m) : 1.0)
                ;
        }
    }

    return phi;
}

void Node::ffieldTest(const int Nr, const int Nth, const int Nph) {
    std::ofstream obsFile, outFile, outAnlFile, coeffsFile, coeffsRotFile;
    obsFile.open("config/obss.txt");
    outFile.open("out/ff.txt");
    outAnlFile.open("out/ffAnl.txt");
    coeffsFile.open("out/mpolecoeffs.txt");
    coeffsRotFile.open("out/mpolecoeffs_rot.txt");

    // coeffsFile << setprecision(9) << scientific;
    // coeffsRotFile << setprecision(9) << scientific;

    std::vector<vec3d> obss;
    for (int ir = 0; ir < Nr; ++ir){
        double r = 5.0*(ir+1.0)*rootLeng;
        for (int ith = 0; ith < Nth; ++ith) {
            // double th = PI / 2.0; 
            double th = PI * ith / static_cast<double>(Nth);
            for (int iph = 0; iph < Nph; ++iph) {
                double ph = 2.0 * PI * iph / static_cast<double>(Nph);
                auto obs = fromSph(vec3d(r, th, ph));
                obss.push_back(obs);
                obsFile << obs << '\n';
            }
        }
    }

    const int order = Node::getExpansionOrder();
    for (int p = 1; p <= order; ++p) {
        Node::setExpansionOrder(p);

        cout << " Computing upward pass... ( Expansion order = " << p << " )\n";
        auto start = chrono::high_resolution_clock::now();
        
        buildMpoleCoeffs();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration_ms = end - start;
        cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

        for (const auto& obs : obss) {
            auto phi = getPhiFromBranchMpole(obs,0);
            outFile << phi.real() << " ";
        }
        outFile << '\n';
        printMpoleCoeffs(coeffsFile);
        if (p < order) resetNode();
    }

    for (const auto& obs : obss) 
        outAnlFile << getDirectPhi(obs) << " ";
    outAnlFile << '\n';

}

/*
void Node::mpoleToLocalTest() {
    int depth = 0;
    auto node = getRandNode(depth);
    while (!node->getParticles().size())
        node = getRandNode(depth);

    auto iList = node->getInteractionList();

    cout << " This node has " << node->getParticles().size() << " particles";
    cout << " and an interaction list size of " << iList.size() << '\n';

    // auto leafNode = dynamic_pointer_cast<Leaf>(node);

    ofstream outFile, outAnlFile;
    outFile.open("out/local.txt");
    outAnlFile.open("out/localAnl.txt");

    //outFile << setprecision(15) << scientific;
    //outAnlFile << setprecision(15) << scientific;

    const int order = Node::getExpansionOrder();
    for (int p = order; p <= order; ++p) {
        Node::setExpansionOrder(p);

        cout << " Computing upward pass...   (" << " Expansion order = " << p << " )\n";
        buildMpoleCoeffs();

        cout << " Computing downward pass... (" << " Expansion order = " << p << " )\n";
        buildLocalCoeffs();

        for (const auto& obs : node->getParticles()) {
            auto dR = toSph(obs->getPos() - center);

            double r = dR[0], th = dR[1], ph = dR[2];
            cmplx phi(0, 0);

            for (int l = 0; l <= order; ++l) {
                cout << (node->getLocalCoeffs())[l].transpose();
                
                realVec legendreCosCoeffs;
                for (int m = 0; m <= l; ++m)
                    legendreCosCoeffs.push_back(legendreCos(th, pair2i(l, m)));

                for (int m = -l; m <= l; ++m)
                    phi += (node->getLocalCoeffs())[l][m+l] * pow(r, l) *
                    legendreCosCoeffs[std::abs(m)] * expI(static_cast<double>(m)*ph);
            }
            outFile << phi.real() << ' ';
        }
        outFile << '\n';

        if (p < order) resetNode();
    }


    //cout << iList.size() << '\n';
    //
    //auto iListBase = node->getBase()->getInteractionList();
    //cout << iListBase.size() << '\n';
    //iList.insert(iList.end(), iListBase.begin(), iListBase.end());

    for (const auto& obs : node->getParticles()) {
        double phi = 0;
        for (const auto& iNode : iList)
            phi += iNode->getDirectPhi(obs->getPos());

        outAnlFile << phi << ' ';
    }
    outAnlFile << '\n';

}*/

void Node::nfieldTest() {
    using namespace std;

    ofstream outFile, outAnlFile;
    outFile.open("out/nf.txt" 
        // , std::ios::app
    );
    outAnlFile.open("out/nfAnl.txt");

    outFile << setprecision(9) << scientific;
    outAnlFile << setprecision(9) << scientific;

    const int order = Node::getExpansionOrder();
    for (int p = order; p <= order; ++p) {
        Node::setExpansionOrder(p);

        cout << " Computing upward pass...   (" << " Expansion order = " << p << " )\n";
        auto start = chrono::high_resolution_clock::now();

        buildMpoleCoeffs();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration_ms = end - start;
        cout << "   Elapsed time: " << duration_ms.count() << " ms\n";

        cout << " Propagating exponential coeffs...\n";
        start = chrono::high_resolution_clock::now();

        propagateExpCoeffs();

        end = chrono::high_resolution_clock::now();
        duration_ms = end - start;
        cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
        cout << "   Elapsed time (M2X): " << t_M2X.count() << " ms\n";
        cout << "   Elapsed time (X2X): " << t_X2X.count() << " ms\n";

        cout << " Computing downward pass...\n";
        start = chrono::high_resolution_clock::now();

        buildLocalCoeffs();

        end = chrono::high_resolution_clock::now();
        duration_ms = end - start;
        cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
        cout << "   Elapsed time (X2L): " << t_X2L.count() << " ms\n";
        cout << "   Elapsed time (L2L): " << t_L2L.count() << " ms\n";

        for (const auto& src : particles)
            src->printPhi(outFile);
        outFile << '\n';

        if (p < order) resetNode();



    }

    cout << " Computing pairwise..." << endl;
    auto start = chrono::high_resolution_clock::now();

    auto phis = getDirectPhis();
    for (const auto& phi : phis)
        outAnlFile << phi << ' ';
    outAnlFile << '\n';

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration_ms = end - start;
    cout << "   Elapsed time: " << duration_ms.count() << " ms\n";
}