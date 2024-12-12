//
//  GrandSubproblem.cpp
//  DefenderAttackerDefender
//
//  Created by Carolina Riascos Alvarez on 2021-03-26.
//  Copyright © 2021 Carolina Riascos Alvarez. All rights reserved.
//

#include "GrandSubproblem.hpp"
void Problem::GrandSubProbMaster(vector<Cycles>&Cycles2ndStage, vector<Chain>&Chains2ndStage, vector<IndexGrandSubSol>&SolFirstStage){
    tStartMP2ndPH = clock();
    // Create model
    GrandSubProb = IloModel(env);
    AtLeastOneFails = IloRangeArray(env);
    
    
    SampleCols2ndStage(Chains2ndStage, Cycles2ndStage, SolFirstStage);
    //Get selected vertices
    vector<int>ListSelVertices = GetSelVertices(SolFirstStage);
    
    //Create arc variables
    mapArcs.clear();
    arc = NumVar2D(env,AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        arc[i] = IloNumVarArray(env, AdjacencyList[i].getSize(), 0, 1, ILOINT);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            mapArcs[make_pair(i, AdjacencyList[i][j] - 1)] = j;
            SetName2(arc[i][j], "arc", i + 1, AdjacencyList[i][j]);
            auto it = ArcsinChCyTHP.find(make_pair(i, AdjacencyList[i][j] - 1));
            if (it == ArcsinChCyTHP.end()){
                arc[i][j].setUB(0);
            }
        }
    }
    //Create vertex variables
    vertex = IloNumVarArray(env, Nodes, 0, 1, ILOINT);
    for (int i = 0; i < Nodes; i++){
        SetName1Index( vertex[i], "vertex",i + 1);
        auto it = NodeCyChTPH.find(i);
        if (it == NodeCyChTPH.end()){
            vertex[i].setUB(0);
        }
    }
    
    //CONSTRAINTS
    string name;
    const char* cName;

    IloExpr sumVertices (env, 0);
    for (int i = 0; i < Nodes; i++){
        sumVertices+=  vertex[i];
    }
    name = "VtxSum";
    cName = name.c_str();
    VxtBudget = IloRange(env, -IloInfinity, sumVertices, MaxVertexFailures, cName);
    GrandSubProb.add(VxtBudget);
    
    IloExpr sumArcs (env, 0);
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            sumArcs+= arc[i][j];
        }
    }
    name = "ArcSum";
    cName = name.c_str();
    ArcBudget = IloRange(env, -IloInfinity, sumArcs, MaxArcFailures, cName);
    GrandSubProb.add(ArcBudget);
    
    //Do not pick a failed arc adjacent to a failed vertex
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        IloExpr in (env,0);  IloExpr out (env,0);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            out+= arc[i][j];
        }
        for (int j = 0; j < PredMap[i].size(); j++){
            in+= arc[PredMap[i][j].first][PredMap[i][j].second];
        }
        GrandSubProb.add(IloRange(env, -IloInfinity, in + out + 2*vertex[i], 2, cName));
        in.end();
        out.end();
    }
    
    //Pairwise revision
    //PairwiseRevision(ListSelVertices);
    
    //First Get at least one vertex/arc fails
    RecoSolCovering.push_back(vector<double>());
    IloNumArray tcysol(env, Cycles2ndStage.size());
    IloNumArray tchsol(env, Chains2ndStage.size());
    double w = 0;
    for (int i = 0; i < Cycles2ndTo3rd.size(); i++){
        tcysol[Cycles2ndTo3rd[i]] = 1;
        w+=Cycles2ndStage[Cycles2ndTo3rd[i]].get_Many();
        RecoSolCovering.back().push_back(Cycles2ndStage[Cycles2ndTo3rd[i]].get_Many());
    }
    for (int i = 0; i < Chains2ndTo3rd.size(); i++){
        tchsol[Chains2ndTo3rd[i]] = 1;
        w+=Chains2ndStage[Chains2ndTo3rd[i]].AccumWeight;
        RecoSolCovering.back().push_back(Chains2ndStage[Chains2ndTo3rd[i]].AccumWeight);
    }
    RecoTotalWCovering.push_back(w);
    sort(RecoSolCovering.back().begin(), RecoSolCovering.back().end(), sortdouble);
    
    if (THP_Method == "Covering"){
        GetAtLeastOneFails(tcysol, tchsol);
    }
    else{
        GetAtLeastOneFailsTwo(tcysol, tchsol);
    }
    
    //cplexGrandSubP.exportModel("GrandSubP.lp");
    //HeuristicsStart2ndPH(Cycles2ndTo3rd, Chains2ndTo3rd,ListSelVertices);
    LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
    cplexGrandSubP = IloCplex(GrandSubProb);
    //cplexGrandSubP.setParam(IloCplex::Param::TimeLimit, LeftTime);
    cplexGrandSubP.setParam(IloCplex::Param::Threads, 1);
    cplexGrandSubP.setOut(env.getNullStream());
    Ite2ndS = 0;
    bool runH = Heuristcs2ndPH();
    if (runH == true){
        runHeuristicstrue++;
        //Build solution
        vertex_sol = IloNumArray(env, Nodes);
        for (int j = 0; j < scenarioHeuristics.size(); j++){
            if (scenarioHeuristics[j].first == - 1){
                vertex_sol[scenarioHeuristics[j].second] = 1;
            }
        }
        arc_sol = IloNumArray2 (env, AdjacencyList.getSize());
        for (int f = 0; f < arc_sol.getSize(); f++){
            arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
        }
        for (int j = 0; j < scenarioHeuristics.size(); j++){
            if (scenarioHeuristics[j].first != - 1){
                for (int i = 0; i < AdjacencyList[scenarioHeuristics[j].first].getSize(); i++){
                    if (AdjacencyList[scenarioHeuristics[j].first][i] == scenarioHeuristics[j].second + 1){
                        arc_sol[scenarioHeuristics[j].first][i] = 1;
                        break;
                    }
                }
            }
        }
        tTotalMP2ndPH += (clock() - tStartMP2ndPH)/double(CLOCKS_PER_SEC);
        THPMIP(Cycles2ndStage, Chains2ndStage, ListSelVertices);
    //}else{
    //    cplexGrandSubP.solve();
    //    if (cplexGrandSubP.getStatus() == IloAlgorithm::Infeasible){
    //        cout << "S.O.S. This should not happen." << endl;
    //    }
    //    else{
    //        vertex_sol = IloNumArray(env, Nodes);
    //        cplexGrandSubP.getValues(vertex_sol,vertex);
    //
    //        arc_sol = IloNumArray2 (env, AdjacencyList.getSize());
    //        for (int f = 0; f < arc_sol.getSize(); f++){
    //            arc_sol[f] = IloNumArray(env, AdjacencyList[f].getSize());
    //            cplexGrandSubP.getValues(arc_sol[f],arc[f]);
    //        }
    //        tTotalMP2ndPH += (clock() - tStartMP2ndPH)/double(CLOCKS_PER_SEC);
    //        THPMIP(Cycles2ndStage, Chains2ndStage, ListSelVertices);
    //    }
    }
    else{
        Print2ndStage("WrongExit");
    }

}
void Problem::PairwiseRevision(vector<int>&ListSelVertices){
    if (RecoursePolicy == "Full"){
        for (int i = 0; i < ListSelVertices.size(); i++){
            for (int j = 0; j < ListSelVertices.size(); j++){
                if (i != j){
                    bool ans = false;
                    vector<int>v;
                    v.push_back(ListSelVertices[j]);
                    vector<pair<int,int>>fv;
                    ans = UnMVtxdueToVtx(v, fv, vector<int>(), make_pair(-1,ListSelVertices[i]));
                    if (ans == false){
                        //cout << vertex[ListSelVertices[i]].getName() << " + " << vertex[ListSelVertices[j]].getName() << endl;
                        GrandSubProb.add(vertex[ListSelVertices[i]] + vertex[ListSelVertices[j]] <= 1);
                    }
                }
            }
        }
    }
}
void Problem::selOnce(map<int,vector<int>>CycleNodeSPH, map<int,vector<int>>ChainNodeSPH){
    
    for (int i = 0; i < Nodes; i++){
        IloExpr expr(env, 0);
        for (int j = 0; j < CycleNodeSPH[i].size(); j++){
            expr+= cyvar[CycleNodeSPH[i][j]];
        }
        for (int j = 0; j < ChainNodeSPH[i].size(); j++){
            expr+= chvar[ChainNodeSPH[i][j]];
        }
        //if (CycleNodeSPH[i].size() > 0 || ChainNodeSPH[i].size() > 0){
            string name = "SelOnce_" + to_string(i + 1);
            //cout << expr << endl;
            Once.add(IloRange(env, -IloInfinity, expr, 1, name.c_str()));
        //}
    }
    GrandSubProb.add(Once);
}
void Problem::IfNoFailuresCY(map<int,int>& Cycles2ndTo3rd, int i){
    string name;

    IloExpr ExprArcs (env, 0);
    IloExpr ExprVtx (env, 0);
    int idx = Cycles2ndTo3rd[i];
    for (int j = 0; j < Cycles2ndStage[idx].get_c().size(); j++){
        ExprVtx+= vertex[Cycles2ndStage[idx].get_c()[j]];
        int u = Cycles2ndStage[idx].get_c()[j];
        if (j == Cycles2ndStage[idx].get_c().size() - 1){
            int v = Cycles2ndStage[idx].get_c()[0];
            int s = mapArcs[make_pair(u,v)];
            ExprArcs+= arc[u][s];
        }
        else{
            int v = Cycles2ndStage[idx].get_c()[j + 1];
            int s = mapArcs[make_pair(u,v)];
            ExprArcs+= arc[u][s];
        }
    }
    //cout << cyvar[i] + ExprVtx + ExprArcs << endl;
    ExprVtx=ExprVtx/(2*MaxVertexFailures);
    ExprArcs=ExprArcs/(2*MaxArcFailures);
    name = "SelIfActive_CY_" + to_string(i + 1);
    //cout << cyvar[i] + ExprVtx + ExprArcs << endl;
    IfNoFailures_CY.add(IloRange(env, -IloInfinity, cyvar[i] + ExprVtx + ExprArcs, 1, name.c_str()));
    GrandSubProb.add(IfNoFailures_CY[IfNoFailures_CY.getSize() - 1]);
    ExprArcs.end();
    ExprVtx.end();
}
void Problem::IfNoFailuresCH(map<int,int>& Chains2ndTo3rd, int i){
    IloExpr ExprArcs (env, 0);
    IloExpr ExprVtx (env, 0);
    int idx = Chains2ndTo3rd[i];
    for (int j = 0; j < Chains2ndStage[idx].Vnodes.size(); j++){
        int u = Chains2ndStage[idx].Vnodes[j].vertex;
        ExprVtx+= vertex[u];
        if (j <= Chains2ndStage[idx].Vnodes.size() - 2){
            int v = Chains2ndStage[idx].Vnodes[j + 1].vertex;
            int s = mapArcs[make_pair(u,v)];
            ExprArcs+= arc[u][s];
        }
    }
    ExprVtx=ExprVtx/(2*MaxVertexFailures);
    ExprArcs=ExprArcs/(2*MaxArcFailures);
    string name = "SelIfActive_CH_" + to_string(i + 1);
    //cout << chvar[i] + ExprVtx + ExprArcs << endl;
    IfNoFailures_CH.add(IloRange(env, -IloInfinity, chvar[i] + ExprVtx + ExprArcs, 1, name.c_str()));
    GrandSubProb.add(IfNoFailures_CH[IfNoFailures_CH.getSize() - 1]);
    ExprArcs.end();
    ExprVtx.end();
}
void Problem::PrintSolSNP(IloNumArray vertex_sol, IloNumArray2 arc_sol){
    cout << "Iteration: " << RepSolCounter << endl;
    for (int i = 0; i < RobustSolTHP.size(); i++){
        cout << "c." << i << " :  ";
        for (int j = 0; j < RobustSolTHP[i].get_c().size(); j++){
            cout << RobustSolTHP[i].get_c()[j] << "\t";
        }
        cout << endl;
    }
    cout << "Vertices failed: ";
    for (int i = 0; i < vertex_sol.getSize(); i++){
        if (vertex_sol[i] > 0.9){
            cout << i << endl;
        }
    }
    cout << "Arcs failed: ";
    for (int i = 0; i < arc_sol.getSize(); i++){
        for (int j = 0; j < arc_sol[i].getSize(); j++){
            if (arc_sol[i][j] > 0.9){
                cout << i << "->" << AdjacencyList[i][j] - 1 << "\t";
            }
        }
    }
    
}
void Problem::UpdateSNPSol(IloNumArray& r_sol, IloNum GrandSubObj){
    RobustObjTPH = GrandSubObj;
    RobustSolTHP.clear();
    for (int i = 0; i < r_sol.getSize(); i++){
        if (r_sol[i] > 0.9){
            //RobustSolTHP.push_back(Cycles(GrandProbSol[i].get_cc(), GrandProbSol[i].get_w()));
        }
    }
}
int FindPosVector(vector<int> array, int value){
    int pos = -1;
    for (int j = 0; j < array.size(); j++){
        if (array[j] == value) return j;
    }
    return pos;
}
vector<Cycles> Problem::Get2ndStageCycles (vector<IndexGrandSubSol>& GrandProbSol, string policy){
    //Create List of all vertices in 1st-stage sol: vinFirstStageSol
    vector<int> vinFirstStageSol;
    if (policy == "Among" || policy == "Full"){
        for (int i = 0; i < GrandProbSol.size(); i++){
            for (int j = 0; j < GrandProbSol[i].get_cc().size(); j++){
                vinFirstStageSol.push_back(GrandProbSol[i].get_cc()[j]);
            }
        }
    }
    
    //Create List of all candidate vertices: ListVertices
    vector<int>ListVertices;
    if (policy == "Among"){
        ListVertices = vinFirstStageSol;
    }
    else if (policy == "Full"){
        for (int i = 0; i < Pairs; i++){
            ListVertices.push_back(i);
        }
    }
    
    vector<Cycles>RecoCycles;
    if (policy == "Among"){
        //Call Find Cycles
        RecoCycles = AmongPolicy(vinFirstStageSol);
    }
    else if (policy == "Full"){
        //Call Find Cycles
        RecoCycles = AllPolicy(vinFirstStageSol);
    }
    else{//BackArcs recourse
        for (int i = 0; i < GrandProbSol.size(); i++){
            //Call Find Cycles
            vector<Cycles>auxCycles;
            vector<int>aux = GrandProbSol[i].get_cc();
            auxCycles = BackRecoursePolicy(aux);
            for (int j = 0; j < auxCycles.size(); j++){
                RecoCycles.push_back(auxCycles[j]);
            }
        }
    }

    
    return RecoCycles;
}
vector<Cycles> Problem::BackRecoursePolicy(vector<int>&ListVertices){
    vector<Cycles> NewListCCs;
    CycleNodeTPH.clear();
    ArcsinCyclesTHP.clear();
    
        IloNumArray2 AdjaList (env, AdjacencyList.getSize());

        for (int j = 0; j < AdjacencyList.getSize(); j++){
            AdjaList[j] = IloNumArray(env);
        }
        for (int j = 0; j < ListVertices.size(); j++){
            for (int l = 0; l < AdjacencyList[ListVertices[j]].getSize(); l++){
                int neighbour = AdjacencyList[ListVertices[j]][l] - 1;
                int pos = FindPosVector(ListVertices, neighbour);
                if (pos != -1) {
                    AdjaList[ListVertices[j]].add(neighbour + 1);
                }
            }
        }
        
        //Find inner cycles
        vector<Cycles> NewList;
        for (int j = 0; j < ListVertices.size(); j++){
            int origin = ListVertices[j];
            NewList = SubCycleFinder(env, AdjaList, origin);
            for (int k = 0; k < NewList.size(); k++){
                NewListCCs.push_back(NewList[k]);
                NewListCCs.back().set_Many(int(NewList[k].get_c().size()));
            }
            //Remove origin from AdjaList
            AdjaList[origin].clear();
        }
        AdjaList.clear();
    
    return NewListCCs;
}
vector<Cycles> Problem::AmongPolicy(vector<int>&ListVertices){
    vector<Cycles> NewListCCs;
    CycleNodeTPH.clear();
    ArcsinCyclesTHP.clear();
    
    IloNumArray2 AdjaList (env, AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        AdjaList[i] = IloNumArray(env);
    }
    
    
    //Fill in Adjacency List
    for (int i = 0; i < ListVertices.size(); i++){
        //cout << endl << LookFor[i] << "\t" << ":";
        for (int j = 0; j < AdjacencyList[ListVertices[i]].getSize(); j++){
            int pos = FindPosVector(ListVertices, AdjacencyList[ListVertices[i]][j] - 1);
            if (pos != -1){
                AdjaList[ListVertices[i]].add(AdjacencyList[ListVertices[i]][j]);
            }
        }
    }
    
    //Find among-cycles
    vector<Cycles> NewList;
    for (int i = 0; i < ListVertices.size(); i++){
        int origin = ListVertices[i];
        NewList = SubCycleFinder(env, AdjaList, origin);
        for (int k = 0; k < NewList.size(); k++){
            NewListCCs.push_back(NewList[k]);
            NewListCCs.back().set_Many(int(NewList[k].get_c().size()));
        }
        //Remove from Adja List
        AdjaList[origin].clear();
    }
    
    return NewListCCs;
}
vector<Cycles> Problem::AllPolicy(vector<int>&ListVertices){
    vector<Cycles> NewListCCs;
    CycleNodeTPH.clear();
    ArcsinCyclesTHP.clear();

    //Make a copy of Adjacency List
    IloNumArray2 AdjaList(env,AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        AdjaList[i] = IloNumArray(env, AdjacencyList[i].getSize());
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            AdjaList[i][j] = AdjacencyList[i][j];
        }
    }
    
    //Find all cycles
    vector<Cycles> NewList;
    for (int i = 0; i < ListVertices.size(); i++){
        int origin = ListVertices[i];
        LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
        if (LeftTime < 0){
            tTotalFindingCyCh+= (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
            Print2ndStage("TimeOut");
        }
        NewList = SubCycleFinder(env, AdjaList, origin);
        for (int k = 0; k < NewList.size(); k++){
            NewListCCs.push_back(NewList[k]);
            int counter = 0;
            for (int l = 0; l < NewList[k].get_c().size(); l++){
                int pos = FindPosVector(ListVertices, NewList[k].get_c()[l]);
                if (pos != -1) counter++;
            }
            NewListCCs.back().set_Many(counter);
        }
        //Remove from Adja List
        AdjaList[origin].clear();
    }
    return NewListCCs;
}
vector<Chain> Problem::Get2ndStageChains (vector<IndexGrandSubSol>& GrandProbSol, string policy){
    
    //Create List of all vertices in 1st-stage sol: vinFirstStageSol
    vector<int> vinFirstStageSol;
    if (policy == "Among" || policy == "Full"){
        for (int i = 0; i < GrandProbSol.size(); i++){
            for (int j = 0; j < GrandProbSol[i].get_cc().size(); j++){
                vinFirstStageSol.push_back(GrandProbSol[i].get_cc()[j]);
            }
        }
    }
    
    //Create List of all candidate vertices: ListVertices
    vector<int>ListVertices;
    if (policy == "Among"){
        ListVertices = vinFirstStageSol;
    }
    else if (policy == "Full"){
        for (int i = 0; i < Nodes; i++){
            ListVertices.push_back(i);
        }
    }
    
    VertexinSolChain = vector<vChain>();
    vector<Chain>RecoChains;
    if (policy == "Full"){
        vector<int>ChainStarters;
        for (int i = Pairs; i < Nodes; i++) ChainStarters.push_back(i);
        //Call InitializeVertexinSolChain
        InitializeVertexinSolChain(ListVertices, VertexinSolChain, AdjacencyList);
        //Call Find Chains
        RecoChains = FindChains(VertexinSolChain, vinFirstStageSol, ChainStarters, false);
    }
    else if (policy == "Among"){
        //Call InitializeVertexinSolChain
        InitializeVertexinSolChain(ListVertices, VertexinSolChain, AdjacencyList);
        //Call Find Chains
        RecoChains = FindChains(VertexinSolChain, vinFirstStageSol, ListVertices, false);
    }
    else{//BackArcs recourse
        for (int i = 0; i < GrandProbSol.size(); i++){
            //Call InitializeVertexinSolChain
            vector<int> aux = GrandProbSol[i].get_cc();
            InitializeVertexinSolChain(aux, VertexinSolChain, AdjacencyList);
            //Call Find Chains
            vector<Chain>auxChains;
            auxChains = FindChains(VertexinSolChain, aux, aux, false);
            for (int j = 0; j < auxChains.size(); j++){
                RecoChains.push_back(auxChains[j]);
            }
        }
    }
    
    
    return RecoChains;
}
void Problem::InitializeVertexinSolChain(vector<int>&ListVertices,vector<vChain>& VertexinSolChain, IloNumArray2 AdjaList){
    vector<int> null(1, -1);
    for (int j = 0; j < ListVertices.size(); j++){
        vector<int> veci;
        for (int l = 0; l < AdjaList[ListVertices[j]].getSize(); l++){
            int neighbour = AdjaList[ListVertices[j]][l] - 1;
            int pos = FindPosVector(ListVertices, neighbour);
            if (pos != -1) {
                veci.push_back(neighbour);
            }
        }
        
        if (veci.size() > 0){
            VertexinSolChain.push_back(vChain(ListVertices[j], veci));
            VertexinSolChain.back().it = VertexinSolChain[ListVertices[j]].veci.begin();
            VertexinSolChain.back().itEnd = VertexinSolChain[ListVertices[j]].veci.end();
        }
        else{
            VertexinSolChain.push_back(vChain(ListVertices[j], null));
            VertexinSolChain.back().it = VertexinSolChain[ListVertices[j]].veci.begin();
            VertexinSolChain.back().itEnd = VertexinSolChain[ListVertices[j]].veci.end();
        }
    }
}
vector<Chain> Problem::FindChains(vector<vChain>& VertexinSolChain, vector<int>& vinFirstStageSol, vector<int>& ListVertices, bool onlyOne){
    vector<vector<int>>Chains;
    vector<int>null(0);
    int u,v;
    bool isin = false;
    vector<vChain>::iterator itVinChain = VertexinSolChain.begin();
    vector<vChain>::iterator itVinChainAux;
    vector<Chain>PPChains;
    bool justPopedBacked = false;
    
    for (itVinChain; itVinChain != VertexinSolChain.end(); itVinChain++){//By construction altruistic donors come first
        if (v2inFirstStageSol(ListVertices, itVinChain->vertex) == true){
            PPChains.push_back(Chain(*itVinChain));
            LeftTime = TimeLimit - (clock() - ProgramStart)/double(CLOCKS_PER_SEC);
            if (LeftTime < 0){
                tTotalFindingCyCh += (clock() - tStart2ndS)/double(CLOCKS_PER_SEC);
                Print2ndStage("TimeOut");
            }
            while(PPChains.back().Vnodes.size() >  0){//next neighbor
                //Increase iterator
                if (PPChains.back().Vnodes.back().FirstPass == true){
                    PPChains.back().Vnodes.back().FirstPass = false;
                }else{PPChains.back().Vnodes.back().it++;}
                u = PPChains.back().Vnodes.back().vertex;
                v = *(PPChains.back().Vnodes.back().it);
                if (PPChains.back().Vnodes.back().veci.size() == 1 && v == -1){
                    PPChains.back().Vnodes.back().it = PPChains.back().Vnodes.back().itEnd;
                }
                //Find pointer to v in VertexinSolChain
                itVinChainAux = VertexinSolChain.begin();
                for (itVinChainAux; itVinChainAux != VertexinSolChain.end(); itVinChainAux++)
                    if (itVinChainAux->vertex == v){
                        break;
                    }
                //If chain is not too long or we have reached the end of some vertex's neighbors
                if (PPChains.back().Vnodes.size() + 1 <= ChainLength  +  1 && PPChains.back().Vnodes.back().it != PPChains.back().Vnodes.back().itEnd){
                    //If vertex not already in chain
                    if (v2AlreadyinChain(PPChains.back().Vnodes, v) == false){
                    //If vertex v is in VertexinSolChain add vertex to current chain, augment AccumWeight, store current chain
                        isin = v2inFirstStageSol(vinFirstStageSol, v);
                        if (isin == false && PPChains.back().Vnodes.size() + 1 == ChainLength  +  1){
                            //Do nothing
                        }
                        else{
                            if (isin == true && PPChains.back().AccumWeight >= 1 && justPopedBacked == false){// Update weight only if v is in the 1st-stage solution
                                Chain aux = PPChains.back();
                                PPChains.push_back(aux);
                                PPChains.back().AccumWeight++;
                            }
                            else {
                                if (isin == true){
                                    PPChains.back().AccumWeight++;
                                }
                            }
                            //Add vertex to current chain
                            PPChains.back().Vnodes.push_back(*itVinChainAux);
                            justPopedBacked = false;
                            if (isin == true && onlyOne == true){
                                break;
                            }
                        }
                    }
                }
                else{
                    if (PPChains.back().AccumWeight > 0){
                        isin = v2inFirstStageSol(vinFirstStageSol, PPChains.back().Vnodes.back().vertex);
                        if (isin == false){
                            PPChains.back().Vnodes.pop_back();
                            Chain aux = PPChains.back();
                            if (justPopedBacked == false) PPChains.push_back(aux);
                        }
                        else{
                            Chain aux = PPChains.back();
                            if (justPopedBacked == false) PPChains.push_back(aux);
                            PPChains.back().Vnodes.pop_back();
                            PPChains.back().AccumWeight--;
                        }
                    }
                    else{
                        PPChains.back().Vnodes.pop_back();
                    }
                    justPopedBacked = true;
                }
            }
            if (onlyOne == true) break;
            PPChains.erase(PPChains.end() - 1);
        }
    }
    
    if (Chains.size() > 0){
        if (Chains[0].size() == 1){
            cout << "S.O.S" << endl;
        }
    }
    return PPChains;
}
bool v2inFirstStageSol(vector<int>sol, int v){
    for (int i = 0; i < sol.size(); i++){
        if (v == sol[i]) return true;
    }
    return false;
}
bool v2AlreadyinChain(vector<vChain> v1, int v2){
    for (int j = 0; j < v1.size(); j++){
        if (v1[j].vertex == v2){
            return true;
        }
    }
    return false;
}
void FindNewNeighbor(vector<Chain>& PPChains){
    while (PPChains.back().Vnodes.back().it != PPChains.back().Vnodes.back().veci.end()){
        PPChains.back().Vnodes.back().it++;
        if (PPChains.back().Vnodes.back().it == PPChains.back().Vnodes.back().itEnd){
            PPChains.back().Vnodes.pop_back();
        }
        else{
            break;
        }
        if (PPChains.back().Vnodes.size() == 0) break;
    }
}
bool sortChains(Chain& c1, Chain& c2){
    return (c1.AccumWeight > c2.AccumWeight);
}
bool sortCycles(Cycles& c1, Cycles& c2){
    return (c1.get_Many() > c2.get_Many());
}
bool IsSameCycle(vector<int> c1, vector<int> c2){
    if (c1.size() == c2.size()){
        for (int j = 0; j < c1.size(); j++){
            if (c1[j] != c2[j]) return false;
        }
    }
    else{
        return false;
    }
    return true;
}
bool IsSameChain(vector<int> c1, vector<vChain> c2){
    if (c1.size() == c2.size()){
        for (int j = 0; j < c1.size(); j++){
            if (c1[j] != c2[j].vertex) return false;
        }
    }
    else{
        return false;
    }
    return true;
}
void Problem::SampleCols2ndStage(vector<Chain>&Chains, vector<Cycles>&Cycles, vector<IndexGrandSubSol>&SolFirstStage){
    //Organize vector by AccumWeight and HowMany
    sort(Chains.begin(), Chains.end(), sortChains);
    sort(Cycles.begin(), Cycles.end(), sortCycles);
    ChainNodeTPH.clear();
    CycleNodeTPH.clear();
    NodeCyChTPH.clear();
    ArcsinChCyTHP.clear();
    ArcsinChainsTHP.clear();
    ArcsinCyclesTHP.clear();
    
//    //Fill in ChainNodeTPH
    for (int i = 0; i < Chains.size(); i++){
        for (int j = 0; j < Chains[i].Vnodes.size(); j++){
            ChainNodeTPH[Chains[i].Vnodes[j].vertex].push_back(i);
            NodeCyChTPH[Chains[i].Vnodes[j].vertex].push_back(i);
            if (j <= Chains[i].Vnodes.size() - 2){
                ArcsinChainsTHP[make_pair(Chains[i].Vnodes[j].vertex, Chains[i].Vnodes[j + 1].vertex)].push_back(i);
                ArcsinChCyTHP[make_pair(Chains[i].Vnodes[j].vertex, Chains[i].Vnodes[j + 1].vertex)].push_back(i);
            }
        }
    }

    //Fill in CycleNodeTPH
    for (int i = 0; i < Cycles.size(); i++){
        for (int j = 0; j < Cycles[i].get_c().size(); j++){
            CycleNodeTPH[Cycles[i].get_c()[j]].push_back(i);
            NodeCyChTPH[Cycles[i].get_c()[j]].push_back(i);
            if (j <= Cycles[i].get_c().size() - 2){
                ArcsinCyclesTHP[make_pair(Cycles[i].get_c()[j], Cycles[i].get_c()[j + 1])].push_back(i);
                ArcsinChCyTHP[make_pair(Cycles[i].get_c()[j], Cycles[i].get_c()[j + 1])].push_back(i);
            }
            else{
                ArcsinCyclesTHP[make_pair(Cycles[i].get_c()[j], Cycles[i].get_c()[0])].push_back(i);
                ArcsinChCyTHP[make_pair(Cycles[i].get_c()[j], Cycles[i].get_c()[0])].push_back(i);
            }
        }
    }
    
    CycleNodeSPH.clear();
    ChainNodeSPH.clear();
    
    
    int counter = 0;
    for (int j = 0; j < SolFirstStage.size(); j++){
        if (SolFirstStage[j].get_cc()[0] < Pairs){// It is a cycle
            for (int i = Cycles.size() - 1; i >= 0; i--){
                //If Cycle[i] is in SolFirstStage
                if (IsSameCycle(SolFirstStage[j].get_cc(), Cycles[i].get_c()) == true){
                    Cycles2ndTo3rd[counter] = i;
                    Cycles3rdTo2nd[i] = counter;
                    for (int j = 0; j < Cycles[i].get_c().size(); j++){
                        CycleNodeSPH[Cycles[i].get_c()[j]].push_back(counter);
                    }
                    counter++;
                    break;//cycle found
                    //if (counter == 3) break;
                }
            }
        }
    }

    counter = 0;
    for (int j = 0; j < SolFirstStage.size(); j++){
        if (SolFirstStage[j].get_cc()[0] >= Pairs){// It is a chain
            for (int i = 0; i < Chains.size(); i++){
                //If Cycle[i] is in SolFirstStage
                if (IsSameChain(SolFirstStage[j].get_cc(), Chains[i].Vnodes) == true){
                    Chains2ndTo3rd[counter] = i;
                    Chains3rdTo2nd[i] = counter;
                    for (int j = 0; j < Chains[i].Vnodes.size(); j++){
                        int u = Chains[i].Vnodes[j].vertex;
                        ChainNodeSPH[u].push_back(counter);
                    }
                    counter++;
                    break;
                    //if (counter == 3) break;
                }
            }
        }
    }
}
vector<int>Problem::GetSelVertices(vector<IndexGrandSubSol>&SolFirstStage){
    vector<int>ListSelVertices;
    for (int i = 0; i < SolFirstStage.size(); i++){
        for (int j = 0; j < SolFirstStage[i].get_cc().size(); j++){
            if (SolFirstStage[i].get_cc()[j] < Pairs){
                ListSelVertices.push_back(SolFirstStage[i].get_cc()[j]);
            }
        }
    }
    
    return ListSelVertices;
}

