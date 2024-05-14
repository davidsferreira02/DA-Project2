#ifndef PROJ1_GRAPH_H
#define PROJ1_GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <unordered_set>

template <class T>
class Edge;

#define INF std::numeric_limits<double>::max()

/************************* Vertex  **************************/

template <class T>
class Vertex {
public:
    Vertex(T in);
    bool operator<(Vertex<T> & vertex) const;

    T getInfo() const;
    std::vector<Edge<T> *> getAdj() const;
    bool isVisited() const;
    bool isProcessing() const;
    unsigned int getIndegree() const;
    double getDist() const;
    Edge<T> *getPath() const;
    std::vector<Edge<T> *> getIncoming() const;

    void setInfo(T info);
    void setVisited(bool visited);
    void setProcesssing(bool processing);
    void setIndegree(unsigned int indegree);
    void setDist(double dist);
    void setPath(Edge<T> *path);
    Edge<T> * addEdge(Vertex<T> *dest, double w);
    bool removeEdge(T in);
    void removeOutgoingEdges();
    std::pair<Edge<T>*, double> removeOutgoingEdge();




protected:
    T info;                // info node
    std::vector<Edge<T> *> adj;  // outgoing edges

    // auxiliary fields
    bool visited = false; // used by DFS, BFS, ...
    bool processing = false;
    unsigned int indegree;
    double dist = 0;
    Edge<T> *path = nullptr;

    std::vector<Edge<T> *> incoming; // incoming edges

    void deleteEdge(Edge<T> *edge);
};

/********************** Edge  ****************************/

template <class T>
class Edge {
public:
    Edge(Vertex<T> *orig, Vertex<T> *dest, double w);

    Vertex<T> * getDest() const;
    double getWeight() const;
    bool isSelected() const;
    Vertex<T> * getOrig() const;
    Edge<T> *getReverse() const;
    double getFlow() const;

    void setSelected(bool selected);
    void setReverse(Edge<T> *reverse);
    void setFlow(double flow);

    void setWeight(double w);

protected:
    Vertex<T> * dest; // destination vertex
    double weight; // edge weight, can also be used for capacity

    // auxiliary fields
    bool selected = false;

    // used for bidirectional edges
    Vertex<T> *orig;
    Edge<T> *reverse = nullptr;

    double flow;
};

template <class T>
void Edge<T>::setWeight(double w) {
    weight = w;
}

/********************** Graph  ****************************/

template <class T>
class Graph {
public:
    ~Graph();
    /*
    * Auxiliary function to find a vertex with a given the content.
    */
    Vertex<T> *findVertex(const T &in) const;
    /*
     *  Adds a vertex with a given content or info (in) to a graph (this).
     *  Returns true if successful, and false if a vertex with that content already exists.
     */
    bool addVertex(const T &in);
    bool removeVertex(const T &in);



    /*
     * Adds an edge to a graph (this), given the contents of the source and
     * destination vertices and the edge weight (w).
     * Returns true if successful, and false if the source or destination vertex does not exist.
     */
    bool addEdge(const T &sourc, const T &dest, double w);


    std::vector<Edge<T>*> primMST(int startIdx) {
        std::vector<Edge<T>*> mst;
        if (vertexSet.empty() || startIdx < 0 || startIdx >= vertexSet.size()) return mst;

        std::priority_queue<std::pair<double, Edge<T>*>, std::vector<std::pair<double, Edge<T>*>>, std::greater<>> pq;
        std::vector<bool> inMST(vertexSet.size(), false);

        // Start from the specified vertex
        inMST[startIdx] = true;
        for (auto edge : vertexSet[startIdx]->getAdj()) {
            pq.push({edge->getWeight(), edge});
        }

        while (!pq.empty()) {
            auto [weight, edge] = pq.top();
            pq.pop();
            Vertex<T>* dest = edge->getDest();

            int destIdx = findVertexIdx(dest->getInfo());
            if (destIdx != -1 && !inMST[destIdx]) {
                mst.push_back(edge);
                inMST[destIdx] = true;

                for (auto nextEdge : dest->getAdj()) {
                    if (!inMST[findVertexIdx(nextEdge->getDest()->getInfo())]) {
                        pq.push({nextEdge->getWeight(), nextEdge});
                    }
                }
            }
        }

        return mst;
    }





// To start the traversal


    Edge<T>* findEdge(const T& source, const T& dest) const {
        Vertex<T>* srcVertex = findVertex(source);
        Vertex<T>* destVertex = findVertex(dest);
        if (srcVertex == nullptr || destVertex == nullptr) {
            return nullptr; // Source or destination vertex doesn't exist
        }
        // Search for the edge in the source vertex's adjacency list
        for (auto edge : srcVertex->getAdj()) {
            if (edge->getDest() == destVertex) {
                return edge;
            }
        }
        return nullptr; // Edge not found
    }
    bool removeEdge(const T &source, const T &dest);
    bool addBidirectionalEdge(const T &sourc, const T &dest, double w);

    int getNumVertex() const;
    std::vector<Vertex<T> *> getVertexSet() const;

    std:: vector<T> dfs() const;
    std:: vector<T> dfs(const T & source) const;
    void dfsVisit(Vertex<T> *v,  std::vector<T> & res) const;
    std::vector<T> bfs(const T & source) const;

    bool isDAG() const;
    bool dfsIsDAG(Vertex<T> *v) const;
    std::vector<T> topsort() const;
protected:
    std::vector<Vertex<T> *> vertexSet;    // vertex set

    double ** distMatrix = nullptr;   // dist matrix for Floyd-Warshall
    int **pathMatrix = nullptr;   // path matrix for Floyd-Warshall

    /*
     * Finds the index of the vertex with a given content.
     */
    int findVertexIdx(const T &in) const;
};



void deleteMatrix(int **m, int n);
void deleteMatrix(double **m, int n);



/************************* Vertex  **************************/

template <class T>
Vertex<T>::Vertex(T in): info(in) {}
/*
 * Auxiliary function to add an outgoing edge to a vertex (this),
 * with a given destination vertex (d) and edge weight (w).
 */
template <class T>
Edge<T> * Vertex<T>::addEdge(Vertex<T> *d, double w) {
    auto newEdge = new Edge<T>(this, d, w);
    adj.push_back(newEdge);
    d->incoming.push_back(newEdge);
    return newEdge;
}


/*
 * Auxiliary function to remove an outgoing edge (with a given destination (d))
 * from a vertex (this).
 * Returns true if successful, and false if such edge does not exist.
 */
template <class T>
bool Vertex<T>::removeEdge(T in) {
    bool removedEdge = false;
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        Vertex<T> *dest = edge->getDest();
        if (dest->getInfo() == in) {
            it = adj.erase(it);
            deleteEdge(edge);
            removedEdge = true; // allows for multiple edges to connect the same pair of vertices (multigraph)
        }
        else {
            it++;
        }
    }
    return removedEdge;
}

/*
 * Auxiliary function to remove an outgoing edge of a vertex.
 */
template <class T>
void Vertex<T>::removeOutgoingEdges() {
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        it = adj.erase(it);
        deleteEdge(edge);
    }
}

template <class T>
bool Vertex<T>::operator<(Vertex<T> & vertex) const {
    return this->dist < vertex.dist;
}

template <class T>
T Vertex<T>::getInfo() const {
    return this->info;
}

template <class T>
std::vector<Edge<T>*> Vertex<T>::getAdj() const {
    return this->adj;
}

template <class T>
bool Vertex<T>::isVisited() const {
    return this->visited;
}

template <class T>
bool Vertex<T>::isProcessing() const {
    return this->processing;
}

template <class T>
unsigned int Vertex<T>::getIndegree() const {
    return this->indegree;
}

template <class T>
double Vertex<T>::getDist() const {
    return this->dist;
}

template <class T>
Edge<T> *Vertex<T>::getPath() const {
    return this->path;
}

template <class T>
std::vector<Edge<T> *> Vertex<T>::getIncoming() const {
    return this->incoming;
}

template <class T>
void Vertex<T>::setInfo(T in) {
    this->info = in;
}

template <class T>
void Vertex<T>::setVisited(bool visited) {
    this->visited = visited;
}

template <class T>
void Vertex<T>::setProcesssing(bool processing) {
    this->processing = processing;
}

template <class T>
void Vertex<T>::setIndegree(unsigned int indegree) {
    this->indegree = indegree;
}

template <class T>
void Vertex<T>::setDist(double dist) {
    this->dist = dist;
}

template <class T>
void Vertex<T>::setPath(Edge<T> *path) {
    this->path = path;
}

template <class T>
void Vertex<T>::deleteEdge(Edge<T> *edge) {
    Vertex<T> *dest = edge->getDest();
    // Remove the corresponding edge from the incoming list
    auto it = dest->incoming.begin();
    while (it != dest->incoming.end()) {
        if ((*it)->getOrig()->getInfo() == info) {
            it = dest->incoming.erase(it);
        }
        else {
            it++;
        }
    }
    delete edge;
}
template <class T>
std::pair<Edge<T>*, double> Vertex<T>::removeOutgoingEdge(){
    Edge<T>* edge;
    double weight;
    for (auto it = adj.begin(); it != adj.end(); ++it) {
        if ((*it)->getWeight() != 0) {
            edge = *it;
            weight = edge->getWeight();
            edge->setWeight(0);
            break;
        }
    }

    return std::make_pair(edge, weight);
}



/********************** Edge  ****************************/

template <class T>
Edge<T>::Edge(Vertex<T> *orig, Vertex<T> *dest, double w): orig(orig), dest(dest), weight(w) {}

template <class T>
Vertex<T> * Edge<T>::getDest() const {
    return this->dest;
}

template <class T>
double Edge<T>::getWeight() const {
    return this->weight;
}

template <class T>
Vertex<T> * Edge<T>::getOrig() const {
    return this->orig;
}

template <class T>
Edge<T> *Edge<T>::getReverse() const {
    return this->reverse;
}

template <class T>
bool Edge<T>::isSelected() const {
    return this->selected;
}

template <class T>
double Edge<T>::getFlow() const {
    return flow;
}

template <class T>
void Edge<T>::setSelected(bool selected) {
    this->selected = selected;
}

template <class T>
void Edge<T>::setReverse(Edge<T> *reverse) {
    this->reverse = reverse;
}

template <class T>
void Edge<T>::setFlow(double flow) {
    this->flow = flow;
}

/********************** Graph  ****************************/

template <class T>
int Graph<T>::getNumVertex() const {
    return vertexSet.size();
}

template <class T>
std::vector<Vertex<T> *> Graph<T>::getVertexSet() const {
    return vertexSet;
}

/*
 * Auxiliary function to find a vertex with a given content.
 */
template <class T>
Vertex<T> * Graph<T>::findVertex(const T &in) const {
    for (auto v : vertexSet)
        if (v->getInfo() == in)
            return v;
    return nullptr;
}

/*
 * Finds the index of the vertex with a given content.
 */
template <class T>
int Graph<T>::findVertexIdx(const T &in) const {
    for (unsigned i = 0; i < vertexSet.size(); i++)
        if (vertexSet[i]->getInfo() == in)
            return i;
    return -1;
}
/*
 *  Adds a vertex with a given content or info (in) to a graph (this).
 *  Returns true if successful, and false if a vertex with that content already exists.
 */
template <class T>
bool Graph<T>::addVertex(const T &in) {
    if (findVertex(in) != nullptr)
        return false;
    vertexSet.push_back(new Vertex<T>(in));
    return true;
}

/*
 *  Removes a vertex with a given content (in) from a graph (this), and
 *  all outgoing and incoming edges.
 *  Returns true if successful, and false if such vertex does not exist.
 */
template <class T>
bool Graph<T>::removeVertex(const T &in) {
    for (auto it = vertexSet.begin(); it != vertexSet.end(); it++) {
        if ((*it)->getInfo() == in) {
            auto v = *it;
            v->removeOutgoingEdges();
            for (auto u : vertexSet) {
                u->removeEdge(v->getInfo());
            }
            vertexSet.erase(it);
            delete v;
            return true;
        }
    }
    return false;
}

/*
 * Adds an edge to a graph (this), given the contents of the source and
 * destination vertices and the edge weight (w).
 * Returns true if successful, and false if the source or destination vertex does not exist.
 */
template <class T>
bool Graph<T>::addEdge(const T &sourc, const T &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    v1->addEdge(v2, w);
    return true;
}

/*
 * Removes an edge from a graph (this).
 * The edge is identified by the source (sourc) and destination (dest) contents.
 * Returns true if successful, and false if such edge does not exist.
 */
template <class T>
bool Graph<T>::removeEdge(const T &sourc, const T &dest) {
    Vertex<T> * srcVertex = findVertex(sourc);
    if (srcVertex == nullptr) {
        return false;
    }
    return srcVertex->removeEdge(dest);
}

template <class T>
bool Graph<T>::addBidirectionalEdge(const T &sourc, const T &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    auto e1 = v1->addEdge(v2, w);
    auto e2 = v2->addEdge(v1, w);
    e1->setReverse(e2);
    e2->setReverse(e1);
    return true;
}

/****************** DFS ********************/

/*
 * Performs a depth-first search (dfs) traversal in a graph (this).
 * Returns a vector with the contents of the vertices by dfs order.
 */
template <class T>
std::vector<T> Graph<T>::dfs() const {
    std::vector<T> res;
    for (auto v : vertexSet)
        v->setVisited(false);
    for (auto v : vertexSet)
        if (!v->isVisited())
            dfsVisit(v, res);
    return res;
}

/*
 * Performs a depth-first search (dfs) in a graph (this) from the source node.
 * Returns a vector with the contents of the vertices by dfs order.
 */
template <class T>
std::vector<T> Graph<T>::dfs(const T & source) const {
    std::vector<int> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }
    // Set that no vertex has been visited yet
    for (auto v : vertexSet) {
        v->setVisited(false);
    }
    // Perform the actual DFS using recursion
    dfsVisit(s, res);

    return res;
}

/*
 * Auxiliary function that visits a vertex (v) and its adjacent, recursively.
 * Updates a parameter with the list of visited node contents.
 */
template <class T>
void Graph<T>::dfsVisit(Vertex<T> *v, std::vector<T> & res) const {
    v->setVisited(true);
    res.push_back(v->getInfo());
    for (auto & e : v->getAdj()) {
        auto w = e->getDest();
        if (!w->isVisited()) {
            dfsVisit(w, res);
        }
    }
}

/****************** BFS ********************/
/*
 * Performs a breadth-first search (bfs) in a graph (this), starting
 * from the vertex with the given source contents (source).
 * Returns a vector with the contents of the vertices by bfs order.
 */
template <class T>
std::vector<T> Graph<T>::bfs(const T & source) const {
    std::vector<int> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }

    // Set that no vertex has been visited yet
    for (auto v : vertexSet) {
        v->setVisited(false);
    }

    // Perform the actual BFS using a queue
    std::queue<Vertex<T> *> q;
    q.push(s);
    s->setVisited(true);
    while (!q.empty()) {
        auto v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for (auto & e : v->getAdj()) {
            auto w = e->getDest();
            if ( ! w->isVisited()) {
                q.push(w);
                w->setVisited(true);
            }
        }
    }
    return res;
}

/****************** isDAG  ********************/
/*
 * Performs a depth-first search in a graph (this), to determine if the graph
 * is acyclic (acyclic directed graph or DAG).
 * During the search, a cycle is found if an edge connects to a vertex
 * that is being processed in the stack of recursive calls (see theoretical classes).
 * Returns true if the graph is acyclic, and false otherwise.
 */

template <class T>
bool Graph<T>::isDAG() const {
    for (auto v : vertexSet) {
        v->setVisited(false);
        v->setProcesssing(false);
    }
    for (auto v : vertexSet) {
        if (! v->isVisited()) {
            if ( ! dfsIsDAG(v) ) return false;
        }
    }
    return true;
}


template <class T>
bool Graph<T>::dfsIsDAG(Vertex<T> *v) const {
    v->setVisited(true);
    v->setProcesssing(true);
    for (auto e : v->getAdj()) {
        auto w = e->getDest();
        if (w->isProcessing()) return false;
        if (! w->isVisited()) {
            if (! dfsIsDAG(w)) return false;
        }
    }
    v->setProcesssing(false);
    return true;
}

template<class T>
std::vector<T> Graph<T>::topsort() const {
    std::vector<int> res;

    for (auto v : vertexSet) {
        v->setIndegree(0);
    }
    for (auto v : vertexSet) {
        for (auto e : v->getAdj()) {
            unsigned int indegree = e->getDest()->getIndegree();
            e->getDest()->setIndegree(indegree + 1);
        }
    }

    std::queue<Vertex<T> *> q;
    for (auto v : vertexSet) {
        if (v->getIndegree() == 0) {
            q.push(v);
        }
    }

    while( !q.empty() ) {
        Vertex<T> * v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for(auto e : v->getAdj()) {
            auto w = e->getDest();
            w->setIndegree(w->getIndegree() - 1);
            if(w->getIndegree() == 0) {
                q.push(w);
            }
        }
    }

    if ( res.size() != vertexSet.size() ) {
        //std::cout << "Impossible topological ordering!" << std::endl;
        res.clear();
        return res;
    }

    return res;
}

/**
 * @brief Tests a vertex and visits it if necessary.
 *
 * This function tests a vertex to determine if it needs to be visited
 * based on the residual capacity of an edge. If the vertex has not been
 * visited and the residual capacity is greater than 0, the vertex is
 * marked as visited, its path is set, and it is added to the queue.
 *
 * @tparam T The type of data stored in the graph vertices.
 * @param q A reference to the queue used for BFS traversal.
 * @param e The edge being considered.
 * @param w The vertex being tested and possibly visited.
 * @param residual The residual capacity of the edge.
 */
template <class T>
void testAndVisit(std::queue< Vertex<T>*> &q, Edge<T> *e, Vertex<T> *w, double residual) {

    if (! w->isVisited() && residual > 0) {

        w->setVisited(true);
        w->setPath(e);
        q.push(w);
    }
}

/**
 * @brief Finds an augmenting path in the graph using BFS.
 *
 * This function finds an augmenting path in the graph using BFS traversal.
 * It marks vertices as visited and explores edges to find paths from the
 * source to the target vertex.
 *
 * @tparam T The type of data stored in the graph vertices.
 * @param g A pointer to the graph.
 * @param s The source vertex.
 * @param t The target vertex.
 * @return True if an augmenting path is found, false otherwise.
 */
template <class T>
bool findAugmentingPath(Graph<T> *g, Vertex<T> *s, Vertex<T> *t) {

    for(auto v : g->getVertexSet()) {
        v->setVisited(false);
    }

    s->setVisited(true);
    std::queue<Vertex<T> *> q;
    q.push(s);

    while( ! q.empty() && ! t->isVisited()) {
        auto v = q.front();
        q.pop();

        for(auto e: v->getAdj()) {
            testAndVisit(q, e, e->getDest(), e->getWeight() - e->getFlow());
        }

        for(auto e: v->getIncoming()) {
            testAndVisit(q, e, e->getOrig(), e->getFlow());
        }
    }

    return t->isVisited();
}

/**
 * @brief Finds the minimum residual capacity along an augmenting path.
 *
 * This function finds the minimum residual capacity along an augmenting
 * path from the source to the target vertex.
 *
 * @tparam T The type of data stored in the graph vertices.
 * @param s The source vertex.
 * @param t The target vertex.
 * @return The minimum residual capacity along the augmenting path.
 */
template <class T>
double findMinResidualAlongPath(Vertex<T> *s, Vertex<T> *t) {
    double f = INF;

    for (auto v = t; v != s; ) {
        auto e = v->getPath();
        if (e->getDest() == v) {
            f = std::min(f, e->getWeight() - e->getFlow());
            v = e->getOrig();
        }
        else {
            f = std::min(f, e->getFlow());
            v = e->getDest();
        }
    }

    return f;
}

/**
 * @brief Augments flow along an augmenting path in the graph.
 *
 * This function augments flow along an augmenting path in the graph based
 * on the minimum residual capacity found along the path.
 *
 * @tparam T The type of data stored in the graph vertices.
 * @param s The source vertex.
 * @param t The target vertex.
 * @param f The minimum residual capacity along the augmenting path.
 */
template <class T>
void augmentFlowAlongPath(Vertex<T> *s, Vertex<T> *t, double f) {

    for (auto v = t; v != s; ) {
        auto e = v->getPath();
        double flow = e->getFlow();
        if (e->getDest() == v) {
            e->setFlow(flow + f);
            v = e->getOrig();
        }
        else {
            e->setFlow(flow - f);
            v = e->getDest();
        }
    }
}



inline void deleteMatrix(int **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete [] m[i];
        delete [] m;
    }
}

inline void deleteMatrix(double **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete [] m[i];
        delete [] m;
    }
}

/**
 * @brief Finds the subset of vertices affected by the removal of a given vertex.
 *
 * This function is a theoretical proposal for a graph algorithm that we will use in the presentation. It finds the subset of vertices
 * that are affected by the removal of a specified vertex from the graph. This function is not
 * utilized in the project and is provided for theoretical discussion purposes.
 *
 * @tparam T The type of data held by vertices in the graph.
 * @param g Pointer to the graph structure.
 * @param removedVertex Pointer to the vertex being removed.
 * @return An unordered set of vertices affected by the removal of the specified vertex.
 */
template <class T>
std::unordered_set<Vertex<T>*> findAffectedSubset(Graph<T>* g, Vertex<T>* removedVertex) {
    std::unordered_set<Vertex<T>*> affectedSubset;
    std::queue<Vertex<T>*> q;
    q.push(removedVertex);

    while (!q.empty()) {
        Vertex<T>* v = q.front();
        q.pop();

        affectedSubset.insert(v);
        v->setVisited(true);

        for (auto edge : v->getAdj()) {
            Vertex<T>* w = edge->getDest();
            if (w && !w->isVisited()) {
                q.push(w);
            }
        }
    }

    return affectedSubset;
}


/**
 * @brief Discharges excess flow from a node in the graph.
 *
 * This function is a theoretical proposal for discharging excess flow from a node in the graph.
 * It is part of the relabel-to-front algorithm for calculating maximum flow incrementally.
 * This function is not utilized in the project and is provided for theoretical discussion purposes.
 *
 * @tparam T The type of data held by vertices in the graph.
 * @param graph Reference to the graph structure.
 * @param u Pointer to the node being discharged.
 * @param source Pointer to the source station.
 * @param sink Pointer to the sink station.
 * @param preflow Map containing preflow values for nodes.
 * @param height Map containing height values for nodes.
 * @param excess Map containing excess flow values for nodes.
 * @param activeNodes Queue containing active nodes in the algorithm.
 */
template <class T>
void discharge(Graph<T>& graph, T* u, T* source, T* sink, std::unordered_map<T*, double>& preflow,
               std::unordered_map<T*, double>& height, std::unordered_map<T*, double>& excess,
               std::queue<T*>& activeNodes) {
    while (excess[u] > 0) {
        bool pushed = false;
        for (auto edge : graph.getVertex(u)->getAdj()) {
            T* v = edge->getDest()->getInfo();
            if (height[u] == height[v] + 1 && edge->getFlow() < edge->getWeight()) {
                double delta = std::min(excess[u], edge->getWeight() - edge->getFlow());
                edge->setFlow(edge->getFlow() + delta);
                preflow[u] -= delta;
                preflow[v] += delta;
                excess[u] -= delta;
                excess[v] += delta;
                if (v != sink && v != source && excess[v] > 0) {
                    activeNodes.push(v);
                }
                pushed = true;
                break;
            }
        }
        if (!pushed) {
            relabel(graph, u, source, sink, preflow, height, excess);
        }
    }
}

/**
 * @brief Relabels a node in the graph to maintain the feasibility of the preflow.
 *
 * This function is a theoretical proposal for relabeling a node in the graph to maintain the
 * feasibility of the preflow. It is part of the relabel-to-front algorithm for calculating
 * maximum flow incrementally. This function is not utilized in the project and is provided for
 * theoretical discussion purposes.
 *
 * @tparam T The type of data held by vertices in the graph.
 * @param graph Reference to the graph structure.
 * @param u Pointer to the node being relabeled.
 * @param source Pointer to the source station.
 * @param sink Pointer to the sink station.
 * @param preflow Map containing preflow values for nodes.
 * @param height Map containing height values for nodes.
 * @param excess Map containing excess flow values for nodes.
 */
template <class T>
void relabel(Graph<T>& graph, T* u, T* source, T* sink, std::unordered_map<T*, double>& preflow,
             std::unordered_map<T*, double>& height, std::unordered_map<T*, double>& excess) {
    double minHeight = std::numeric_limits<double>::max();
    for (auto edge : graph.getVertex(u)->getAdj()) {
        T* v = edge->getDest()->getInfo();
        if (edge->getFlow() < edge->getWeight()) {
            minHeight = std::min(minHeight, height[v]);
        }
    }
    height[u] = minHeight + 1;
}


template <class T>
Graph<T>::~Graph() {
    deleteMatrix(distMatrix, vertexSet.size());
    deleteMatrix(pathMatrix, vertexSet.size());
}
#endif //PROJ1_GRAPH_H